// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"time"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

// global data
var global struct {

	// multiprocessing data
	Rank     int   // my rank in distributed cluster
	Nproc    int   // number of processors
	Root     bool  // am I root? (i.e. myrank == 0)
	Distr    bool  // distributed simulation with more than one mpi processor
	Verbose  bool  // verbose == root
	WspcStop []int // stop flags [nprocs]
	WspcInum []int // workspace of integer numbers [nprocs]

	// auxiliar structures
	DynCoefs *DynCoefs // dynamic coefficients

	// simulation and materials
	Sim *inp.Simulation // simulation data
	Mdb *inp.MatDb      // materials database
}

// End must be called and the end to flush log file
func End() {
	inp.FlushLog()
}

// Start initialises 'global' and starts logging
func Start(simfilepath string, erasefiles, verbose bool) {

	// multiprocessing data
	global.Rank = 0
	global.Nproc = 1
	global.Root = true
	global.Distr = false
	if mpi.IsOn() {
		global.Rank = mpi.Rank()
		global.Nproc = mpi.Size()
		global.Root = global.Rank == 0
		global.Distr = global.Nproc > 1
	}
	global.Verbose = verbose
	if !global.Root {
		global.Verbose = false
	}
	if global.Distr {
		global.WspcStop = make([]int, global.Nproc)
		global.WspcInum = make([]int, global.Nproc)
	}

	// auxiliar structures
	global.DynCoefs = new(DynCoefs)

	// simulation and  materials
	dolog := true
	global.Sim = inp.ReadSim(simfilepath, erasefiles, dolog)
	global.Mdb = inp.ReadMat(global.Sim.Data.Matfile, dolog)

	// fix show residual flag
	if !global.Root {
		global.Sim.Data.ShowR = false
	}
}

// Run runs FE simulation
func Run() {

	// alloc domains
	var domains []*Domain
	for _, reg := range global.Sim.Regions {
		domains = append(domains, NewDomain(reg))
	}

	// make sure to call linear solver clean up routines upon exit
	defer func() {
		for _, d := range domains {
			d.LinSol.Clean()
			if global.Verbose {
				utl.Pfgrey(" . . . linear solver cleaned . . .\n")
			}
		}
	}()

	// current time and output time
	t := 0.0
	tout := 0.0

	// message
	if global.Verbose {
		cpu_time := time.Now()
		defer func() {
			utl.Pfblue2("cpu time = %v\n", time.Now().Sub(cpu_time))
		}()
		defer func() {
			utl.Pf("\nfinal time = %g\n", t)
		}()
	}

	// loop over stages
	for stgidx, stg := range global.Sim.Stages {

		// time incrementers
		Dt := stg.Control.DtFunc
		DtOut := stg.Control.DtoFunc
		tf := stg.Control.Tf
		tout = t + DtOut.F(t, nil)
		tidx := 0

		// set stage
		for _, d := range domains {
			d.SetStage(stgidx, global.Sim.Stages[stgidx])
			d.Sol.T = t
			if Stop(d.Out(tidx), "Out") {
				return
			}
		}
		tidx += 1

		// time loop
		var Δt, Δtout float64
		var lasttimestep bool
		for t < tf {

			// time increment
			Δt = Dt.F(t, nil)
			if t+Δt >= tf {
				Δt = tf - t
				lasttimestep = true
			}
			if Δt < global.Sim.Solver.DtMin {
				return
			}

			// time update
			t += Δt
			for _, d := range domains {
				d.Sol.T = t
			}
			Δtout = DtOut.F(t, nil)
			//utl.PfYel(">>>  t=%g  Δt=%g  t+Δt=%g  tf=%g  Δtout=%g\n", t, Δt, t+Δt, tf, Δtout)

			// message
			if global.Verbose {
				//time.Sleep(10000000)
				if !global.Sim.Data.ShowR {
					utl.PrintTimeLong(t)
				}
			}

			// run iterations
			for _, d := range domains {
				stop := run_iterations(t, Δt, d)
				if stop {
					return
				}
			}

			// perform output
			if t >= tout || lasttimestep {
				//utl.Pfblue2("tout = %v\n", tout)
				for _, d := range domains {
					if Stop(d.Out(tidx), "Out") {
						return
					}
				}
				tout += Δtout
				tidx += 1
				//utl.Pforan(">>>  tout=%g  tidx=%d\n", tout, tidx)
			}
		}
	}
	return
}

// run_iterations solves the nonlinear problem
func run_iterations(t, Δt float64, d *Domain) (stop bool) {

	// zero accumulated increments
	la.VecFill(d.Sol.ΔY, 0)

	// calculate global starred vectors and interpolate starred variables from nodes to integration points
	d.star_vars(Δt)

	// auxiliary variables
	var it int
	var largFb, largFb0, Lδu float64

	// message
	if global.Sim.Data.ShowR {
		utl.Pfyel("\n%13s%4s%23s%23s\n", "t", "it", "largFb", "Lδu")
	}
	defer func() {
		if global.Sim.Data.ShowR {
			utl.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}
	}()

	// iterations
	for it = 0; it < global.Sim.Solver.NmaxIt; it++ {

		// assemble right-hand side vector (fb) with negative of residuals
		la.VecFill(d.Fb, 0)
		for _, e := range d.Elems {
			if Stop(e.AddToRhs(d.Fb, d.Sol), "right-hand side vector") {
				return true // stop
			}
		}

		// join all fb
		if global.Distr {
			mpi.AllReduceSum(d.Fb, d.Wb) // this must be done here because there might be nodes sharing boundary conditions
		}

		// point natural boundary conditions; e.g. concentrated loads
		d.PtNatBcs.AddToRhs(d.Fb, t)

		// essential boundary conditioins; e.g. constraints
		d.EssenBcs.AddToRhs(d.Fb, d.Sol)

		// find largest absolute component of fb
		largFb = la.VecLargest(d.Fb, 1)
		if it == 0 {
			// store largest absolute component of fb
			largFb0 = largFb
		} else {
			// check convergence on Lf0
			if largFb < global.Sim.Solver.FbTol*largFb0 { // converged on fb
				break
			}
		}

		// check convergence on fb_min
		if largFb < global.Sim.Solver.FbMin { // converged with smallest value of fb
			break
		}

		// assemble Jacobian matrix
		do_asm_fact := (it == 0 || !global.Sim.Data.CteTg)
		if do_asm_fact {

			// assemble element matrices
			d.Kb.Start()
			for _, e := range d.Elems {
				if Stop(e.AddToKb(d.Kb, d.Sol, it == 0), "assembly") {
					return true // stop
				}
			}

			// join A and tr(A) matrices into Kb
			if global.Root {
				d.Kb.PutMatAndMatT(&d.EssenBcs.A)
			}

			// initialise linear solver
			if d.InitLSol {
				d.LinSol.InitR(d.Kb, global.Sim.LinSol.Symmetric, global.Sim.LinSol.Verbose, global.Sim.LinSol.Timing)
				d.InitLSol = false
			}

			// perform factorisation
			if Stop(d.LinSol.Fact(), "factorisation") {
				return true // stop
			}
		}

		// solve for wb := δyb
		if Stop(d.LinSol.SolveR(d.Wb, d.Fb, false), "solve") {
			return true // stop
		}

		// update primary variables (y) and Lagrange multipliers (λ)
		for i := 0; i < d.Ny; i++ {
			d.Sol.Y[i] += d.Wb[i]  // y += δy
			d.Sol.ΔY[i] += d.Wb[i] // ΔY += δy
		}
		for i := 0; i < d.Nlam; i++ {
			d.Sol.L[i] += d.Wb[d.Ny+i] // λ += δλ
		}

		// backup / restore
		if it == 0 {
			// create backup copy of all secondary variables
			for _, e := range d.ElemIntvars {
				e.BackupIvs()
			}
		} else {
			// recover last converged state from backup copy
			for _, e := range d.ElemIntvars {
				e.RestoreIvs()
			}
		}

		// update secondary variables
		for _, e := range d.Elems {
			if Stop(e.Update(d.Sol), "update") {
				return true // stop
			}
		}

		// compute RMS norm of δu and check convegence on δu
		Lδu = la.VecRmsErr(d.Wb[:d.Ny], global.Sim.Solver.Atol, global.Sim.Solver.Rtol, d.Sol.Y[:d.Ny])

		// message
		if global.Sim.Data.ShowR {
			utl.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}

		// stop if converged on δu
		if Lδu < global.Sim.Solver.Itol {
			break
		}
	}

	// check if iterations diverged
	if it == global.Sim.Solver.NmaxIt {
		utl.PfMag("max number of iterations reached: it = %d\n", it)
		return true // stop
	}

	// success
	return
}
