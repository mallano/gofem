// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"path/filepath"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gofem/msolid"

	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/mpi"
)

// Global holds global data
var Global struct {

	// multiprocessing data
	Rank     int   // my rank in distributed cluster
	Nproc    int   // number of processors
	Root     bool  // am I root? (i.e. myrank == 0)
	Distr    bool  // distributed simulation with more than one mpi processor
	Verbose  bool  // verbose == root
	WspcStop []int // stop flags [nprocs]
	WspcInum []int // workspace of integer numbers [nprocs]

	// simulation and materials
	Sim *inp.Simulation // simulation data
	Mdb *inp.MatDb      // materials database

	// auxiliar structures
	DynCoefs *DynCoefs // dynamic coefficients

	// for debugging
	DebugKb func(d *Domain, firstIt bool)
}

// End must be called and the end to flush log file
func End() {
	inp.FlushLog()
}

// Start initialises 'global' and starts logging
func Start(simfilepath string, erasefiles, verbose bool) (startisok bool) {

	// multiprocessing data
	Global.Rank = 0
	Global.Nproc = 1
	Global.Root = true
	Global.Distr = false
	if mpi.IsOn() {
		Global.Rank = mpi.Rank()
		Global.Nproc = mpi.Size()
		Global.Root = Global.Rank == 0
		Global.Distr = Global.Nproc > 1
	}
	Global.Verbose = verbose
	if !Global.Root {
		Global.Verbose = false
	}
	Global.WspcStop = make([]int, Global.Nproc)
	Global.WspcInum = make([]int, Global.Nproc)

	// simulation and  materials
	dir := filepath.Dir(simfilepath)
	fn := filepath.Base(simfilepath)
	Global.Sim = inp.ReadSim(dir, fn, erasefiles)
	LogErrCond(Global.Sim == nil, "ReadSim failed\n")
	if Stop() {
		return
	}
	Global.Mdb = inp.ReadMat(Global.Sim.Data.FnameDir, Global.Sim.Data.Matfile)
	LogErrCond(Global.Mdb == nil, "ReadMat failed\n")
	if Stop() {
		return
	}

	// fix show residual flag
	if !Global.Root {
		Global.Sim.Data.ShowR = false
	}

	// auxiliar structures
	Global.DynCoefs = new(DynCoefs)
	if !Global.DynCoefs.Init(&Global.Sim.Solver) {
		return
	}

	// success
	return true
}

// Run runs FE simulation
func Run() (runisok bool) {

	// alloc domains
	var domains []*Domain
	for _, reg := range Global.Sim.Regions {
		dom := NewDomain(reg)
		if dom == nil {
			break
		}
		domains = append(domains, dom)
	}
	if Stop() {
		return
	}

	// make sure to call linear solver clean up routines upon exit
	defer func() {
		for _, d := range domains {
			if !d.InitLSol {
				d.LinSol.Clean()
			}
		}
	}()

	// current time and output time
	t := 0.0
	tout := 0.0
	tidx := 0

	// message
	if Global.Verbose {
		//cpu_time := time.Now()
		defer func() {
			//io.Pfblue2("cpu time   = %v\n", time.Now().Sub(cpu_time))
		}()
		defer func() {
			//io.Pfcyan("\nfinal time = %g\n", t)
		}()
	}

	// summary
	nstg := len(Global.Sim.Stages)
	var sum Summary
	sum.TidxIni = make([]int, nstg)
	defer func() {
		sum.NumTidx = tidx
		SaveSum(&sum)
	}()

	// loop over stages
	for stgidx, stg := range Global.Sim.Stages {

		// summary
		sum.TidxIni[stgidx] = tidx

		// time incrementers
		Dt := stg.Control.DtFunc
		DtOut := stg.Control.DtoFunc
		tf := stg.Control.Tf
		tout = t + DtOut.F(t, nil)

		// set stage
		for _, d := range domains {
			if !d.SetStage(stgidx, Global.Sim.Stages[stgidx]) {
				break
			}
			d.Sol.T = t
			if !d.Out(tidx) {
				break
			}
		}
		if Stop() {
			return
		}
		tidx += 1

		// log models
		mconduct.LogModels()
		mreten.LogModels()
		mporous.LogModels()
		msolid.LogModels()

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
			if Δt < Global.Sim.Solver.DtMin {
				return true
			}

			// dynamic coefficients
			Global.DynCoefs.CalcBoth(Δt)

			// time update
			t += Δt
			for _, d := range domains {
				d.Sol.T = t
			}
			Δtout = DtOut.F(t, nil)
			//io.PfYel(">>>  t=%g  Δt=%g  t+Δt=%g  tf=%g  Δtout=%g\n", t, Δt, t+Δt, tf, Δtout)

			// message
			if Global.Verbose {
				//time.Sleep(100000000)
				if !Global.Sim.Data.ShowR {
					//io.Pf("time       = %g\r", t)
				}
			}

			// run iterations
			for _, d := range domains {
				if !run_iterations(t, Δt, d) {
					return
				}
			}

			// perform output
			if t >= tout || lasttimestep {
				//io.Pfblue2("tout = %v\n", tout)
				for _, d := range domains {
					if !d.Out(tidx) {
						break
					}
				}
				if Stop() {
					return
				}
				tout += Δtout
				tidx += 1
				//io.Pforan(">>>  tout=%g  tidx=%d\n", tout, tidx)
			}
		}
	}
	return true
}

// run_iterations solves the nonlinear problem
func run_iterations(t, Δt float64, d *Domain) (ok bool) {

	// zero accumulated increments
	la.VecFill(d.Sol.ΔY, 0)

	// calculate global starred vectors and interpolate starred variables from nodes to integration points
	d.star_vars(Δt)

	// auxiliary variables
	var it int
	var largFb, largFb0, Lδu float64

	// message
	if Global.Sim.Data.ShowR {
		io.Pfyel("\n%13s%4s%23s%23s\n", "t", "it", "largFb", "Lδu")
	}
	defer func() {
		if Global.Sim.Data.ShowR {
			io.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}
	}()

	// iterations
	for it = 0; it < Global.Sim.Solver.NmaxIt; it++ {

		// assemble right-hand side vector (fb) with negative of residuals
		la.VecFill(d.Fb, 0)
		for _, e := range d.Elems {
			if !e.AddToRhs(d.Fb, d.Sol) {
				break
			}
		}
		if Stop() {
			return
		}

		// join all fb
		if Global.Distr {
			mpi.AllReduceSum(d.Fb, d.Wb) // this must be done here because there might be nodes sharing boundary conditions
		}

		la.PrintVec("fb", d.Fb, "%13.10f ", false)
		panic("stop")

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
			if largFb < Global.Sim.Solver.FbTol*largFb0 { // converged on fb
				break
			}
		}

		// check convergence on fb_min
		if largFb < Global.Sim.Solver.FbMin { // converged with smallest value of fb
			break
		}

		// assemble Jacobian matrix
		do_asm_fact := (it == 0 || !Global.Sim.Data.CteTg)
		if do_asm_fact {

			// assemble element matrices
			d.Kb.Start()
			for _, e := range d.Elems {
				if !e.AddToKb(d.Kb, d.Sol, it == 0) {
					break
				}
			}
			if Stop() {
				return
			}

			// debug
			if Global.DebugKb != nil {
				Global.DebugKb(d, it == 0)
			}

			// join A and tr(A) matrices into Kb
			if Global.Root {
				d.Kb.PutMatAndMatT(&d.EssenBcs.A)
			}

			// initialise linear solver
			if d.InitLSol {
				d.LinSol.InitR(d.Kb, Global.Sim.LinSol.Symmetric, Global.Sim.LinSol.Verbose, Global.Sim.LinSol.Timing)
				d.InitLSol = false
			}

			// perform factorisation
			LogErr(d.LinSol.Fact(), "factorisation")
			if Stop() {
				return
			}
		}

		// solve for wb := δyb
		LogErr(d.LinSol.SolveR(d.Wb, d.Fb, false), "solve")
		if Stop() {
			return
		}

		// update primary variables (y)
		for i := 0; i < d.Ny; i++ {
			d.Sol.Y[i] += d.Wb[i]  // y += δy
			d.Sol.ΔY[i] += d.Wb[i] // ΔY += δy
		}
		if !Global.Sim.Data.Steady {
			for _, I := range d.T1eqs {
				d.Sol.Dydt[I] = Global.DynCoefs.β1*d.Sol.Y[I] - d.Sol.Psi[I]
			}
			for _, I := range d.T2eqs {
				d.Sol.Dydt[I] = Global.DynCoefs.α4*d.Sol.Y[I] - d.Sol.Chi[I]
				d.Sol.D2ydt2[I] = Global.DynCoefs.α1*d.Sol.Y[I] - d.Sol.Zet[I]
			}
		}

		// update Lagrange multipliers (λ)
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
			if !e.Update(d.Sol) {
				break
			}
		}
		if Stop() {
			return
		}

		// compute RMS norm of δu and check convegence on δu
		Lδu = la.VecRmsErr(d.Wb[:d.Ny], Global.Sim.Solver.Atol, Global.Sim.Solver.Rtol, d.Sol.Y[:d.Ny])

		// message
		if Global.Sim.Data.ShowR {
			io.Pf("%13.6e%4d%23.15e%23.15e\n", t, it, largFb, Lδu)
		}

		// stop if converged on δu
		if Lδu < Global.Sim.Solver.Itol {
			break
		}
	}

	// check if iterations diverged
	if it == Global.Sim.Solver.NmaxIt {
		io.PfMag("max number of iterations reached: it = %d\n", it)
		return
	}

	// success
	return true
}
