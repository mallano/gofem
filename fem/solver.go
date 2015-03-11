// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"
	"path/filepath"
	"time"

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

	// simulation, materials, meshes and convenience variables
	Sim    *inp.Simulation // simulation data
	Ndim   int             // space dimension
	Dirout string          // directory for output of results
	Fnkey  string          // filename key; e.g. mysim.sim => mysim
	Enc    string          // encoder; e.g. "gob" or "json"
	Stat   bool            // save residuals in summary
	LogBcs bool            // log essential and ptnatural boundary conditions
	Debug  bool            // debug flag

	// auxiliar structures
	DynCoefs *DynCoefs    // dynamic coefficients
	HydroSt  *HydroStatic // computes hydrostatic states

	// for debugging
	DebugKb func(d *Domain, it int) // debug Kb callback function
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

	// simulation and convenience variables
	dir := filepath.Dir(simfilepath)
	fn := filepath.Base(simfilepath)
	Global.Sim = inp.ReadSim(dir, fn, erasefiles)
	LogErrCond(Global.Sim == nil, "ReadSim failed\n")
	if Stop() {
		return
	}
	Global.Ndim = Global.Sim.Ndim
	Global.Dirout = Global.Sim.Data.DirOut
	Global.Fnkey = Global.Sim.Data.FnameKey
	Global.Enc = Global.Sim.Data.Encoder
	Global.Stat = Global.Sim.Data.Stat
	Global.LogBcs = Global.Sim.Data.LogBcs
	Global.Debug = Global.Sim.Data.Debug

	// fix show residual flag
	if !Global.Root {
		Global.Sim.Data.ShowR = false
	}

	// auxiliar structures
	Global.DynCoefs = new(DynCoefs)
	if !Global.DynCoefs.Init(&Global.Sim.Solver) {
		return
	}
	Global.HydroSt = new(HydroStatic)
	Global.HydroSt.Init()

	// success
	return true
}

// Run runs FE simulation
func Run() (runisok bool) {

	// plot functions
	if Global.Sim.PlotF != nil && Global.Root {
		Global.Sim.Functions.PlotAll(Global.Sim.PlotF, Global.Dirout, Global.Fnkey)
	}

	// alloc domains
	var domains []*Domain
	for _, reg := range Global.Sim.Regions {
		dom := NewDomain(reg, Global.Distr)
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

	// summary of outputs; e.g. with output times
	cputime := time.Now()
	var sum Summary
	sum.OutTimes = []float64{t}
	defer func() {
		sum.Save()
		if Global.Verbose && !Global.Debug {
			io.Pf("\nfinal t  = %v\n", t)
			io.Pfblue2("cpu time = %v\n", time.Now().Sub(cputime))
		}
	}()

	// loop over stages
	for stgidx, stg := range Global.Sim.Stages {

		// time incrementers
		Dt := stg.Control.DtFunc
		DtOut := stg.Control.DtoFunc
		tf := stg.Control.Tf
		tout = t + DtOut.F(t, nil)

		// set stage
		for _, d := range domains {
			if LogErrCond(!d.SetStage(stgidx, Global.Sim.Stages[stgidx], Global.Distr), "SetStage failed") {
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

		// skip stage?
		if stg.Skip {
			continue
		}

		// time loop
		ndiverg := 0 // number of steps diverging
		md := 1.0    // time step multiplier if divergence control is on
		var Δt, Δtout float64
		var lasttimestep bool
		for t < tf {

			// check for continued divergence
			if ndiverg >= Global.Sim.Solver.NdvgMax {
				LogErrCond(true, "maximum number of steps diverging reached: %d", ndiverg)
				return
			}

			// time increment
			Δt = Dt.F(t, nil) * md
			if t+Δt >= tf {
				Δt = tf - t
				lasttimestep = true
			}
			if Δt < Global.Sim.Solver.DtMin {
				return true
			}

			// dynamic coefficients
			if LogErr(Global.DynCoefs.CalcBoth(Δt), "cannot compute dynamic coefficients") {
				return
			}

			// time update
			t += Δt
			for _, d := range domains {
				d.Sol.T = t
			}
			Δtout = DtOut.F(t, nil)

			// message
			if Global.Verbose {
				if !Global.Sim.Data.ShowR && !Global.Debug {
					io.PfWhite("time     = %g\r", t)
				}
			}

			// for all domains
			for _, d := range domains {

				// backup solution if divergence control is on
				if Global.Sim.Solver.Diverg {
					d.backup()
				}

				// run iterations
				diverging, ok := run_iterations(t, Δt, d, &sum)
				if !ok {
					return
				}

				// restore solution and reduce time step if divergence control is on
				if Global.Sim.Solver.Diverg {
					if diverging {
						if Global.Verbose {
							io.Pfred(". . . iterations diverging . . .\n")
						}
						d.restore()
						ndiverg += 1
						md *= 0.5
						continue
					}
					ndiverg = 0
					md = 1.0
				}
			}

			// perform output
			if t >= tout || lasttimestep {
				sum.OutTimes = append(sum.OutTimes, t)
				for _, d := range domains {
					//if true {
					if false {
						debug_print_p_results(d)
					}
					if false {
						debug_print_up_results(d)
					}
					if !d.Out(tidx) {
						break
					}
				}
				if Stop() {
					return
				}
				tout += Δtout
				tidx += 1
			}
		}
	}
	return true
}

// run_iterations solves the nonlinear problem
func run_iterations(t, Δt float64, d *Domain, sum *Summary) (diverging, ok bool) {

	// zero accumulated increments
	la.VecFill(d.Sol.ΔY, 0)

	// calculate global starred vectors and interpolate starred variables from nodes to integration points
	if LogErr(d.star_vars(Δt), "cannot compute starred variables") {
		return
	}

	// auxiliary variables
	var it int
	var largFb, largFb0, Lδu float64
	var prevFb, prevLδu float64

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

		// point natural boundary conditions; e.g. concentrated loads
		d.PtNatBcs.AddToRhs(d.Fb, t)

		// essential boundary conditioins; e.g. constraints
		d.EssenBcs.AddToRhs(d.Fb, d.Sol)

		// debug
		if Global.Debug {
			//la.PrintVec("fb", d.Fb[:d.Ny], "%13.10f ", false)
			//panic("stop")
		}

		// find largest absolute component of fb
		largFb = la.VecLargest(d.Fb, 1)

		// save residual
		if Global.Stat {
			sum.Resids.Append(it, largFb)
		}

		// check largFb value
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

		// check divergence on fb
		if it > 1 && Global.Sim.Solver.Diverg {
			if largFb > prevFb {
				diverging = true
				break
			}
		}
		prevFb = largFb

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
				Global.DebugKb(d, it)
			}

			// join A and tr(A) matrices into Kb
			if Global.Root {
				d.Kb.PutMatAndMatT(&d.EssenBcs.A)
			}

			// initialise linear solver
			if d.InitLSol {
				if LogErr(d.LinSol.InitR(d.Kb, Global.Sim.LinSol.Symmetric, Global.Sim.LinSol.Verbose, Global.Sim.LinSol.Timing), "cannot initialise linear solver") {
					return
				}
				d.InitLSol = false
			}

			// perform factorisation
			LogErr(d.LinSol.Fact(), "factorisation")
			if Stop() {
				return
			}
		}

		// debug
		//KK := d.Kb.ToMatrix(nil).ToDense()
		//la.PrintMat("KK", KK, "%20.10f", false)
		//panic("stop")

		// solve for wb := δyb
		LogErr(d.LinSol.SolveR(d.Wb, d.Fb, false), "solve")
		if Stop() {
			return
		}

		// debug
		if Global.Debug {
			//la.PrintVec("wb", d.Wb[:d.Ny], "%13.10f ", false)
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

		// check divergence on Lδu
		if it > 1 && Global.Sim.Solver.Diverg {
			if Lδu > prevLδu {
				diverging = true
				break
			}
		}
		prevLδu = Lδu
	}

	// check if iterations diverged
	if it == Global.Sim.Solver.NmaxIt {
		io.PfMag("max number of iterations reached: it = %d\n", it)
		return
	}

	// success
	ok = true
	return
}

func debug_print_p_results(d *Domain) {
	io.Pf("\ntime = %23.10f\n", d.Sol.T)
	for _, v := range d.Msh.Verts {
		n := d.Vid2node[v.Id]
		eqpl := n.GetEq("pl")
		var pl float64
		if eqpl >= 0 {
			pl = d.Sol.Y[eqpl]
		}
		if math.Abs(pl) < 1e-13 {
			pl = 0
		}
		io.Pf("%3d : pl=%23.10v\n", v.Id, pl)
	}
}

func debug_print_up_results(d *Domain) {
	io.Pf("\ntime = %23.10f\n", d.Sol.T)
	for _, v := range d.Msh.Verts {
		n := d.Vid2node[v.Id]
		eqpl := n.GetEq("pl")
		equx := n.GetEq("ux")
		equy := n.GetEq("uy")
		var pl, ux, uy float64
		if eqpl >= 0 {
			pl = d.Sol.Y[eqpl]
		}
		if equx >= 0 {
			ux = d.Sol.Y[equx]
		}
		if equy >= 0 {
			uy = d.Sol.Y[equy]
		}
		if math.Abs(pl) < 1e-13 {
			pl = 0
		}
		if math.Abs(ux) < 1e-13 {
			ux = 0
		}
		if math.Abs(uy) < 1e-13 {
			uy = 0
		}
		io.Pf("%3d : pl=%23.10v ux=%23.10f uy=%23.10f\n", v.Id, pl, ux, uy)
	}
}
