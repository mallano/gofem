// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"time"

	"github.com/cpmech/gosl/io"
)

// Stat handles statistics, timing and summary
var Stat struct {

	// essential
	Sum      *Summary  // the summary
	Detailed bool      // detailed stat is on
	TimeCpu  time.Time // first (cpu) time
	TimeIni  time.Time // initial reference time

	// TODO: can go to summary (per stage?)
	TimeFb  time.Duration // time during computation and assembly of fb
	TimeKb  time.Duration // time during computation and assembly of Kb
	TimeLis time.Duration // time during linear solver call (solution)
	TimeUpd time.Duration // time during updated of element's internal values

	// derived
	cstgidx int // current stage index
	ctidx   int // current output time index
}

func StatEnd(t float64) {
	if Stat.Sum != nil {
		SaveSum(Stat.Sum)
	}
	if Global.Verbose && !Global.Debug {
		io.Pfcyan("\nfinal time = %g\n", t)
		io.Pfblue2("cpu time   = %v\n", time.Now().Sub(Stat.TimeCpu))
	}
}

func StatInit(stgidx, tidx int) {
	if Stat.Sum == nil {
		Stat.Sum = new(Summary)
		Stat.Sum.Times = []float64{0}
		Stat.Sum.StgTidx = make([]int, len(Global.Sim.Stages))
	}
	Stat.Sum.StgTidx[stgidx] = tidx
}

func StatSimtime(t float64) {
	Stat.Sum.Times = append(Stat.Sum.Times, t)
}

// StatInitIter initialises stat for iterations benchmarking
func StatInitIter(stgidx, tidx int) {
	if Stat.Detailed {
		if Stat.Sum.Resids == nil {
		}
		Stat.cstgidx = stgidx
		Stat.ctidx = tidx
	}
}

// StatIniTime resets time
func StatResetTime() {
	Stat.TimeIni = time.Now()
}

// StatFb add info on fb calculation
func StatFb(largFb float64) {
	Stat.TimeFb = time.Now().Sub(Stat.TimeIni)
	if Stat.Detailed {
	}
}

// StatKb adds info on Kb calculation
func StatKb() {
	Stat.TimeKb = time.Now().Sub(Stat.TimeIni)
}

// StatLinSol adds info on linear solution
func StatLinsol() {
	Stat.TimeLis = time.Now().Sub(Stat.TimeIni)
}

// StatUpdate adds info on update
func StatUpdate() {
	Stat.TimeUpd = time.Now().Sub(Stat.TimeIni)
}
