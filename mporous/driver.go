// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/utl"
)

// Driver run simulations with models for porous media
type Driver struct {

	// input
	Mdl *Model // porous model

	// settings
	Silent  bool    // do not show error messages
	CheckD  bool    // do check consistent matrix
	UseDfwd bool    // use DerivFwd (forward differences) instead of DerivCen (central differences) when checking D
	TolCcb  float64 // tolerance to check Ccb
	TolCcd  float64 // tolerance to check Ccd
	VerD    bool    // verbose check of D

	// results
	Res []*State // results
}

// Init initialises driver
func (o *Driver) Init(mdl *Model) (err error) {
	o.Mdl = mdl
	o.TolCcb = 1e-7
	o.TolCcd = 1e-7
	o.VerD = true
	o.CheckD = true
	return
}

// Run runs simulation
func (o *Driver) Run(Pc []float64) (err error) {

	// allocate results arrays
	np := len(Pc)
	o.Res = make([]*State, np)
	o.Res[0] = new(State)

	// initialise first state
	pg, divus := 0.0, 0.0
	err = o.Mdl.InitState(o.Res[0], -Pc[0], pg, divus)
	if err != nil {
		return
	}

	// auxiliary
	derivfcn := num.DerivCen
	if o.UseDfwd {
		derivfcn = num.DerivFwd
	}

	// update states
	var pcOld, pcNew, Δpc, tmp, Ccb, Ccbtmp, Ccd float64
	var stmp State
	for i := 1; i < np; i++ {

		// increment
		pcOld = Pc[i-1]
		pcNew = Pc[i]
		Δpc = pcNew - pcOld

		// update
		o.Res[i] = o.Res[i-1].GetCopy()
		err = o.Mdl.Update(o.Res[i], -Δpc, pg, divus)
		if err != nil {
			return
		}

		// check consistent moduli
		if o.CheckD {

			// check Ccb
			Ccb, err = o.Mdl.Ccb(o.Res[i])
			if err != nil {
				return
			}
			dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, pcNew = pcNew, x
				Δpc = pcNew - pcOld
				stmp.Set(o.Res[i-1])
				o.Mdl.Update(&stmp, -Δpc, pg, divus)
				res, pcNew = stmp.Sl, tmp
				return
			}, pcNew)
			utl.AnaNum(utl.Sf("Ccb @ %.3f,%.4f", pcNew, o.Res[i].Sl), o.TolCcb, Ccb, dnum, o.VerD)

			// check Ccd
			Ccd, err = o.Mdl.Ccd(o.Res[i])
			if err != nil {
				return
			}
			dnum = derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, pcNew = pcNew, x
				Δpc = pcNew - pcOld
				stmp.Set(o.Res[i-1])
				o.Mdl.Update(&stmp, -Δpc, pg, divus)
				Ccbtmp, _ = o.Mdl.Ccb(&stmp)
				res, pcNew = Ccbtmp, tmp
				return
			}, pcNew)
			utl.AnaNum(utl.Sf("Ccd @ %.3f,%.4f", pcNew, o.Res[i].Sl), o.TolCcd, Ccd, dnum, o.VerD)
		}
	}
	return
}
