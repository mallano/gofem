// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
)

// HydroStatic computes water pressure (pl) and intrinsic liquid density (ρL)
// based on the following model
//
//    ρL = ρL0 + Cl・pl   thus   dρL/dpl = Cl
//
//    Z(z) = zmax + T・(z - zmax)   with 0 ≤ T ≤ 1
//    dZ   = (z - zmax)・dT
//    dpl  = ρL(pl)・g・(-dZ)
//    dpl  = ρL(pl)・g・(zmax - z)・dT
//    Δz   = zmax - z
//
//            / dpl/dT \   / ρL(pl)・g・Δz \
//    dY/dT = |        | = |               |
//            \ dρL/dT /   \ Cl・dpl/dT    /
//
type HydroStatic struct {
	zwater float64
	ρL0    float64
	g      float64
	Cl     float64
	fcn    ode.Cb_fcn
	Jac    ode.Cb_jac
	sol    ode.ODE
}

func (o *HydroStatic) Init() {

	// basic data
	o.zwater = Global.Sim.WaterLevel
	o.ρL0 = Global.Sim.WaterRho0
	o.g = Global.Sim.Gfcn.F(0, nil)
	o.Cl = o.ρL0 / Global.Sim.WaterBulk

	// x := {pl, ρl}
	o.fcn = func(f []float64, x float64, y []float64, args ...interface{}) error {
		Δz := args[0].(float64)
		//ρL := o.ρL0
		ρL := y[1]
		f[0] = ρL * o.g * Δz // dpl/dT
		f[1] = o.Cl * f[0]   // dρL/dT
		return nil
	}

	o.Jac = func(dfdy *la.Triplet, x float64, y []float64, args ...interface{}) error {
		if dfdy.Max() == 0 {
			dfdy.Init(2, 2, 4)
		}
		Δz := args[0].(float64)
		dfdy.Start()
		dfdy.Put(0, 0, 0)
		dfdy.Put(0, 1, o.g*Δz)
		dfdy.Put(1, 0, 0)
		dfdy.Put(1, 1, o.Cl*o.g*Δz)
		return nil
	}

	silent := true
	o.sol.Init("Radau5", 2, o.fcn, o.Jac, nil, nil, silent)
}

func (o HydroStatic) Calc(z float64) (pl, ρL float64, err error) {
	Δz := o.zwater - z
	y := []float64{0, o.ρL0} // pl0, ρL0
	err = o.sol.Solve(y, 0, 1, 1, false, Δz)
	if err != nil {
		err = chk.Err("HydroStatic failed when calculating pressure using ODE solver: %v", err)
		return
	}
	return y[0], y[1], nil
}

// SetHydroSt sets the initial state to a hydrostatic condition
func (o *Domain) SetHydroSt(stg *inp.Stage) (ok bool) {

	// check for hydrost data
	hst := stg.HydroSt
	if hst == nil {
		return true
	}

	// set Sol
	ndim := Global.Ndim
	for _, n := range o.Nodes {
		z := n.Vert.C[ndim-1]
		dof := n.GetDof("pl")
		if dof != nil {
			pl, _, err := Global.HydroSt.Calc(z)
			if LogErr(err, "hydrost: cannot compute pl") {
				return
			}
			o.Sol.Y[dof.Eq] = pl
		}
	}

	// set elements
	var err error
	for _, e := range o.ElemIntvars {

		// get element's integration points data
		ele := e.(Elem)
		d := ele.OutIpsData()
		nip := len(d)

		// build map with pressures @ ips
		pl := make([]float64, nip)
		ρL := make([]float64, nip)
		for i, ip := range d {
			z := ip.X[ndim-1]
			pl[i], ρL[i], err = Global.HydroSt.Calc(z)
			if LogErr(err, "hydrost: cannot compute pl and ρL") {
				return
			}
		}
		ivs := map[string][]float64{"pl": pl, "ρL": ρL}

		// set element's states
		if LogErrCond(!e.SetIvs(ivs), "hydrostatic: element's internal values setting failed") {
			return
		}
	}

	// success
	log.Printf("dom: initial hydrostatic state set")
	return true
}
