// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

func Test_sigini01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sigini01")

	// start simulation
	if !Start("data/sigini01.sim", true, chk.Verbose) {
		tst.Errorf("test failed\n")
		return
	}
	defer End()

	// run simulation
	if !Run() {
		tst.Errorf("test failed\n")
		return
	}

	// allocate domain
	distr := false
	d := NewDomain(Global.Sim.Regions[0], distr)
	if !d.SetStage(0, Global.Sim.Stages[0], distr) {
		tst.Errorf("SetStage failed\n")
		return
	}

	// check displacements
	tolu := 1e-16
	for _, n := range d.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{d.Sol.Y[eqx], d.Sol.Y[eqy]}
		chk.Vector(tst, "u", tolu, u, nil)
	}

	// analytical solution
	qnV, qnH := -100.0, -50.0
	ν := 0.25
	σx, σy := qnH, qnV
	σz := ν * (σx + σy)
	σref := []float64{σx, σy, σz, 0}

	// check stresses
	e := d.Elems[0].(*ElemU)
	tols := 1e-13
	for idx, _ := range e.IpsElem {
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		chk.Vector(tst, "σ", tols, σ, σref)
	}
}

func Test_sigini02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("sigini02")

	// start simulation
	if !Start("data/sigini02.sim", true, chk.Verbose) {
		tst.Errorf("test failed\n")
		return
	}
	defer End()

	// run simulation
	if !Run() {
		tst.Errorf("test failed\n")
		return
	}

	// allocate domain
	distr := false
	d := NewDomain(Global.Sim.Regions[0], distr)
	if !d.SetStage(0, Global.Sim.Stages[0], distr) {
		tst.Errorf("SetStage failed\n")
		return
	}

	// read results
	sum := Global.Sum
	sum.Read()
	io.Pforan("sum = %+v\n", sum)
	ntout := len(sum.OutTimes)
	d.In(ntout - 1)

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
		&fun.Prm{N: "qnH0", V: -20},
		&fun.Prm{N: "qnV0", V: -20},
		&fun.Prm{N: "qnH", V: -50},
		&fun.Prm{N: "qnV", V: -100},
	})

	// check displacements
	t := d.Sol.T
	tolu := 1e-16
	for _, n := range d.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{d.Sol.Y[eqx], d.Sol.Y[eqy]}
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	e := d.Elems[0].(*ElemU)
	tols := 1e-13
	for idx, ip := range e.IpsElem {
		x := e.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		io.Pforan("σ = %v\n", σ)
		sol.CheckStress(tst, t, σ, x, tols)
	}
}
