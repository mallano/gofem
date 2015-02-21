// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gosl/chk"
)

func Test_sigini01(tst *testing.T) {

	verbose()
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
	d := NewDomain(Global.Sim.Regions[0])
	if !d.SetStage(0, Global.Sim.Stages[0]) {
		tst.Errorf("SetStage failed\n")
		return
	}

	// read results
	sum := ReadSum()
	d.In(sum.NumTidx - 1)

	// solution
	var sol ana.CteStressPstrain
	sol.Init(nil)

	// check displacements
	t := d.Sol.T
	tolu := 1e-17
	for _, n := range d.Nodes {
		eqx := n.GetEq("ux")
		eqy := n.GetEq("uy")
		u := []float64{d.Sol.Y[eqx], d.Sol.Y[eqy]}
		sol.CheckDispl(tst, t, u, n.Vert.C, tolu)
	}

	// check stresses
	e := d.Elems[0].(*ElemU)
	tols := 1e-14
	for idx, ip := range e.IpsElem {
		x := e.Cell.Shp.IpRealCoords(e.X, ip)
		σ := e.States[idx].Sig
		sol.CheckStress(tst, t, σ, x, tols)
	}
}
