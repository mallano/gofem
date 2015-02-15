// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_beam01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("beam01")

	// domain
	if !Start("data/beam01.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
	}
	defer End()
	dom := NewDomain(Global.Sim.Regions[0])
	if dom == nil {
		tst.Errorf("test failed\n")
	}
	if !dom.SetStage(0, Global.Sim.Stages[0]) {
		tst.Errorf("test failed\n")
	}

	// nodes and elements
	utl.IntAssert(len(dom.Nodes), 2)
	utl.IntAssert(len(dom.Elems), 1)

	// check dofs
	for _, nod := range dom.Nodes {
		utl.IntAssert(len(nod.Dofs), 3)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{0, 1})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5})

	// check solution arrays
	ny := 6
	nλ := 3
	nyb := ny + nλ
	utl.IntAssert(len(dom.Sol.Y), ny)
	utl.IntAssert(len(dom.Sol.Dydt), 0)
	utl.IntAssert(len(dom.Sol.D2ydt2), 0)
	utl.IntAssert(len(dom.Sol.Psi), 0)
	utl.IntAssert(len(dom.Sol.Zet), 0)
	utl.IntAssert(len(dom.Sol.Chi), 0)
	utl.IntAssert(len(dom.Sol.L), nλ)
	utl.IntAssert(len(dom.Sol.ΔY), ny)

	// check linear solver arrays
	utl.IntAssert(len(dom.Fb), nyb)
	utl.IntAssert(len(dom.Wb), nyb)

	// check umap
	e := dom.Elems[0].(*Beam)
	utl.Pforan("e = %v\n", e.Umap)
	utl.CompareInts(tst, "umap", e.Umap, []int{0, 1, 2, 3, 4, 5})

	// constraints
	utl.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_ux_eqs []int // constrained ux equations [sorted]
	var ct_uy_eqs []int // constrained uy equations [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		utl.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		utl.Pforan("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "ux":
			ct_ux_eqs = append(ct_ux_eqs, eq)
		case "uy":
			ct_uy_eqs = append(ct_uy_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	utl.CompareInts(tst, "constrained ux equations", ct_ux_eqs, []int{0})
	utl.CompareInts(tst, "constrained uy equations", ct_uy_eqs, []int{1, 4})
}
