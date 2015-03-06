// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_beam01(tst *testing.T) {

	chk.PrintTitle("beam01")

	// domain
	if !Start("data/beam01.sim", true, chk.Verbose) {
		tst.Errorf("test failed\n")
	}
	defer End()
	distr := false
	dom := NewDomain(Global.Sim.Regions[0], distr)
	if dom == nil {
		tst.Errorf("test failed\n")
	}
	if !dom.SetStage(0, Global.Sim.Stages[0], distr) {
		tst.Errorf("test failed\n")
	}

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 2)
	chk.IntAssert(len(dom.Elems), 1)

	// check dofs
	for _, nod := range dom.Nodes {
		chk.IntAssert(len(nod.Dofs), 3)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "nids", nids, []int{0, 1})
	chk.Ints(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5})

	// check solution arrays
	ny := 6
	nλ := 3
	nyb := ny + nλ
	chk.IntAssert(len(dom.Sol.Y), ny)
	chk.IntAssert(len(dom.Sol.Dydt), 0)
	chk.IntAssert(len(dom.Sol.D2ydt2), 0)
	chk.IntAssert(len(dom.Sol.Psi), 0)
	chk.IntAssert(len(dom.Sol.Zet), 0)
	chk.IntAssert(len(dom.Sol.Chi), 0)
	chk.IntAssert(len(dom.Sol.L), nλ)
	chk.IntAssert(len(dom.Sol.ΔY), ny)

	// check linear solver arrays
	chk.IntAssert(len(dom.Fb), nyb)
	chk.IntAssert(len(dom.Wb), nyb)

	// check umap
	e := dom.Elems[0].(*Beam)
	io.Pforan("e = %v\n", e.Umap)
	chk.Ints(tst, "umap", e.Umap, []int{0, 1, 2, 3, 4, 5})

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), nλ)
	var ct_ux_eqs []int // constrained ux equations [sorted]
	var ct_uy_eqs []int // constrained uy equations [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pforan("key=%v eq=%v\n", c.Key, eq)
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
	chk.Ints(tst, "constrained ux equations", ct_ux_eqs, []int{0})
	chk.Ints(tst, "constrained uy equations", ct_uy_eqs, []int{1, 4})
}
