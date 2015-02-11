// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_bh16a(tst *testing.T) {

	/*  solid bracket with thickness = 0.25
	*
	*          1     -10                connectivity:
	*   (-100) o'-,__                    eid :  verts
	*          |     '-,__ 3   -10         0 : 0, 2, 3
	*          |        ,'o-,__            1 : 3, 1, 0
	*          |  1   ,'  |    '-,__ 5     2 : 2, 4, 5
	*          |    ,'    |  3   ,-'o      3 : 5, 3, 2
	*          |  ,'  0   |   ,-'   |
	*          |,'        |,-'   2  |   constraints:
	*   (-100) o----------o---------o    -100 : fixed on x and y
	*          0          2         4
	 */

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("bh16a")

	// domain
	if !Start("data/bh16.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
	}
	defer End()
	dom := NewDomain(global.Sim.Regions[0])
	if dom == nil {
		tst.Errorf("test failed\n")
	}
	if !dom.SetStage(0, global.Sim.Stages[0]) {
		tst.Errorf("test failed\n")
	}

	// nodes and elements
	utl.IntAssert(len(dom.Nodes), 6)
	utl.IntAssert(len(dom.Elems), 4)

	// check dofs
	for _, nod := range dom.Nodes {
		utl.IntAssert(len(nod.dofs), 2)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{0, 2, 3, 1, 4, 5})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11})

	// check solution arrays
	ny := 6 * 2
	nλ := 4
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
	umaps := [][]int{
		{0, 1, 2, 3, 4, 5},
		{4, 5, 6, 7, 0, 1},
		{2, 3, 8, 9, 10, 11},
		{10, 11, 4, 5, 2, 3},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemU)
		utl.Pforan("e%d.umap = %v\n", e.Cell.Id, e.Umap)
		utl.CompareInts(tst, "umap", e.Umap, umaps[i])
	}

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
	utl.CompareInts(tst, "constrained ux equations", ct_ux_eqs, []int{0, 6})
	utl.CompareInts(tst, "constrained uy equations", ct_uy_eqs, []int{1, 7})
}

func Test_bh16b(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("bh16b")

	// run simulation
	if !Start("data/bh16.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
	}
	defer End()
	if !Run() {
		tst.Errorf("test failed\n")
	}

	// check
	skipK := false
	tolK := 1e-12
	tolu := 1e-15
	tols := 1e-12
	verb := true
	TestingCompareResultsU(tst, "data/bh16.sim", "cmp/bh16.cmp", tolK, tolu, tols, skipK, verb)
}

func Test_bh14(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			utl.CallerInfo(3)
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("bh14")

	// run simulation
	if !Start("data/bh14.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
	}
	defer End()
	if !Run() {
		tst.Errorf("test failed\n")
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-15
	tols := 1e-17
	verb := true
	TestingCompareResultsU(tst, "data/bh14.sim", "cmp/bh14.cmp", tolK, tolu, tols, skipK, verb)
}
