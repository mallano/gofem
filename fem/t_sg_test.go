// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_sg52a(tst *testing.T) {

	/* Smith & Griffths (5th ed) Figure 5.2 p173
	 *
	 *          0.25       0.5      0.25 kN/m
	 *            ↓         ↓         ↓
	 *    ---    ▷0---------1---------2
	 *     |      |       ,'|       ,'|   E = 1e6 kN/m²
	 *     |      |  0  ,'  |  2  ,'  |   ν = 0.3
	 *     |      |   ,'    |   ,'    |
	 *            | ,'   1  | ,'  3   |   connectivity:
	 *    1 m    ▷3'--------4'--------5     0 : 1 0 3
	 *            |       ,'|       ,'|     1 : 3 4 1
	 *     |      |  4  ,'  |  6  ,'  |     2 : 2 1 4
	 *     |      |   ,'    |   ,'    |     3 : 4 5 2
	 *     |      | ,'   5  | ,'   7  |     4 : 4 3 6
	 *    ---    ▷6'--------7'--------8     5 : 6 7 4
	 *            △         △         △     6 : 5 4 7
	 *                                      7 : 7 8 5
	 *            |------- 1 m -------|
	 */

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sg52a")

	// domain
	Start("data/sg52.sim", true, !utl.Tsilent)
	defer End()
	dom := NewDomain(global.Sim.Regions[0])
	dom.SetStage(0, global.Sim.Stages[0])

	// nodes and elements
	utl.IntAssert(len(dom.Nodes), 9)
	utl.IntAssert(len(dom.Elems), 8)

	// check dofs
	for _, nod := range dom.Nodes {
		utl.IntAssert(len(nod.dofs), 2)
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{1, 0, 3, 4, 2, 5, 6, 7, 8})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17})

	// check solution arrays
	ny := 9 * 2
	nλ := 6
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
		{8, 9, 0, 1, 6, 7},
		{6, 7, 10, 11, 8, 9},
		{6, 7, 4, 5, 12, 13},
		{12, 13, 14, 15, 6, 7},
		{10, 11, 6, 7, 14, 15},
		{14, 15, 16, 17, 10, 11},
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
			tst.Fatalf("key %s is incorrect", c.Key)
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	utl.CompareInts(tst, "constrained ux equations", ct_ux_eqs, []int{2, 4, 12})
	utl.CompareInts(tst, "constrained uy equations", ct_uy_eqs, []int{13, 15, 17})

	// point loads
	utl.IntAssert(len(dom.PtNatBcs.Bcs), 3)
	utl.StrAssert(dom.PtNatBcs.Bcs[0].Key, "fuy")
	utl.StrAssert(dom.PtNatBcs.Bcs[1].Key, "fuy")
	utl.StrAssert(dom.PtNatBcs.Bcs[2].Key, "fuy")
	utl.IntAssert(dom.PtNatBcs.Bcs[0].Eq, 3)
	utl.IntAssert(dom.PtNatBcs.Bcs[1].Eq, 1)
	utl.IntAssert(dom.PtNatBcs.Bcs[2].Eq, 9)
}

func Test_sg52b(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sg52b")

	// run simulation
	Start("data/sg52.sim", true, !utl.Tsilent)
	defer End()
	Run()

	// check
	skipK := false
	tolK := 1e-9
	tolu := 1e-17
	tols := 1.56e-15
	verb := true
	TestingCompareResultsU(tst, "data/sg52.sim", "cmp/sg52.cmp", tolK, tolu, tols, skipK, verb)
}

func Test_sg57(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sg57")

	// run simulation
	Start("data/sg57.sim", true, !utl.Tsilent)
	defer End()
	Run()

	// check
	skipK := false
	tolK := 0.35
	tolu := 2e-9
	tols := 0.0002
	verb := true
	TestingCompareResultsU(tst, "data/sg57.sim", "cmp/sg57.cmp", tolK, tolu, tols, skipK, verb)
}

func Test_sg511(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sg511")

	// run simulation
	Start("data/sg511.sim", true, !utl.Tsilent)
	defer End()
	Run()

	// check
	skipK := false
	tolK := 0.1
	tolu := 3e-14
	tols := 1.56e-7
	verb := true
	TestingCompareResultsU(tst, "data/sg511.sim", "cmp/sg511.cmp", tolK, tolu, tols, skipK, verb)
}

func Test_sg515(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sg515")

	// run simulation
	Start("data/sg515.sim", true, !utl.Tsilent)
	defer End()
	Run()

	// check
	skipK := false
	tolK := 0.15
	tolu := 3e-13
	tols := 3e-8
	verb := true
	TestingCompareResultsU(tst, "data/sg515.sim", "cmp/sg515.cmp", tolK, tolu, tols, skipK, verb)
}

func Test_sg517(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sg517")

	// run simulation
	Start("data/sg517.sim", true, !utl.Tsilent)
	defer End()
	Run()

	// check
	skipK := false
	tolK := 0.0036
	tolu := 1e-6
	tols := 1e-4
	verb := true
	TestingCompareResultsU(tst, "data/sg517.sim", "cmp/sg517.cmp", tolK, tolu, tols, skipK, verb)
}
