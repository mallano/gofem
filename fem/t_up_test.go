// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func test_up01a(tst *testing.T) {

	/* this tests simulates seepage flow along a column
	 * by reducing the initial hydrostatic pressure at
	 * at the bottom of the column
	 *
	 *      Nodes / Tags                       Equations
	 *                              ux uy pl               ux uy pl
	 *     8 o----o----o 9 (-5)     53 54 55  o----o----o  50 51 52
	 *       |   14    |             .  .  .  |  58 59  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *    21 o    o    o 22 (-6)    60 61  .  o    o    o  56 57  .
	 *       |   26    |             .  .  .  |  62 63  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     6 o----o----o 7 (-4)     39 40 41  o----o----o  36 37 38
	 *       |   13    |             .  .  .  |  44 45  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *    19 |    o    o 20 (-6)    46 47  .  |    o    o  42 43  .
	 *       |   25    |             .  .  .  |  48 49  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     4 o----o----o 5 (-3)     25 26 27  o----o----o  22 23 24
	 *       |   12    |             .  .  .  |  30 31  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *    17 o    o    o 18 (-6)    32 33  .  o    o    o  28 29  .
	 *       |   24    |             .  .  .  |  34 35  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     2 o----o----o 3 (-2)      9 10 11  o----o----o   6  7  8
	 *       |   11    |             .  .  .  |  16 17  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *    15 o    o    o 16 (-6)    18 19     o    o    o  14 15
	 *       |   23    |             .  .  .  |  20 21  |   .  .  .
	 *       |         |             .  .  .  |         |   .  .  .
	 *     0 o----o----o 1 (-1)      0  1  2  o----o----o   3  4  5
	 *           10                              12 13
	 */

	// capture errors and flush log
	defer End()

	verbose()
	chk.PrintTitle("up01a")

	// start simulation
	if !Start("data/up01.sim", true, chk.Verbose) {
		chk.Panic("cannot start simulation")
	}

	// domain
	dom := NewDomain(Global.Sim.Regions[0])
	if dom == nil {
		chk.Panic("cannot allocate new domain")
	}

	// set stage
	if !dom.SetStage(0, Global.Sim.Stages[0]) {
		chk.Panic("cannot set stage")
	}

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 27)
	chk.IntAssert(len(dom.Elems), 4)

	// nodes with pl
	nods_with_pl := map[int]bool{0: true, 2: true, 4: true, 6: true, 8: true, 1: true, 3: true, 5: true, 7: true, 9: true}

	// check dofs
	for _, nod := range dom.Nodes {
		if nods_with_pl[nod.Vert.Id] {
			chk.IntAssert(len(nod.Dofs), 3)
			chk.StrAssert(nod.Dofs[0].Key, "ux")
			chk.StrAssert(nod.Dofs[1].Key, "uy")
			chk.StrAssert(nod.Dofs[2].Key, "pl")
		} else {
			chk.IntAssert(len(nod.Dofs), 2)
			chk.StrAssert(nod.Dofs[0].Key, "ux")
			chk.StrAssert(nod.Dofs[1].Key, "uy")
		}
	}

	// check equations
	nids, eqs := get_nids_eqs(dom)
	chk.Ints(tst, "eqs", eqs, utl.IntRange(10*3+17*2))
	chk.Ints(tst, "nids", nids, []int{
		0, 1, 3, 2, 10, 16, 11, 15, 23,
		5, 4, 18, 12, 17, 24,
		7, 6, 20, 13, 19, 25,
		9, 8, 22, 14, 21, 26,
	})

	// check pmap
	Pmaps := [][]int{
		{2, 5, 8, 11},
		{11, 8, 24, 27},
		{27, 24, 38, 41},
		{41, 38, 52, 55},
	}
	Umaps := [][]int{
		{0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21},
		{9, 10, 6, 7, 22, 23, 25, 26, 16, 17, 28, 29, 30, 31, 32, 33, 34, 35},
		{25, 26, 22, 23, 36, 37, 39, 40, 30, 31, 42, 43, 44, 45, 46, 47, 48, 49},
		{39, 40, 36, 37, 50, 51, 53, 54, 44, 45, 56, 57, 58, 59, 60, 61, 62, 63},
	}
	for i, ele := range dom.Elems {
		e := ele.(*ElemUP)
		io.Pfpink("%2d : Pmap = %v\n", e.Id(), e.P.Pmap)
		io.Pfpink("%2d : Umap = %v\n", e.Id(), e.U.Umap)
		chk.Ints(tst, "Pmap", e.P.Pmap, Pmaps[i])
		chk.Ints(tst, "Umap", e.U.Umap, Umaps[i])
	}

	// constraints
	chk.IntAssert(len(dom.EssenBcs.Bcs), 9*2+2+3)
	var ct_ux_eqs []int // equations with ux prescribed [sorted]
	var ct_uy_eqs []int // equations with uy prescribed [sorted]
	var ct_pl_eqs []int // equations with pl prescribed [sorted]
	for _, c := range dom.EssenBcs.Bcs {
		chk.IntAssert(len(c.Eqs), 1)
		eq := c.Eqs[0]
		io.Pfgrey("key=%v eq=%v\n", c.Key, eq)
		switch c.Key {
		case "ux":
			ct_ux_eqs = append(ct_ux_eqs, eq)
		case "uy":
			ct_uy_eqs = append(ct_uy_eqs, eq)
		case "pl":
			ct_pl_eqs = append(ct_pl_eqs, eq)
		default:
			tst.Errorf("key %s is incorrect", c.Key)
		}
	}
	sort.Ints(ct_ux_eqs)
	sort.Ints(ct_uy_eqs)
	sort.Ints(ct_pl_eqs)
	chk.Ints(tst, "equations with ux prescribed", ct_ux_eqs, []int{0, 3, 6, 9, 14, 18, 22, 25, 28, 32, 36, 39, 42, 46, 50, 53, 56, 60})
	chk.Ints(tst, "equations with uy prescribed", ct_uy_eqs, []int{1, 4, 13})
	chk.Ints(tst, "equations with pl prescribed", ct_pl_eqs, []int{2, 5})

	// initial values @ nodes
	io.Pforan("initial values @ nodes\n")
	for _, nod := range dom.Nodes {
		z := nod.Vert.C[1]
		eq := nod.Dofs[0].Eq
		chk.Scalar(tst, io.Sf("pl @ %g", z), 1e-17, dom.Sol.Y[eq], 100-10*z)
	}

	// intial values @ integration points
	io.Pforan("initial values @ integration points\n")
	for _, ele := range dom.Elems {
		e := ele.(*ElemP)
		for idx, ip := range e.IpsElem {
			s := e.States[idx]
			z := e.Shp.IpRealCoords(e.X, ip)[1]
			chk.Scalar(tst, io.Sf("sl @ %g", z), 1e-17, s.Sl, 1)
			chk.Scalar(tst, io.Sf("pl @ %g", z), 1e-13, s.Pl, 100-10*z)
		}
	}
}
