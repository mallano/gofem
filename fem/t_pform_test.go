// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_frees01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("frees01")

	// start simulation
	if !Start("data/frees01.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
		return
	}

	// make sure to flush log
	defer End()

	// domain
	dom := NewDomain(global.Sim.Regions[0])
	if dom == nil {
		tst.Errorf("test failed\n")
		return
	}

	// set stage
	if !dom.SetStage(0, global.Sim.Stages[0]) {
		tst.Errorf("test failed\n")
		return
	}

	// nodes and elements
	utl.IntAssert(len(dom.Nodes), 62)
	utl.IntAssert(len(dom.Elems), 15)

	// vertices with "fl"
	seepverts := map[int]bool{3: true, 45: true, 7: true, 49: true, 11: true, 53: true, 15: true, 57: true, 19: true, 61: true, 23: true}

	// check dofs
	var seepeqs []int
	for _, nod := range dom.Nodes {
		if seepverts[nod.vert.Id] {
			utl.IntAssert(len(nod.dofs), 2)
			seepeqs = append(seepeqs, nod.dofs[1].Eq)
		} else {
			utl.IntAssert(len(nod.dofs), 1)
		}
	}
	sort.Ints(seepeqs)
	utl.Pforan("seepeqs = %v\n", seepeqs)
	utl.CompareInts(tst, "seepeqs", seepeqs, []int{14, 16, 19, 30, 32, 43, 45, 56, 58, 69, 71})

	// check Fmap
	e2 := dom.Elems[2].(*ElemP)
	utl.CompareInts(tst, "e2.Fmap", e2.Fmap, []int{14, 16, 19})
	e5 := dom.Elems[5].(*ElemP)
	utl.CompareInts(tst, "e5.Fmap", e5.Fmap, []int{16, 30, 32})
	e8 := dom.Elems[8].(*ElemP)
	utl.CompareInts(tst, "e8.Fmap", e8.Fmap, []int{30, 43, 45})
	e11 := dom.Elems[11].(*ElemP)
	utl.CompareInts(tst, "e11.Fmap", e11.Fmap, []int{43, 56, 58})
	e14 := dom.Elems[14].(*ElemP)
	utl.CompareInts(tst, "e14.Fmap", e14.Fmap, []int{56, 69, 71})
}
