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

func Test_frees01a(tst *testing.T) {

	// finalise analysis process and catch errors
	defer End()

	//verbose()
	chk.PrintTitle("frees01a")

	// start simulation
	if !Start("data/frees01.sim", true, chk.Verbose) {
		chk.Panic("cannot start FE simulation")
	}

	// domain
	dom := NewDomain(Global.Sim.Regions[0])
	if dom == nil {
		chk.Panic("cannot run FE simulation")
	}

	// set stage
	if !dom.SetStage(0, Global.Sim.Stages[0]) {
		chk.Panic("cannot set stage\n")
	}

	// nodes and elements
	chk.IntAssert(len(dom.Nodes), 62)
	chk.IntAssert(len(dom.Elems), 15)

	// vertices with "fl"
	seepverts := map[int]bool{3: true, 45: true, 7: true, 49: true, 11: true, 53: true, 15: true, 57: true, 19: true, 61: true, 23: true}

	// check dofs
	var seepeqs []int
	for _, nod := range dom.Nodes {
		if seepverts[nod.Vert.Id] {
			chk.IntAssert(len(nod.Dofs), 2)
			seepeqs = append(seepeqs, nod.Dofs[1].Eq)
		} else {
			chk.IntAssert(len(nod.Dofs), 1)
		}
	}
	sort.Ints(seepeqs)
	io.Pforan("seepeqs = %v\n", seepeqs)
	chk.Ints(tst, "seepeqs", seepeqs, []int{14, 16, 19, 30, 32, 43, 45, 56, 58, 69, 71})

	// check Fmap
	e2 := dom.Elems[2].(*ElemP)
	chk.Ints(tst, "e2.Fmap", e2.Fmap, []int{14, 16, 19})
	e5 := dom.Elems[5].(*ElemP)
	chk.Ints(tst, "e5.Fmap", e5.Fmap, []int{16, 30, 32})
	e8 := dom.Elems[8].(*ElemP)
	chk.Ints(tst, "e8.Fmap", e8.Fmap, []int{30, 43, 45})
	e11 := dom.Elems[11].(*ElemP)
	chk.Ints(tst, "e11.Fmap", e11.Fmap, []int{43, 56, 58})
	e14 := dom.Elems[14].(*ElemP)
	chk.Ints(tst, "e14.Fmap", e14.Fmap, []int{56, 69, 71})
}
