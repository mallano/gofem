// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func get_nids_eqs(dom *Domain) (nids, eqs []int) {
	for _, nod := range dom.Nodes {
		nids = append(nids, nod.vert.Id)
		for _, dof := range nod.dofs {
			eqs = append(eqs, dof.Eq)
		}
	}
	return
}

func Test_fourlayers01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("fourlayers01")

	if !Start("data/fourlayers.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
	}
	defer End()
	dom := NewDomain(global.Sim.Regions[0])
	if dom == nil {
		tst.Errorf("test failed\n")
	}

	utl.Pforan("stage # 0\n")
	if !dom.SetStage(0, global.Sim.Stages[0]) {
		tst.Errorf("test failed\n")
	}
	nids, eqs := get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{1, 2, 14, 12, 0, 10})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 12, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11})

	utl.Pforan("stage # 1\n")
	if !dom.SetStage(1, global.Sim.Stages[1]) {
		tst.Errorf("test failed\n")
	}
	nids, eqs = get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{10, 12, 9, 6, 1, 2, 14, 0, 8})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 2, 3, 19, 4, 5, 20, 6, 7, 8, 9, 18, 10, 11, 12, 13, 14, 15, 16, 17})

	utl.Pforan("stage # 2\n")
	if !dom.SetStage(2, global.Sim.Stages[2]) {
		tst.Errorf("test failed\n")
	}
	nids, eqs = get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{10, 12, 9, 6, 1, 2, 14, 0, 7, 11, 8, 13})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 2, 3, 25, 4, 5, 26, 6, 7, 8, 9, 24, 10, 11, 12, 13, 14, 15, 16, 17, 27, 18, 19, 20, 21, 22, 23})

	utl.Pforan("stage # 3\n")
	if !dom.SetStage(3, global.Sim.Stages[3]) {
		tst.Errorf("test failed\n")
	}
	nids, eqs = get_nids_eqs(dom)
	utl.CompareInts(tst, "nids", nids, []int{7, 13, 5, 4, 10, 12, 9, 6, 1, 2, 14, 11, 3, 0, 8})
	utl.CompareInts(tst, "eqs", eqs, []int{0, 1, 33, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 31, 12, 13, 32, 14, 15, 16, 17, 30, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29})
}
