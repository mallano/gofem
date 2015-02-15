// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/utl"
)

func Test_node01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("node01")

	// new Vertice
	vert := inp.Vert{0, -1, []float64{0.0, 0.0, 0.0}}
	utl.Pforan("Vert = %v\n", vert)

	// new list of Dofs
	dofs := []*Dof{{"ux", 0}, {"uy", 1}, {"uz", 2}}

	// new Node with previous Dofs
	node0 := Node{dofs, &vert}
	utl.Pforan("node0= %v\n", node0)

	// getting dofs by key
	dof0 := node0.GetDof("ux")
	dof1 := node0.GetDof("uy")
	dof2 := node0.GetDof("uz")
	utl.Pforan("dof0= %v\n", dof0)
	utl.Pforan("dof1= %v\n", dof1)
	utl.Pforan("dof2= %v\n", dof2)

	utl.CompareStrs(tst, "Dof Ukey", []string{dof0.Key, dof1.Key, dof2.Key}, []string{"ux", "uy", "uz"})

	// getting Eq number by key
	eq0 := node0.GetEq("ux")
	eq1 := node0.GetEq("uy")
	eq2 := node0.GetEq("uz")
	utl.Pforan("eq0= %v\n", eq0)
	utl.Pforan("eq1= %v\n", eq1)
	utl.Pforan("eq2= %v\n", eq2)

	utl.CompareInts(tst, "Dof Eq", []int{eq0, eq1, eq2}, []int{0, 1, 2})
}
