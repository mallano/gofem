// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_porous01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("porous01")

	mdb := ReadMat("data", "porous.mat")
	if mdb == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pf("porous.mat just read:\n%v\n", mdb)

	mat := mdb.Get("porous1")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pforan("mat = %+v\n", mat)

	cnd := mdb.GroupGet("porous1", "c")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pfcyan("cnd = %+v\n", cnd)

	lrm := mdb.GroupGet("porous1", "l")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pforan("lrm = %+v\n", lrm)

	por := mdb.GroupGet("porous1", "p")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pfcyan("por = %+v\n", por)

	sld := mdb.GroupGet("porous1", "s")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pforan("sld = %+v\n", sld)
}

func Test_porous02(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("porous02")

	mdb := ReadMat("data", "porous.mat")
	if mdb == nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pf("porous.mat just read:\n%v\n", mdb)

	cnd, lrm, por, err := mdb.GroupGet3("porous1", "c", "l", "p")
	if err != nil {
		tst.Errorf("test failed\n")
		return
	}
	utl.Pfcyan("cnd = %+v\n", cnd)
	utl.Pforan("lrm = %+v\n", lrm)
	utl.Pfcyan("por = %+v\n", por)
}
