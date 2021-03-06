// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_porous01(tst *testing.T) {

	chk.PrintTitle("porous01")

	mdb := ReadMat("data", "porous.mat")
	if mdb == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pf("porous.mat just read:\n%v\n", mdb)

	mat := mdb.Get("porous1")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pforan("mat = %+v\n", mat)

	cnd := mdb.GroupGet("porous1", "c")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pfcyan("cnd = %+v\n", cnd)

	lrm := mdb.GroupGet("porous1", "l")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pforan("lrm = %+v\n", lrm)

	por := mdb.GroupGet("porous1", "p")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pfcyan("por = %+v\n", por)

	sld := mdb.GroupGet("porous1", "s")
	if mat == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pforan("sld = %+v\n", sld)
}

func Test_porous02(tst *testing.T) {

	chk.PrintTitle("porous02")

	mdb := ReadMat("data", "porous.mat")
	if mdb == nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pf("porous.mat just read:\n%v\n", mdb)

	cnd, lrm, por, err := mdb.GroupGet3("porous1", "c", "l", "p")
	if err != nil {
		tst.Errorf("test failed\n")
		return
	}
	io.Pfcyan("cnd = %+v\n", cnd)
	io.Pforan("lrm = %+v\n", lrm)
	io.Pfcyan("por = %+v\n", por)
}
