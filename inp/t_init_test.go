// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"os"
	"testing"

	"github.com/cpmech/gosl/utl"
)

func init() {
	utl.Tsilent = true
}

func Test_msh01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mSome error has happened:[0m\n", err)
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("msh01")

	msh := ReadMsh("data/bh16.msh", !utl.Tsilent)
	utl.Pforan("%v\n", msh)
}

func Test_sim01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mSome error has happened:[0m\n", err)
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("sim01")

	sim := ReadSim("data/bh16.sim", true, !utl.Tsilent)
	if !utl.Tsilent {
		sim.GetInfo(os.Stdout)
		utl.Pf("\n")
	}
}

func Test_mat01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mSome error has happened:[0m\n", err)
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("mat01")

	mdb1 := ReadMat("data/bh.mat", !utl.Tsilent)
	utl.Pforan("bh.mat just read:\n%v\n", mdb1)

	fn := "test_bh.mat"
	utl.WriteFileSD("/tmp/gofem/inp", fn, mdb1.String())

	mdb2 := ReadMat("/tmp/gofem/inp/"+fn, !utl.Tsilent)
	utl.Pfblue2("\n%v\n", mdb2)
}

func Test_mat02(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mSome error has happened:[0m\n", err)
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("mat02 (conversion)")

	convertsymbols := true
	MatfileOld2New("/tmp/gofem/inp", "new_layers.mat", "data/old_layers.mat", convertsymbols)

	mdb := ReadMat("/tmp/gofem/inp/new_layers.mat", !utl.Tsilent)
	utl.Pfblue2("%v\n", mdb)
}
