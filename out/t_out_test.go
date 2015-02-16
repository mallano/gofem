// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_out01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("out01")

	stageIdx := 0
	regionIdx := 0
	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"

	defer With(datadir+"p01.sim", stageIdx, regionIdx)()
	Tplot("pl", &At{2.5, 0}, nil)
	Tplot("pl", &At{2.5, 10}, nil)
	Tplot("sl", &At{2.5, 10}, nil)

	nnod := 27
	nele := 4
	nip := 4
	utl.IntAssert(len(Dom.Nodes), nnod)
	utl.IntAssert(len(Ipoints), nele*nip)
	utl.IntAssert(len(TplotData), 2)
	utl.CompareStrs(tst, "TplotKeys", TplotKeys, []string{"pl", "sl"})

	err := Show()
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}
}

func Test_out02(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("out02")

	stageIdx := 0
	regionIdx := 0

	defer With("data/p01.sim", stageIdx, regionIdx)()
	Tplot("pl", &At{2.5, 10}, nil)
	Show()
}
