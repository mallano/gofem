// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_mdl01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("mdl01")

	// info
	simfnk := "mdl01"
	matname := "mat1"
	getnew := false
	example := true

	// conductivity model
	cnd := mconduct.GetModel(simfnk, matname, "m1", getnew)
	err := cnd.Init(cnd.GetPrms(example))
	if err != nil {
		tst.Errorf("mconduct.Init failed: %v\n", err)
		return
	}

	// liquid retention model
	lrm := mreten.GetModel(simfnk, matname, "ref-m1", getnew)
	err = lrm.Init(lrm.GetPrms(example))
	if err != nil {
		tst.Errorf("mreten.Init failed: %v\n", err)
		return
	}

	// porous model
	mdl := GetModel(simfnk, matname, getnew)
	err = mdl.Init(mdl.GetPrms(example), cnd, lrm)
	if err != nil {
		tst.Errorf("mporous.Init failed: %v\n", err)
		return
	}

	// initial and final values
	pc0 := -5.0
	sl0 := 1.0
	pcf := 20.0

	// plot lrm
	doplot := true
	if doplot {
		npts := 41
		plt.Reset()
		mreten.Plot(mdl.Lrm, pc0, sl0, pcf, npts, "'b.-'", "'r+-'", "ref-m1_drying")
	}

	// state A
	var A StateLG
	pl0 := -5.0
	err = mdl.InitState(&A, pl0, 0)
	if err != nil {
		tst.Errorf("mporous.InitState failed: %v\n", err)
		return
	}

	// state B
	var B StateLG
	pl0 = -10.0
	err = mdl.InitState(&B, pl0, 0)
	if err != nil {
		tst.Errorf("mporous.InitState failed: %v\n", err)
		return
	}

	// show graph
	if doplot {
		plt.PlotOne(A.Pg-A.Pl, A.Sl, "'ro', label='A'")
		plt.PlotOne(B.Pg-B.Pl, B.Sl, "'go', label='B'")
		mreten.PlotEnd(true)
	}
}
