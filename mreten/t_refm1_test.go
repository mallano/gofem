// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_refm1a(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("refm1a")

	mdl := GetModel("testsim", "mat1", "ref-m1", false)
	mdl.Init(mdl.GetPrms())

	pc0 := -5.0
	sl0 := 1.0
	pcf := 20.0
	npts := 41

	doplot := true
	if doplot {
		plt.Reset()
		Plot(mdl, pc0, sl0, pcf, npts, "'b.-'", "'r+-'", "ref-m1")
	}

	npts = 11
	tolCc, tolD1, tolD2 := 1e-10, 1e-10, 1e-10
	Check(tst, mdl, pc0, sl0, pcf, npts, tolCc, tolD1, tolD2, true, []float64{}, 1e-7, doplot)

	if doplot {
		PlotEnd(true)
	}
}
