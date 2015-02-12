// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_lin01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("lin01")

	mdl := GetModel("testsim", "mat1", "lin", false)
	mdl.Init(mdl.GetPrms())

	// check derivatives
	doplot := true
	tolCc, tolD1, tolD2 := 1e-10, 1e-10, 1e-10
	Check(tst, mdl, -0.5, 1, 3, 11, tolCc, tolD1, tolD2, true, []float64{0.2}, 1e-7, doplot)

	// plot
	if false {
		Plot(mdl, -0.5, 1, 3, 11, "'b.-'", "'r+-'", "lin")
		PlotEnd(true)
	}
}
