// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_vg01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("vg01")

	model := GetModel("testsim", "mat1", "vg", false)
	model.Init(model.GetPrms())
	mdl := model.(Direct)
	utl.Pforan("mdl = %v\n", mdl)

	// plot
	if false {
		Pc := utl.LinSpace(-1, 3, 101)
		Sl := make([]float64, len(Pc))
		for i, pc := range Pc {
			Sl[i] = mdl.Sl(pc)
		}
		plt.Plot(Pc, Sl, "'b.-', label='vg', clip_on=0")
		plt.Gll("$p_c$", "$s_{\\ell}$", "")
		plt.Cross()
		plt.Show()
	}
}
