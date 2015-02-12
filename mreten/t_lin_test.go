// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
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

	model := GetModel("testsim", "mat1", "lin", false)
	model.Init(fun.Prms{
		&fun.Prm{N: "lam", V: 0.5},
		&fun.Prm{N: "pcae", V: 0.2},
		&fun.Prm{N: "slres", V: 0.1},
	})
	mdl := model.(Direct)
	utl.Pforan("mdl = %v\n", mdl)

	// plot
	if false {
		Pc := utl.LinSpace(-1, 3, 101)
		Sl := make([]float64, len(Pc))
		for i, pc := range Pc {
			Sl[i] = mdl.Sl(pc)
		}
		plt.Plot(Pc, Sl, "'b.-', label='lin', clip_on=0")
		plt.Gll("$p_c$", "$s_{\\ell}$", "")
		plt.Cross()
		plt.Show()
	}
}
