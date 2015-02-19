// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

func Test_vm01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	chk.PrintTitle("vm01")

	// allocate driver
	ndim, pstress := 2, false
	simfnk, modelname := "test", "vm"
	var drv Driver
	err := drv.Init(simfnk, modelname, ndim, pstress, []*fun.Prm{
		&fun.Prm{N: "K", V: 1.5},
		&fun.Prm{N: "G", V: 1},
		&fun.Prm{N: "qy0", V: 2},
		&fun.Prm{N: "H", V: 0.5},
	})
	drv.CheckD = true
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// vm model
	vm := drv.model.(*VonMises)

	// path
	p0 := 0.0
	Δp := 3.0
	Δq := vm.qy0
	ϵ := 1e-3
	DP := []float64{Δp + ϵ, 3, 2, 1, 0}
	DQ := []float64{Δq + ϵ, 4, 2, 1, 3}
	nincs := 1
	niout := 1
	noise := 0.0
	var pth Path
	err = pth.SetPQstrain(ndim, nincs, niout, vm.K, vm.G, p0, DP, DQ, noise)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// run
	err = drv.Run(&pth)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// plot
	//if DPsaveFig {
	//}
}
