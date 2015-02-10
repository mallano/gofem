// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_fileio01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("fileio01")

	// start
	Start("data/bh16.sim", true, !utl.Tsilent)
	defer End()

	// domain A
	domA := NewDomain(global.Sim.Regions[0])
	domA.SetStage(0, global.Sim.Stages[0])
	for i, _ := range domA.Sol.Y {
		domA.Sol.Y[i] = float64(i)
	}
	utl.Pforan("domA.Sol.Y = %v\n", domA.Sol.Y)

	// write file
	tidx := 123
	err := domA.SaveSol(tidx)
	if err != nil {
		tst.Fatalf("%v", err)
		return
	}
	utl.Pfblue2("file %v written\n", out_nod_path(tidx))

	// domain B
	domB := NewDomain(global.Sim.Regions[0])
	domB.SetStage(0, global.Sim.Stages[0])
	utl.Pfpink("domB.Sol.Y (before) = %v\n", domB.Sol.Y)

	// read file
	err = domB.ReadSol(tidx)
	if err != nil {
		tst.Fatalf("%v", err)
		return
	}
	utl.Pfgreen("domB.Sol.Y (after) = %v\n", domB.Sol.Y)

	// check
	utl.CheckVector(tst, "Y", 1e-17, domA.Sol.Y, domB.Sol.Y)
	utl.CheckVector(tst, "dy/dt", 1e-17, domA.Sol.Dydt, domB.Sol.Dydt)
	utl.CheckVector(tst, "d²y/dt²", 1e-17, domA.Sol.D2ydt2, domB.Sol.D2ydt2)
}
