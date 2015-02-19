// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/plt"
)

func Test_mdl01(tst *testing.T) {

	doplot := false
	//doplot := true
	//utl.Tsilent = false
	chk.PrintTitle("mdl01")

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
	lrm_name := "ref-m1"
	//lrm_name := "vg"
	lrm := mreten.GetModel(simfnk, matname, lrm_name, getnew)
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
	//mdl.MEtrial = false
	mdl.ShowR = true

	// initial and final values
	pc0 := -5.0
	sl0 := 1.0
	pcf := 20.0

	// plot lrm
	if doplot {
		npts := 41
		plt.Reset()
		mreten.Plot(mdl.Lrm, pc0, sl0, pcf, npts, "'b.-'", "'r+-'", lrm_name)
	}

	// state A
	pl0 := -5.0
	pg, divus := 0.0, 0.0
	A, err := mdl.NewState(pl0, pg, divus)
	if err != nil {
		tst.Errorf("mporous.NewState failed: %v\n", err)
		return
	}

	// state B
	pl0 = -10.0
	B, err := mdl.NewState(pl0, pg, divus)
	if err != nil {
		tst.Errorf("mporous.NewState failed: %v\n", err)
		return
	}

	// plot A and B points
	if doplot {
		pcA := A.Pg - A.Pl
		pcB := B.Pg - B.Pl
		plt.PlotOne(pcA, A.Sl, "'gs', clip_on=0, label='A', ms=10")
		plt.PlotOne(pcB, B.Sl, "'ks', clip_on=0, label='B'")
	}

	// incremental update
	//Δpl := -20.0 // << problems with this one and VG
	Δpl := -5.0
	n := 23
	iwet := 10
	Pc := make([]float64, n)
	Sl := make([]float64, n)
	Pc[0] = A.Pg - A.Pl
	Sl[0] = A.Sl
	for i := 1; i < n; i++ {
		if i > iwet {
			Δpl = -Δpl
			iwet = n
		}
		err = mdl.Update(A, Δpl, pg, divus)
		if err != nil {
			tst.Errorf("test failed: %v\n", err)
			return
		}
		Pc[i] = A.Pg - A.Pl
		Sl[i] = A.Sl
	}

	// show graph
	if doplot {
		plt.Plot(Pc, Sl, "'ro-', clip_on=0, label='update'")
		mreten.PlotEnd(true)
	}
}
