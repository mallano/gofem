// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_derivs01(tst *testing.T) {

	doplot := false
	//doplot := true
	//utl.Tsilent = false
	chk.PrintTitle("derivs01")

	// info
	simfnk := "derivs01"
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
	//mdl.Ncns = true
	//mdl.Ncns2 = true

	// path
	pc0 := 0.0
	pcf := 20.0
	np := 5
	//P := []float64{10, 5, 20, 0}
	P := []float64{5}
	pth := GetPathCycle(pc0, P, np)
	io.Pforan("pth = %v\n", pth)

	// driver
	var drv Driver
	err = drv.Init(mdl)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}
	err = drv.Run(pth)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	// plot
	if doplot {
		npts := 41
		plt.Reset()
		mreten.Plot(mdl.Lrm, pc0, 1.0, pcf, npts, "'b.-'", "'r+-'", lrm_name)
		n := len(drv.Res)
		Pc := make([]float64, n)
		Sl := make([]float64, n)
		for i, s := range drv.Res {
			Pc[i] = s.Pg - s.Pl
			Sl[i] = s.Sl
		}
		plt.Plot(Pc, Sl, "'ko--', clip_on=0")
		mreten.PlotEnd(true)
	}
}
