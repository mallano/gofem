// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func onequa_solution(tst *testing.T, t float64, dom *fem.Domain, tolu, tolσ float64) {

	// analytical solution
	qnV, qnH := -100.0, -50.0
	E, ν := 1000.0, 0.25
	lx, ly := 1.0, 1.0
	σx, σy := qnH*t, qnV*t
	σz := ν * (σx + σy)
	εx := (σx - ν*(σy+σz)) / E
	εy := (σy - ν*(σz+σx)) / E

	// check displacements
	ux := lx * εx
	uy := ly * εy
	Ycor := []float64{0, 0, ux, 0, ux, uy, 0, uy}
	utl.CheckVector(tst, "Y", tolu, dom.Sol.Y, Ycor)

	// check stresses
	e := dom.Elems[0].(*fem.ElemU)
	σcor := []float64{σx, σy, σz, 0}
	for idx, _ := range e.IpsElem {
		utl.CheckVector(tst, "σ", tolσ, e.States[idx].Sig, σcor)
	}
}

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

	// run FE simulation
	if !fem.Start("data/onequa4.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
		return
	}
	defer fem.End()
	if !fem.Run() {
		tst.Errorf("test failed\n")
		return
	}

	// load results
	if !Start("data/onequa4.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// commands for reading time-series
	Tseries("ux", &IdsOrTags{0, 1, 2, 3}, nil)

	// apply commands
	err := Apply()
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
	}

	// check FE simulation results
	utl.Pforan("T = %v\n", TseriesT)
	utl.Pforan("R = %v\n", TseriesR)
}

// this test needs 'fem' package to be tested first
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

	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	if !Start(datadir+"p01.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// get first ip coordinates
	xip := Ipoints[0].X
	utl.Pfcyan("xip = %v\n", xip)

	// commands for reading time-series
	Tseries("pl", &At{2.5, 0}, nil)
	Tseries("pl", &At{2.5, 10}, nil)
	Tseries("sl", &At{xip[0], xip[1]}, nil)

	// check slices
	nnod := 27
	nele := 4
	nip := 4
	utl.IntAssert(len(Dom.Nodes), nnod)
	utl.IntAssert(len(Ipoints), nele*nip)
	utl.IntAssert(len(TseriesKeys), 2)
	utl.IntAssert(len(TseriesData), 2)
	utl.CompareStrs(tst, "TplotKeys", TseriesKeys, []string{"pl", "sl"})

	// check quantities
	for i, dat := range TseriesData {
		key := TseriesKeys[i]
		utl.Pforan("key=%v => dat=%v\n", key, dat)
		if key == "pl" {
			utl.IntAssert(len(dat.Qts), 2)
			utl.IntAssert(len(dat.Sty), 2)
		}
		if key == "sl" {
			utl.IntAssert(len(dat.Qts), 1)
			utl.IntAssert(len(dat.Sty), 1)
		}
	}
}

// this test needs 'fem' package to be tested first
func Test_out03(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("out03")

	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	if !Start(datadir+"p01.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// commands for reading time-series
	Tseries("pl", &At{2.5, 0}, Styles{{Label: "A", Marker: "o"}})
	Tseries("pl", &At{2.5, 10}, Styles{{Label: "B"}})

	// apply commands
	err := Apply()
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
	}

	// show figure
	if !utl.Tsilent {
		Show(func() {
			plt.SubplotI(Spd["pl"])
			plt.AxisYrange(-10, 110)
			plt.Gll("$t$", "$p_{\\ell}$", "")
		})
	}
}
