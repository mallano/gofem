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

func onequa_solution(t float64) (σx, σy, σz, εx, εy float64) {
	qnV := -100.0
	qnH := -50.0
	E := 1000.0
	ν := 0.25
	σx = qnH * t
	σy = qnV * t
	σz = ν * (σx + σy)
	εx = (σx - ν*(σy+σz)) / E
	εy = (σy - ν*(σz+σx)) / E
	return
}

func onequa_check_u(tst *testing.T, t float64, u, x []float64, tol float64) {

	// analytical solution
	_, _, _, εx, εy := onequa_solution(t)

	// check displacements
	lx := 1.0
	ly := 1.0
	ux := lx * εx * x[0] / lx
	uy := ly * εy * x[1] / ly
	utl.CheckScalar(tst, "ux", tol, u[0], ux)
	utl.CheckScalar(tst, "uy", tol, u[1], uy)
}

func onequa_check_sig(tst *testing.T, t float64, σ, x []float64, tol float64) {

	// analytical solution
	σx, σy, σz, _, _ := onequa_solution(t)

	// check stresses
	utl.CheckScalar(tst, "σx ", tol, σ[0], σx)
	utl.CheckScalar(tst, "σy ", tol, σ[1], σy)
	utl.CheckScalar(tst, "σz ", tol, σ[2], σz)
	utl.CheckScalar(tst, "σxy", tol, σ[3], 0)
}

func Test_out01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
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
	verts := Verts{0, 1, 2, 3}
	cells := Cells{{0, 0}, {0, 1}, {-1, 2}, {-1, 3}}
	Tseries("ux", verts, nil)
	Tseries("uy", verts, nil)
	Tseries("sx", cells, nil)
	Tseries("sy", cells, nil)
	Tseries("sz", cells, nil)
	Tseries("sxy", cells, nil)

	// check slices
	nnod := 4
	nele := 1
	nip := 4
	utl.IntAssert(len(Dom.Nodes), nnod)
	utl.IntAssert(len(Ipoints), nele*nip)
	utl.IntAssert(len(Cid2ips), 1)
	utl.IntAssert(len(TseriesKeys), 6)
	utl.IntAssert(len(TseriesData), 6)
	utl.IntAssert(len(TseriesKey2idx), 6)
	utl.CompareStrs(tst, "TplotKeys", TseriesKeys, []string{"ux", "uy", "sx", "sy", "sz", "sxy"})
	utl.IntAssert(TseriesKey2idx["ux"], 0)
	utl.IntAssert(TseriesKey2idx["uy"], 1)
	utl.IntAssert(TseriesKey2idx["sx"], 2)
	utl.IntAssert(TseriesKey2idx["sy"], 3)
	utl.IntAssert(TseriesKey2idx["sz"], 4)
	utl.IntAssert(TseriesKey2idx["sxy"], 5)

	// check quantities
	for i, dat := range TseriesData {
		key := TseriesKeys[i]
		if key == "ux" || key == "uy" {
			utl.IntAssert(len(dat.Qts), 4)
			utl.IntAssert(len(dat.Sty), 4)
		}
		if key == "sx" || key == "sy" || key == "sz" || key == "sxy" {
			utl.IntAssert(len(dat.Qts), 4)
			utl.IntAssert(len(dat.Sty), 4)
		}
	}

	// apply commands
	err := Apply()
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
	}

	// check displacements
	tolu := 1e-15
	iux := TseriesKey2idx["ux"]
	iuy := TseriesKey2idx["uy"]
	for i, q := range TseriesData[iux].Qts {
		utl.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			utl.Pfyel("t=%g\n", t)
			u := []float64{
				TseriesRes[iux][i][j],
				TseriesRes[iuy][i][j],
			}
			onequa_check_u(tst, t, u, q.X, tolu)
		}
	}

	// check stresses
	tolσ := 1e-14
	isx := TseriesKey2idx["sx"]
	isy := TseriesKey2idx["sy"]
	isz := TseriesKey2idx["sz"]
	isxy := TseriesKey2idx["sxy"]
	for i, q := range TseriesData[isx].Qts {
		utl.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			utl.Pfyel("t=%g\n", t)
			σ := []float64{
				TseriesRes[isx][i][j],
				TseriesRes[isy][i][j],
				TseriesRes[isz][i][j],
				TseriesRes[isxy][i][j],
			}
			onequa_check_sig(tst, t, σ, q.X, tolσ)
		}
	}

	// show figure
	if !utl.Tsilent {
		Show(nil)
	}
}

func Test_out02(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("out02")

	// run FE simulation
	if !fem.Start("data/twoqua4.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
		return
	}
	defer fem.End()
	if !fem.Run() {
		tst.Errorf("test failed\n")
		return
	}

	// load results
	if !Start("data/twoqua4.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// verts: commands for reading time-series
	verts := Verts{-1}
	Tseries("ux", verts, nil)
	Tseries("uy", verts, nil)

	// get second ip coordinates
	xip := Ipoints[1].X
	utl.Pfcyan("xip = %v\n", xip)

	// cells: commands for reading time-series
	var cells Locator
	var nkeys, nipsel int
	testnumber := 1
	switch testnumber {
	case 1:
		cells = Along{[]float64{xip[0], 0}, []float64{xip[0], 1}}
		nkeys = 6
		nipsel = 4
	default:
		cells = AllCells()
		nkeys = 6
		nipsel = 8
	}
	Tseries("sx", cells, nil)
	Tseries("sy", cells, nil)
	Tseries("sz", cells, nil)
	Tseries("sxy", cells, nil)

	utl.Pforan("keys = %v\n", TseriesKeys)

	// check slices
	nnod := 6
	nele := 2
	nip := 4
	utl.IntAssert(len(Dom.Nodes), nnod)
	utl.IntAssert(len(Ipoints), nele*nip)
	utl.IntAssert(len(Cid2ips), 2)
	utl.IntAssert(len(TseriesKeys), nkeys)
	utl.IntAssert(len(TseriesData), nkeys)
	utl.IntAssert(len(TseriesKey2idx), nkeys)
	utl.CompareStrs(tst, "TplotKeys", TseriesKeys, []string{"ux", "uy", "sx", "sy", "sz", "sxy"})
	utl.IntAssert(TseriesKey2idx["ux"], 0)
	utl.IntAssert(TseriesKey2idx["uy"], 1)
	utl.IntAssert(TseriesKey2idx["sx"], 2)
	utl.IntAssert(TseriesKey2idx["sy"], 3)
	utl.IntAssert(TseriesKey2idx["sz"], 4)
	utl.IntAssert(TseriesKey2idx["sxy"], 5)

	// check quantities
	for i, dat := range TseriesData {
		key := TseriesKeys[i]
		if key == "ux" || key == "uy" {
			utl.IntAssert(len(dat.Qts), 1)
			utl.IntAssert(len(dat.Sty), 1)
		}
		if key == "sx" || key == "sy" || key == "sz" || key == "sxy" {
			utl.IntAssert(len(dat.Qts), nipsel)
			utl.IntAssert(len(dat.Sty), nipsel)
		}
	}

	// apply commands
	err := Apply()
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
	}

	// check displacements
	tolu := 1e-15
	iux := TseriesKey2idx["ux"]
	iuy := TseriesKey2idx["uy"]
	for i, q := range TseriesData[iux].Qts {
		utl.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			utl.Pfyel("t=%g\n", t)
			u := []float64{
				TseriesRes[iux][i][j],
				TseriesRes[iuy][i][j],
			}
			onequa_check_u(tst, t, u, q.X, tolu)
		}
	}

	// check stresses
	tolσ := 1e-13
	isx := TseriesKey2idx["sx"]
	isy := TseriesKey2idx["sy"]
	isz := TseriesKey2idx["sz"]
	isxy := TseriesKey2idx["sxy"]
	for i, q := range TseriesData[isx].Qts {
		utl.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			utl.Pfyel("t=%g\n", t)
			σ := []float64{
				TseriesRes[isx][i][j],
				TseriesRes[isy][i][j],
				TseriesRes[isz][i][j],
				TseriesRes[isxy][i][j],
			}
			onequa_check_sig(tst, t, σ, q.X, tolσ)
		}
	}

	// show figure
	if !utl.Tsilent {
		//Show(nil)
	}
}

// this test needs 'fem' package to be tested first
func test_out03(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("out03")

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
	Tseries("pl", At{2.5, 0}, nil)
	Tseries("pl", At{2.5, 10}, nil)
	Tseries("sl", At{xip[0], xip[1]}, nil)

	// check slices
	nnod := 27
	nele := 4
	nip := 4
	utl.IntAssert(len(Dom.Nodes), nnod)
	utl.IntAssert(len(Ipoints), nele*nip)
	utl.IntAssert(len(TseriesKeys), 2)
	utl.IntAssert(len(TseriesData), 2)
	utl.IntAssert(len(TseriesKey2idx), 2)
	utl.CompareStrs(tst, "TplotKeys", TseriesKeys, []string{"pl", "sl"})
	utl.IntAssert(TseriesKey2idx["pl"], 0)
	utl.IntAssert(TseriesKey2idx["sl"], 1)

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
func test_out04(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("out04")

	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	if !Start(datadir+"p01.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// get first ip coordinates
	ip := Ipoints[0].X

	// commands for reading time-series
	Tseries("pl", &At{2.5, 0}, Styles{{Label: "A", Marker: "o"}})
	Tseries("pl", &At{2.5, 10}, Styles{{Label: "B"}})
	Tseries("sl", &At{ip[0], ip[1]}, Styles{{Label: "ip"}})

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
		})
	}
}
