// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
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
	chk.Scalar(tst, "ux", tol, u[0], ux)
	chk.Scalar(tst, "uy", tol, u[1], uy)
}

func onequa_check_sig(tst *testing.T, t float64, σ, x []float64, tol float64) {

	// analytical solution
	σx, σy, σz, _, _ := onequa_solution(t)

	// check stresses
	chk.Scalar(tst, "σx ", tol, σ[0], σx)
	chk.Scalar(tst, "σy ", tol, σ[1], σy)
	chk.Scalar(tst, "σz ", tol, σ[2], σz)
	chk.Scalar(tst, "σxy", tol, σ[3], 0)
}

func Test_out01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("out01")

	// run FE simulation
	if !fem.Start("data/onequa4.sim", true, chk.Verbose) {
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
	chk.IntAssert(len(Dom.Nodes), nnod)
	chk.IntAssert(len(Ipoints), nele*nip)
	chk.IntAssert(len(Cid2ips), 1)
	chk.IntAssert(len(TseriesKeys), 6)
	chk.IntAssert(len(TseriesData), 6)
	chk.IntAssert(len(TseriesKey2idx), 6)
	chk.Strings(tst, "TplotKeys", TseriesKeys, []string{"ux", "uy", "sx", "sy", "sz", "sxy"})
	chk.IntAssert(TseriesKey2idx["ux"], 0)
	chk.IntAssert(TseriesKey2idx["uy"], 1)
	chk.IntAssert(TseriesKey2idx["sx"], 2)
	chk.IntAssert(TseriesKey2idx["sy"], 3)
	chk.IntAssert(TseriesKey2idx["sz"], 4)
	chk.IntAssert(TseriesKey2idx["sxy"], 5)

	// check quantities
	for i, dat := range TseriesData {
		key := TseriesKeys[i]
		if key == "ux" || key == "uy" {
			chk.IntAssert(len(dat.Qts), 4)
			chk.IntAssert(len(dat.Sty), 4)
		}
		if key == "sx" || key == "sy" || key == "sz" || key == "sxy" {
			chk.IntAssert(len(dat.Qts), 4)
			chk.IntAssert(len(dat.Sty), 4)
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
		io.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			io.Pfyel("t=%g\n", t)
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
		io.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			io.Pfyel("t=%g\n", t)
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
	if chk.Verbose {
		Show(nil)
	}
}

func Test_out02(tst *testing.T) {

	//verbose()
	chk.PrintTitle("out02")

	// run FE simulation
	if !fem.Start("data/twoqua4.sim", true, chk.Verbose) {
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
	io.Pfcyan("xip = %v\n", xip)

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

	io.Pforan("keys = %v\n", TseriesKeys)

	// check slices
	nnod := 6
	nele := 2
	nip := 4
	chk.IntAssert(len(Dom.Nodes), nnod)
	chk.IntAssert(len(Ipoints), nele*nip)
	chk.IntAssert(len(Cid2ips), 2)
	chk.IntAssert(len(TseriesKeys), nkeys)
	chk.IntAssert(len(TseriesData), nkeys)
	chk.IntAssert(len(TseriesKey2idx), nkeys)
	chk.Strings(tst, "TplotKeys", TseriesKeys, []string{"ux", "uy", "sx", "sy", "sz", "sxy"})
	chk.IntAssert(TseriesKey2idx["ux"], 0)
	chk.IntAssert(TseriesKey2idx["uy"], 1)
	chk.IntAssert(TseriesKey2idx["sx"], 2)
	chk.IntAssert(TseriesKey2idx["sy"], 3)
	chk.IntAssert(TseriesKey2idx["sz"], 4)
	chk.IntAssert(TseriesKey2idx["sxy"], 5)

	// check quantities
	for i, dat := range TseriesData {
		key := TseriesKeys[i]
		if key == "ux" || key == "uy" {
			chk.IntAssert(len(dat.Qts), 1)
			chk.IntAssert(len(dat.Sty), 1)
		}
		if key == "sx" || key == "sy" || key == "sz" || key == "sxy" {
			chk.IntAssert(len(dat.Qts), nipsel)
			chk.IntAssert(len(dat.Sty), nipsel)
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
		io.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			io.Pfyel("t=%g\n", t)
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
		io.PfWhite("\nx=%g\n", q.X)
		for j, t := range TseriesTimes {
			io.Pfyel("t=%g\n", t)
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
	if chk.Verbose {
		//Show(nil)
	}
}

// this test needs 'fem' package to be tested first
func test_out03(tst *testing.T) {

	//verbose()
	chk.PrintTitle("out03")

	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	if !Start(datadir+"p01.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// get ip
	n := len(Ipoints)
	ip := Ipoints[n-1].X
	io.Pfcyan("ip = %v\n", ip)

	// commands for reading time-series
	Tseries("pl", At{2.5, 0}, nil)
	Tseries("pl", At{2.5, 10}, nil)
	Tseries("sl", At{ip[0], ip[1]}, nil)

	// check slices
	nnod := 27
	nele := 4
	nip := 4
	chk.IntAssert(len(Dom.Nodes), nnod)
	chk.IntAssert(len(Ipoints), nele*nip)
	chk.IntAssert(len(TseriesKeys), 2)
	chk.IntAssert(len(TseriesData), 2)
	chk.IntAssert(len(TseriesKey2idx), 2)
	chk.Strings(tst, "TplotKeys", TseriesKeys, []string{"pl", "sl"})
	chk.IntAssert(TseriesKey2idx["pl"], 0)
	chk.IntAssert(TseriesKey2idx["sl"], 1)

	// check quantities
	for i, dat := range TseriesData {
		key := TseriesKeys[i]
		io.Pforan("key=%v => dat=\n%v\n", key, dat)
		if key == "pl" {
			chk.IntAssert(len(dat.Qts), 2)
			chk.IntAssert(len(dat.Sty), 2)
		}
		if key == "sl" {
			chk.IntAssert(len(dat.Qts), 1)
			chk.IntAssert(len(dat.Sty), 1)
		}
	}
}

// this test needs 'fem' package to be tested first
func test_out04(tst *testing.T) {

	//verbose()
	chk.PrintTitle("out04")

	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	if !Start(datadir+"p01.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// get last ip
	n := len(Ipoints)
	ip := Ipoints[n-1].X

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
	if chk.Verbose {
		Show(func() {
			plt.SubplotI(Spd["pl"])
			plt.AxisYrange(-10, 110)
		})
	}
}

// this test needs 'fem' package to be tested first
func test_out05(tst *testing.T) {

	chk.PrintTitle("out05")

	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	if !Start(datadir+"p02.sim", 0, 0) {
		tst.Errorf("Start failed\n")
		return
	}
	defer End()

	// get last ip
	n := len(Ipoints)
	ip := Ipoints[n-1].X

	// commands for reading time-series
	Tseries("pl", &At{2.5, 0}, nil)
	Tseries("pl", &At{2.5, 10}, nil)
	Tseries("sl", &At{ip[0], ip[1]}, nil)

	// apply commands
	err := Apply()
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
	}

	//e := Dom.Elems[3].(*fem.ElemP)
	//io.Pforan("sl @ ip = %v\n", e.States[3].Sl)

	// show figure
	if chk.Verbose {
		Show(nil)
	}
}
