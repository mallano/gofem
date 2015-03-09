// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_out01(tst *testing.T) {

	// finalise analysis process and catch errors
	defer func() {
		if err := recover(); err != nil {
			tst.Fail()
			io.PfRed("ERROR: %v\n", err)
		} else {
			fem.End()
		}
	}()

	// test title
	//verbose()
	chk.PrintTitle("out01")

	// run FE simulation
	if !fem.Start("data/onequa4.sim", true, chk.Verbose) {
		chk.Panic("cannot start FE simulation")
	}
	if !fem.Run() {
		chk.Panic("cannot run FE simulation")
	}

	// start analysis process
	Start("data/onequa4.sim", 0, 0)

	// define points
	Define("A B C D", N{0, 1, 2, 3})
	Define("a b c d", P{{0, 0}, {0, 1}, {-1, 2}, {-1, 3}})
	Define("!right side", Along{{1, 0}, {1, 1}})
	io.Pfcyan("Entities = %v\n", R)

	// check slices
	nnod := 4
	nele := 1
	nip := 4
	chk.IntAssert(len(Dom.Nodes), nnod)
	chk.IntAssert(len(Ipoints), nele*nip)
	chk.IntAssert(len(Cid2ips), 1)
	chk.IntAssert(len(Ipkey2ips), 4)
	chk.IntAssert(len(Ipkeys), 4)
	for key, ips := range Ipkey2ips {
		chk.Ints(tst, io.Sf("%s : ips", key), ips, utl.IntRange(4))
	}

	// load results
	LoadResults(nil)

	// check points
	nlabels := []string{"A", "B", "C", "D"}
	for _, l := range nlabels {
		if _, ok := R[l]; !ok {
			chk.Panic("1: %q alias in Entities is not available", l)
		}
	}
	plabels := []string{"a", "b", "c", "d"}
	for _, l := range plabels {
		if _, ok := R[l]; !ok {
			chk.Panic("2: %q alias in Entities is not available", l)
		}
	}
	if _, ok := R["right side"]; !ok {
		chk.Panic("3: %q alias in Entities is not available", "right side")
	}

	// check u-keys
	ukeys := []string{"ux", "uy"}
	for _, l := range nlabels {
		for _, p := range R[l] {
			//io.Pfyel("p = %v\n", p)
			if p == nil {
				chk.Panic("1: p is nil")
			}
			for _, key := range ukeys {
				if _, ok := p.Vals[key]; !ok {
					chk.Panic("%s is not available in point", key)
				}
			}
		}
	}

	// check s-keys
	skeys := fem.StressKeys()
	for _, l := range plabels {
		for _, p := range R[l] {
			//io.Pfgreen("q = %v\n", p)
			if p == nil {
				chk.Panic("2: p is nil")
			}
			for _, key := range skeys {
				if _, ok := p.Vals[key]; !ok {
					chk.Panic("%s is not available in point", key)
				}
			}
		}
	}

	// check GetRes
	uxC := GetRes("ux", "C", 0)
	uxR := GetRes("ux", "right side", -1)
	io.Pforan("uxC = %v\n", uxC)
	io.Pforan("uxR = %v\n", uxR)
	chk.IntAssert(len(uxC), 2)
	idx := len(uxC) - 1
	chk.Vector(tst, "uxR", 1e-17, uxR, []float64{uxC[idx], uxC[idx]})

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
		&fun.Prm{N: "qnH", V: -50},
		&fun.Prm{N: "qnV", V: -100},
	})

	// check displacements
	tolu := 1e-15
	for _, l := range nlabels {
		x := GetCoords(l)
		ux := GetRes("ux", l, 0)
		uy := GetRes("uy", l, 0)
		io.Pforan("ux=%v uy=%v\n", ux, uy)
		for j, t := range T {
			io.Pfyel("t=%g\n", t)
			sol.CheckDispl(tst, t, []float64{ux[j], uy[j]}, x, tolu)
		}
	}

	// check stresses
	tolσ := 1e-14
	for _, l := range plabels {
		x := GetCoords(l)
		sx := GetRes("sx", l, 0)
		sy := GetRes("sy", l, 0)
		sz := GetRes("sz", l, 0)
		sxy := GetRes("sxy", l, 0)
		for j, t := range T {
			io.Pfyel("t=%g\n", t)
			sol.CheckStress(tst, t, []float64{sx[j], sy[j], sz[j], sxy[j]}, x, tolσ)
		}
	}
}

func Test_out02(tst *testing.T) {

	// finalise analysis process and catch errors
	defer func() {
		if err := recover(); err != nil {
			tst.Fail()
			io.PfRed("ERROR: %v\n", err)
		} else {
			fem.End()
		}
	}()

	// test title
	//verbose()
	chk.PrintTitle("out02")

	// run FE simulation
	if !fem.Start("data/twoqua4.sim", true, chk.Verbose) {
		chk.Panic("cannot start FE simulation")
	}
	if !fem.Run() {
		chk.Panic("cannot run FE simulation")
	}

	// start analysis process
	Start("data/twoqua4.sim", 0, 0)

	// get second ip coordinates
	xip := Ipoints[1].X
	io.Pfcyan("xip = %v\n", xip)

	// define points
	Define("A", N{-1})
	Define("ips", Along{{xip[0], 0}, {xip[0], 1}})

	// load results
	LoadResults(nil)

	// solution
	var sol ana.CteStressPstrain
	sol.Init(fun.Prms{
		&fun.Prm{N: "qnH", V: -50},
		&fun.Prm{N: "qnV", V: -100},
	})

	// check displacements
	tolu := 1e-15
	x := GetCoords("A")
	ux := GetRes("ux", "A", 0)
	uy := GetRes("uy", "A", 0)
	io.Pforan("ux=%v uy=%v\n", ux, uy)
	for j, t := range T {
		io.Pfyel("t=%g\n", t)
		sol.CheckDispl(tst, t, []float64{ux[j], uy[j]}, x, tolu)
	}
}
