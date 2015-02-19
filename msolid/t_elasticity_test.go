// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

func check_constants(tst *testing.T, E, ν, Kcor, Gcor, lcor float64) {

	K, G := Calc_K_from_Enu(E, ν), Calc_G_from_Enu(E, ν)
	l := Calc_l_from_Enu(E, ν)

	io.Pf("E = %v\n", E)
	io.Pf("ν = %v\n", ν)
	io.Pf("K = %v\n", K)
	io.Pf("G = %v\n", G)
	io.Pf("l = %v\n", l)

	chk.Scalar(tst, "KfromEν", 1e-17, K, Kcor)
	chk.Scalar(tst, "GfromEν", 1e-17, G, Gcor)
	chk.Scalar(tst, "lfromEν", 1e-17, l, lcor)

	EfromKG, νfromKG := Calc_E_from_KG(K, G), Calc_nu_from_KG(K, G)
	chk.Scalar(tst, "EfromKG", 1e-17, EfromKG, E)
	chk.Scalar(tst, "νfromKG", 1e-17, νfromKG, ν)

	EfromlG, νfromlG := Calc_E_from_lG(l, G), Calc_nu_from_lG(l, G)
	chk.Scalar(tst, "EfromlG", 1e-17, EfromlG, E)
	chk.Scalar(tst, "νfromlG", 1e-17, νfromlG, ν)

	EfromKν, GfromKν := Calc_E_from_Knu(K, ν), Calc_G_from_Knu(K, ν)
	chk.Scalar(tst, "EfromKν", 1e-17, EfromKν, E)
	chk.Scalar(tst, "GfromKν", 1e-17, GfromKν, G)
}

func Test_elast01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	chk.PrintTitle("elast01")

	check_constants(tst, 60, 0.25, 40, 24, 24)
	//check_constants(tst, 600, 0.2, 1000.0/3, 250, 500.0/3.0)
}

func Test_elast02(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	chk.PrintTitle("elast02")

	K, G := 2.0, 3.0/4.0

	ndim, pstress := 2, false
	var ec SmallElasticity
	err := ec.Init(ndim, pstress, []*fun.Prm{
		&fun.Prm{N: "K", V: K},
		&fun.Prm{N: "G", V: G},
	})
	io.Pforan("ec: %v\n", &ec)
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}

	nsig, nalp, nphi, large := 2*ndim, 0, 0, false
	state := NewState(nsig, nalp, nphi, large)

	D := la.MatAlloc(nsig, nsig)
	ec.CalcD(D, state)

	a := K + 4.0*G/3.0
	b := K - 2.0*G/3.0
	c := 2.0 * G
	chk.Matrix(tst, "D", 1e-15, D, [][]float64{
		{a, b, b, 0},
		{b, a, b, 0},
		{b, b, a, 0},
		{0, 0, 0, c},
	})
}
