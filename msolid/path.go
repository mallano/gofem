// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"encoding/json"
	"math"

	"github.com/cpmech/gosl/tsr"
	"github.com/cpmech/gosl/utl"
)

// Path holds data for solid constitutive model simulations
type Path struct {

	// from json
	Sx    []float64 // σx stress components
	Sy    []float64 // σy stress components
	Sz    []float64 // σz components
	Ex    []float64 // εx strain components
	Ey    []float64 // εx strain components
	Ez    []float64 // εz strain components
	UseS  []int     // use stress component
	UseE  []int     // use strain component
	Nincs int       // number of increments
	Niout int       // number of increments for output
	MultS float64   // multiplier for stresses
	MultE float64   // multiplier for strains
	UseMS bool      // use MultS
	UseME bool      // use MultE

	// derived
	ndim int // space dimension
	ncp  int // number of stress components = 2 * ndim
	size int // number of path components
}

// Size returns the number of path components
func (o *Path) Size() int { return o.size }

// SetIsoCompS sets an isotropic compression path (stress driven)
func (o *Path) SetIsoCompS(ndim, nincs, niout int, P []float64) (err error) {

	// constants
	o.Nincs, o.Niout = nincs, niout

	// stress path
	n := len(P)
	o.Sx, o.Sy, o.Sz, o.UseS = make([]float64, n), make([]float64, n), make([]float64, n), make([]int, n)
	for i := 0; i < n; i++ {
		o.Sx[i], o.Sy[i], o.Sz[i], o.UseS[i] = -P[i], -P[i], -P[i], 1
	}

	// set additional information
	return o.init(ndim)
}

// SetPQstrain sets a p-q path with w=1 (compression); but given in terms of strains
func (o *Path) SetPQstrain(ndim, nincs, niout int, K, G, p0 float64, DP, DQ []float64, noise float64) (err error) {

	// constants
	o.Nincs, o.Niout = nincs, niout

	// strain path
	n := 1 + len(DP)
	o.Sx, o.Sy, o.Sz, o.UseS = []float64{-p0}, []float64{-p0}, []float64{-p0}, []int{0}
	o.Ex, o.Ey, o.Ez, o.UseE = make([]float64, n), make([]float64, n), make([]float64, n), make([]int, n)
	var Δp, Δq float64
	nσ := 4
	Δε := make([]float64, nσ)
	axsym := true
	o.UseE[0] = 1
	for i := 1; i < n; i++ {
		Δp, Δq = DP[i-1], DQ[i-1]
		CalcΔεElast(Δε, K, G, Δp, Δq, axsym)
		o.Ex[i], o.Ey[i], o.Ez[i], o.UseE[i] = o.Ex[i-1]+Δε[0], o.Ey[i-1]+Δε[1], o.Ez[i-1]+Δε[2], 1
		o.Ey[i] += noise
		o.Ez[i] -= noise
	}

	// set additional information
	return o.init(ndim)
}

// ReadJson reads json file
func (o *Path) ReadJson(ndim int, fname string) (err error) {

	// read file
	b, err := utl.ReadFile(fname)
	if err != nil {
		return utl.Err(_path_err11, fname, err)
	}

	// decode
	err = json.Unmarshal(b, o)
	if err != nil {
		return utl.Err(_path_err12, fname, err)
	}

	// set additional information
	return o.init(ndim)
}

// ReadTable loads path from datafile in table format
//  Note: n -- number of lines to read. use -1 to read all lines
func (o *Path) ReadTable(ndim, nincs, niout int, fname string, n int, mσ, mε float64, stresspath bool) (err error) {

	// constants
	o.Nincs, o.Niout = nincs, niout

	// read data table
	keys, d, err := utl.ReadTable(fname)
	if err != nil {
		return
	}
	if n < 0 {
		n = len(d[keys[0]])
	}

	// find x-y-z keys
	ex, ey, ez := "ex", "ey", "ez"
	sx, sy, sz := "sx", "sy", "sz"
	for _, key := range keys {
		switch key {
		case "Ex":
			ex, ey, ez = "Ex", "Ey", "Ez"
		case "Ea":
			ex, ey, ez = "Et", "Er", "Ea"
		case "ea":
			ex, ey, ez = "et", "er", "ea"
		case "Sx":
			sx, sy, sz = "Sx", "Sy", "Sz"
		case "Sa":
			sx, sy, sz = "St", "Sr", "Sa"
		case "sa":
			sx, sy, sz = "st", "sr", "sa"
		}
	}

	// set stress path
	if stresspath {
		o.Sx, o.Sy, o.Sz = make([]float64, n), make([]float64, n), make([]float64, n)
		for i := 0; i < n; i++ {
			o.Sx[i], o.Sy[i], o.Sz[i] = mσ*d[sx][i], mσ*d[sy][i], mσ*d[sz][i]
		}
		o.UseS = utl.IntVals(n, 1)

		// set strain path
	} else {
		o.Ex, o.Ey, o.Ez = make([]float64, n), make([]float64, n), make([]float64, n)
		for i := 0; i < n; i++ {
			o.Ex[i], o.Ey[i], o.Ez[i] = mε*d[ex][i], mε*d[ey][i], mε*d[ez][i]
		}
		o.Sx, o.Sy, o.Sz = []float64{mσ * d[sx][0]}, []float64{mσ * d[sy][0]}, []float64{mσ * d[sz][0]}
		o.UseE = utl.IntVals(n, 1)
	}

	// set additional information
	return o.init(ndim)
}

// init initialises states variables after {Sx, Sy, Sz} or {Ex, Ey, Ez} have been set
func (o *Path) init(ndim int) (err error) {

	// constants
	o.ndim = ndim
	o.ncp = 2 * ndim

	// size of slices and flags
	hasS, hasE := false, false
	allS, allE := true, true
	nSx, nSy, nSz := len(o.Sx), len(o.Sy), len(o.Sz)
	nEx, nEy, nEz := len(o.Ex), len(o.Ey), len(o.Ez)
	if nSx > 0 {
		hasS, allE = true, false
	}
	if nSy > 0 {
		hasS, allE = true, false
	}
	if nSz > 0 {
		hasS, allE = true, false
	}
	if nEx > 0 {
		hasE, allS = true, false
	}
	if nEy > 0 {
		hasE, allS = true, false
	}
	if nEz > 0 {
		hasE, allS = true, false
	}

	// check nS
	if nSx != nSy || nSx != nSz {
		return utl.Err(_path_err02, nSx, nSy, nSz)
	}

	// check for initial stresses
	if nSx < 1 || nSy < 1 || nSz < 1 {
		return utl.Err(_path_err03)
	}

	// unset hasS if only initial stress were given
	if nSx == 1 {
		hasS, allE = false, true
		if !hasE {
			return utl.Err(_path_err04)
		}
	}

	// other checks
	o.size = 0 // number of path components
	if hasS {
		o.size = nSx
	}
	if hasE {
		if nEx != nEy || nEx != nEz {
			return utl.Err(_path_err05, nEx, nEy, nEz)
		}
		o.size = nEx
	}
	if hasS && hasE {
		if nEx != nSx || nEy != nSx || nEz != nSx {
			return utl.Err(_path_err06)
		}
	}
	if !allS && !allE {
		if len(o.UseS) != nSx || len(o.UseE) != nSx {
			return utl.Err(_path_err07, len(o.UseS), len(o.UseE), nSx)
		}
	}

	// check size and Nincs
	if o.size < 2 {
		return utl.Err(_path_err08)
	}
	if o.Nincs < 1 {
		o.Nincs = 1
	}

	// multipliers
	if o.MultS < 1e-7 {
		o.MultS = 1
	}
	if o.MultE < 1e-7 {
		o.MultE = 1
	}

	// set use flags
	if allS {
		o.UseS = utl.IntVals(o.size, 1)
		o.UseE = make([]int, o.size)
	}
	if allE {
		o.UseE = utl.IntVals(o.size, 1)
		o.UseS = make([]int, o.size)
	}
	return
}

// CalcΔεElast calculates Δε corresponding to an elastic loading with Δp and Δq
func CalcΔεElast(Δε []float64, K, G float64, Δp, Δq float64, axsym bool) (Δεv, Δεd float64, err error) {
	Δεv = -Δp / K
	Δεd = Δq / (3.0 * G)
	var Δεx, Δεy, Δεz float64
	if axsym { // axisymmetric
		compression := true
		if compression {
			Δεx = Δεv/3.0 + Δεd/2.0
			Δεy = Δεx
			Δεz = Δεv/3.0 - Δεd
		} else {
			Δεx = Δεv/3.0 - Δεd/2.0
			Δεy = Δεv/3.0 + Δεd
			Δεz = Δεx
		}
	} else { // plane-strain with Δεy = Δεx / 2
		c := 9.0 * Δεd * Δεd / (4.0 * Δεv * Δεv)
		α := 0.0
		if math.Abs(c-1.0) > 1e-15 {
			d := 3.0 * (4.0*c - 1.0)
			if d < 0.0 {
				return 0, 0, utl.Err("discriminant < 0:  c=%v  d=%v", c, d)
			}
			α1 := (1.0 + 2.0*c + math.Sqrt(d)) / (2.0 - 2.0*c)
			α2 := (1.0 + 2.0*c - math.Sqrt(d)) / (2.0 - 2.0*c)
			α = α1
			utl.Pfyel("d, α1, α2 = %v, %v, %v\n", d, α1, α2)
		}
		utl.Pfyel("c, α = %v, %v\n", c, α)
		Δεy = Δεv / (1.0 + α)
		Δεx = α * Δεy
		Δεz = 0
	}
	//utl.Pfpink("Δp=%v, Δq=%v => Δεv=%v, Δεd=%v => Δεx=%v, Δεy=%v, Δεz=%v\n", Δp, Δq, Δεv, Δεd, Δεx, Δεy, Δεz)
	Δε[0] = Δεx
	Δε[1] = Δεy
	Δε[2] = Δεz
	Δε[3] = 0
	Δεv_ := tsr.M_εv(Δε)
	Δεd_ := tsr.M_εd(Δε)
	if math.Abs(Δεv-Δεv_) > 1e-15 {
		return 0, 0, utl.Err(_path_err09, Δεv, Δεv_)
	}
	if Δεd < 0 {
		Δεd_ = -Δεd_ // allow negative values
	}
	if math.Abs(Δεd-Δεd_) > 1e-15 {
		return 0, 0, utl.Err(_path_err10, Δεd, Δεd_)
	}
	return
}

// error messages
var (
	_path_err01 = "cannot read file <%s>\n"
	_path_err02 = "all S slices must have the same size. nSx=%d, nSy=%d, nSz=%d\n"
	_path_err03 = "at least one component of Sx,Sy,Sz must be given to initialise the stresses\n"
	_path_err04 = "with only initial stresses given, E slices must be given\n"
	_path_err05 = "all E slices must have the same size. nEx=%d, nEy=%d, nEz=%d\n"
	_path_err06 = "when using S and E slices at the same time, all {S,E} slices must have the same size\n"
	_path_err07 = "when using S and E slices, UseS and UseE must be given (with the same size as S and E). len(UseS)=%d, len(UseE)=%d, n{S,E}=%d\n"
	_path_err08 = "number of path components must be at least 2\n"
	_path_err09 = "failed on Δεv: %v ≠ %v\n"
	_path_err10 = "failed on Δεd: %v ≠ %v\n"
	_path_err11 = "cannot open file %v\n"
	_path_err12 = "cannot unmarshal file %v\n"
)
