// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package ana implements analytical solutions
package ana

import (
	"testing"

	"github.com/cpmech/gosl/fun"

	"github.com/cpmech/gosl/chk"
)

// CteStressPstrain computes the constant-stress solution to a simple
// linear elastic plane-strain problem
//
//              qnV
//        ↓↓↓↓↓↓↓↓↓↓↓↓↓
//      ▷ o-----------o ←
//      ▷ |           | ←
//      ▷ |     E     | ← qnH
//  ly  ▷ |     ν     | ←
//      ▷ |           | ←
//      ▷ o-----------o ←
//        △  △  △  △  △
//             lx
type CteStressPstrain struct {
	qnV float64 // vertical distributed load
	qnH float64 // horizontal distributed load
	E   float64 // Young's modulus
	ν   float64 // Poisson's coefficient
	lx  float64 // x-length
	ly  float64 // y-length
}

// Init initialises this structure
func (o *CteStressPstrain) Init(prms fun.Prms) {

	// default values
	o.qnV = -100.0
	o.qnH = -50.0
	o.E = 1000.0
	o.ν = 0.25
	o.lx = 1.0
	o.ly = 1.0

	// parameters
	for _, p := range prms {
		switch p.N {
		case "qnV":
			o.qnV = p.V
		case "qnH":
			o.qnH = p.V
		case "E":
			o.E = p.V
		case "nu":
			o.ν = p.V
		case "lx":
			o.lx = p.V
		case "ly":
			o.ly = p.V
		}
	}
}

// Solution computes solution
func (o CteStressPstrain) Solution(t float64) (σx, σy, σz, εx, εy float64) {
	σx = o.qnH * t
	σy = o.qnV * t
	σz = o.ν * (σx + σy)
	εx = (σx - o.ν*(σy+σz)) / o.E
	εy = (σy - o.ν*(σz+σx)) / o.E
	return
}

// CheckDispl checks displacements
func (o CteStressPstrain) CheckDispl(tst *testing.T, t float64, u, x []float64, tol float64) {

	// analytical solution
	_, _, _, εx, εy := o.Solution(t)

	// check displacements
	ux := o.lx * εx * x[0] / o.lx
	uy := o.ly * εy * x[1] / o.ly
	chk.Scalar(tst, "ux", tol, u[0], ux)
	chk.Scalar(tst, "uy", tol, u[1], uy)
}

// CheckStress check stresses
func (o CteStressPstrain) CheckStress(tst *testing.T, t float64, σ, x []float64, tol float64) {

	// analytical solution
	σx, σy, σz, _, _ := o.Solution(t)

	// check stresses
	chk.Scalar(tst, "σx ", tol, σ[0], σx)
	chk.Scalar(tst, "σy ", tol, σ[1], σy)
	chk.Scalar(tst, "σz ", tol, σ[2], σz)
	chk.Scalar(tst, "σxy", tol, σ[3], 0)
}
