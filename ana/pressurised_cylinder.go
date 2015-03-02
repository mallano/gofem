// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package ana implements analytical solutions
package ana

import (
	"math"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/num"
)

// CteStressPstrain computes the constant-stress solution to a simple
// linear elastic plane-strain problem
//
//               , - - ,
//           , '         ' ,
//         ,                 ,
//        ,      .-'''-.      ,
//       ,      / ↖ ↑ ↗ \      ,
//       ,     |  ← P →  |     ,
//       ,      \ ↙ ↓ ↘ /      ,
//        ,      `-...-'      ,
//         ,                 ,
//           ,            , '
//             ' - , ,  '

type Hill struct {
	// Data
	a  float64 // Inner radius
	b  float64 // Outer radius
	E  float64 // Young's modulus
	ν  float64 // Poisson's coefficient
	σy float64 // Uniaxial yield stress
	P  float64 // Pressure prescribed on the inner surface

	// Derived data
	coef float64
	Y    float64
	P0   float64
	Plim float64
}

var tmp *Hill // temporary value (necessary in functions)

// Init initialises this structure
func (o *Hill) Init(prms fun.Prms) {
	// default values
	o.a = 100      // [mm]
	o.b = 200      // [mm]
	o.E = 210000.0 // Young modulus
	o.ν = 0.3      // Poisson ratio

	// derived
	o.coef = o.a * o.a / (o.b * o.b)
	o.Y = 2.0 * o.σy / math.Sqrt(3.0)
	o.P0 = o.Y * (1 - o.coef) / 2.0
	o.Plim = o.Y * math.Log(o.b/o.a)

	tmp = o // saves to temporary
}

// Solution computes solution
func (o Hill) Solution(t float64) (σx, σy, σz, εx, εy float64) {
	return
}

// Elastic solution for the radial displacement
func (o Hill) ub_e() float64 {
	return 2.0 * o.P * o.b * (1.0 - o.ν*o.ν) / (o.E/o.coef - o.E)
}

// Plastic solution for the radial displacment
func (o Hill) plas(c float64) (float64, float64) {
	var P float64
	P = o.Y * (math.Log(c/o.a) + (1.0-c*c/(o.b*o.b))/2.0)
	ub := o.Y * c * c * (1.0 - o.ν*o.ν) / (o.E * o.b)
	return P, ub
}

func (o Hill) Getc(P float64) float64 {
	var nls num.NlSolver
	defer nls.Clean()

	// function to be solved
	fx := func(fx, X []float64) (err error) {
		x := X[0]
		fx[0] = P/o.Y - (math.Log(x/o.a) + (1-x*x/(o.b*o.b))/2)
		return
	}
	// derivative function of f
	dfdx := func(dfdx [][]float64, X []float64) (err error) {
		x := X[0]
		dfdx[0][0] = -1.0/x + x/(o.b*o.b)
		return
	}

	Res := make([]float64, 1)
	nls.Init(1, fx, nil, dfdx, true, false, nil)
	nls.Solve(Res, false)
	return Res[0]
}

func (o Hill) sig(r, c float64) (float64, float64) { //_sig
	b := o.b
	Y := o.Y
	var sr, st float64
	if r > c { // elastic
		sr = -Y * c * c * (b*b/(r*r) - 1.0) / (2.0 * b * b)
		st = Y * c * c * (b*b/(r*r) + 1.0) / (2.0 * b * b)
	} else {
		sr = Y * (-0.5 - math.Log(c/r) + c*c/(2.0*b*b))
		st = Y * (0.5 - math.Log(c/r) + c*c/(2.0*b*b))
	}
	return sr, st
}

func (o Hill) Sig(R, C, SR, ST []float64) {
	for i := 0; i < len(R); i++ {
		SR[i], ST[i] = o.sig(R[i], C[i])
	}
}
