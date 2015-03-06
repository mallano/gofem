// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"code.google.com/p/gemlab"
	"code.google.com/p/gosl/utl"
)

func main() {

	var dat gemlab.InData

	B, H, W := 0.6, 3.0, 3.0

	coarse := true
	ufine := false
	dat.Quad = true
	dat.ToQua9 = true
	dat.Nparts = 4
	dat.Pfull = false

	y := H / 2.0
	b := B / 2.0
	l := b / 3.0
	a := b - l/2.0
	m := a + l
	c := a
	L := a + l + c
	e := l / 2.0
	f := a
	g := H - e - f
	n := f + g

	var M, N, O, P, Q, R int
	if coarse {
		M, N, O, P, Q, R = 2, 2, 3, 2, 2, 1
	} else {
		if ufine { // ultrafine
			M, N, O, P, Q, R = 8, 12, 12, 8, 8, 4
		} else {
			M, N, O, P, Q, R = 4, 6, 6, 4, 4, 2
		}
	}

	α, β, γ, δ := -0.6, -0.4, 0.8, -0.2

	dat.Sregs = &gemlab.Sregs{
		//         0  1  2  3  4  5  6  7  8
		[]int{-1, -1, -1, -1, -1, -1, -1, -1, -1}, // region tags
		[]int{N, Q, N, O, N, O, Q, N, Q},          // Nx
		[]int{R, R, O, M, M, M, M, P, P},          // Ny
		[]int{1, 1, 1, 1, 1, 1, 1, 1, 1},          // Nz
		[]float64{0, γ, 0, β, 0, -β, γ, 0, γ},     // Ax
		[]float64{δ, δ, β, 0, 0, 0, 0, α, α},      // Ay
		[]float64{0, 0, 0, 0, 0, 0, 0, 0, 0},      // Az
		[]int{0, 0, 0, 0, 0, 0, 0, 0, 0},          // NlX
		[]int{0, 0, 0, 0, 0, 0, 0, 0, 0},          // NlY
		[]int{0, 0, 0, 0, 0, 0, 0, 0, 0},          // NlZ
		[][]float64{ // points
			{0, 0}, //  0
			{L, 0}, //  1
			{W, 0}, //  2
			{0, g}, //  3
			{L, g}, //  4
			{W, g}, //  5
			{a, n}, //  6
			{m, n}, //  7
			{0, H}, //  8
			{a, H}, //  9
			{m, H}, // 10
			{L, H}, // 11
			{W, H}, // 12
			{0, y}, // 13
			{L, y}, // 14
			{W, y}, // 15
		},
		[][]int{ // regions connectivity
			{0, 1, 14, 13}, //  0
			{1, 2, 15, 14}, //  1
			{3, 4, 7, 6},   //  2
			{3, 6, 9, 8},   //  3
			{6, 7, 10, 9},  //  4
			{7, 4, 11, 10}, //  5
			{4, 5, 12, 11}, //  6
			{13, 14, 4, 3}, //  7
			{14, 15, 5, 4}, //  8
		},
		[][]int{ // boundary tags
			{-10, 0, -15, -13}, //  0
			{-10, -11, -15, 0}, //  1
			{0, 0, 0, 0},       //  2
			{0, 0, -14, -13},   //  3
			{0, 0, -16, 0},     //  4
			{0, 0, -12, 0},     //  5
			{0, -11, -12, 0},   //  6
			{0, 0, 0, -13},     //  7
			{0, -11, 0, 0},     //  8
		},
	}
	dat.Sregs.Draw(".", "indentD2_blocks")

	tol := 1e-5
	dat.VtagsL = &gemlab.VtagsL{
		[]int{-6, -300, -66}, // tags
		[][]float64{ // xxa
			{b, 0},
			{0, H},
			{0, 0},
		},
		[][]float64{ // xxb
			{b, H},
			{b, H},
			{0, H},
		},
		[]float64{tol, tol, tol}, // tols
	}

	tol = 1e-7
	dat.StagsL = &gemlab.StagsL{
		[]int{-14, -12}, // tags
		[][]float64{ // xxa
			{0, H},
			{0, H},
		},
		[][]float64{ // xxb
			{W, H},
			{W, H},
		},
		[][]float64{ // cmin
			{a, H},
			{b, H},
		},
		[][]float64{ // cmax
			{b, H},
			{m, H},
		},
		[]float64{tol, tol}, // tols
	}

	tol = 1e-3
	dat.Vtags = &gemlab.Vtags{
		[]int{-1, -2, -3, -4, -5, -11, -22, -33, -44, -55},
		[][]float64{
			{b, H},       //  -1
			{b, H - 0.1}, //  -2
			{b, 2.7},     //  -3
			{b, H / 2.0}, //  -4
			{b, 0},       //  -5
			{0, H},       // -11
			{0, H - 0.1}, // -12
			{0, 2.7},     // -13
			{0, H / 2.0}, // -14
			{0, 0},       // -15
		},
		[]float64{tol, tol, tol, tol, tol, tol, tol, tol, tol, tol},
	}

	tol = 0.01
	dat.Ctags = &gemlab.Ctags{
		[]int{-2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3},
		[][]float64{
			{0.28, 0.26},
			{0.27, 1.07},
			{0.27, 1.62},
			{0.27, 2.07},
			{0.27, 2.40},
			{0.28, 2.61},
			{0.29, 2.74},
			{0.29, 2.86},
			{0.29, 2.92},
			{0.295, 2.978},
			{0.00, 0.26},
			{0.00, 1.07},
			{0.00, 1.62},
			{0.00, 2.07},
			{0.00, 2.40},
			{0.00, 2.61},
			{0.00, 2.73},
			{0.00, 2.84},
			{0.00, 2.91},
			{0.00, 3.00},
		},
		[]float64{tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol, tol},
	}

	fnk := "square"
	if coarse {
		fnk += "-coarse"
	}
	if ufine {
		fnk += "-ufine"
	}
	if dat.ToQua9 {
		fnk += "-q9"
	}

	if err := gemlab.Generate(fnk, &dat); err != nil {
		utl.PfRed("%v\n", err.Error())
	}
}
