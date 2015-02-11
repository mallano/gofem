// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/utl"
)

func Test_shape01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("Test shape01")

	for name, shape := range factory {

		utl.Pfyel("--------------------------------- %-6s---------------------------------\n", name)

		// check S
		tol := 1e-17
		errS := 0.0
		if name == "tri10" {
			tol = 1e-14
		}
		for n := 0; n < shape.Nverts; n++ {
			rst := []float64{0, 0, 0}
			for i := 0; i < shape.Gndim; i++ {
				rst[i] = shape.NatCoords[i][n]
			}
			shape.Func(shape.S, shape.dSdR, rst[0], rst[1], rst[2], false)
			utl.Pforan("S = %v\n", shape.S)
			for m := 0; m < shape.Nverts; m++ {
				if n == m {
					errS += math.Abs(shape.S[m] - 1.0)
				} else {
					errS += math.Abs(shape.S[m])
				}
			}
		}
		if errS > tol {
			tst.Errorf("%s failed with err = %g\n", name, errS)
			return
		}

		// check dSdR
		tol = 1e-14
		h := 1.0e-1
		S_temp := make([]float64, shape.Nverts)
		if name == "lin5" || name == "tri15" || name == "lin4" || name == "tri10" || name == "qua12" || name == "qua16" {
			tol = 1.0e-10
		}
		for n := 0; n < shape.Nverts; n++ {
			rst := []float64{0, 0, 0}
			for i := 0; i < shape.Gndim; i++ {
				rst[i] = shape.NatCoords[i][n]
			}
			// analytical
			shape.Func(shape.S, shape.dSdR, rst[0], rst[1], rst[2], true)
			// numerical
			for i := 0; i < shape.Gndim; i++ {
				dSndRi, _ := num.DerivCentral(func(x float64, args ...interface{}) (Sn float64) {
					rst_temp := []float64{rst[0], rst[1], rst[2]}
					rst_temp[i] = x
					shape.Func(S_temp, nil, rst_temp[0], rst_temp[1], rst_temp[2], false)
					Sn = S_temp[n]
					return
				}, rst[i], h)
				utl.Pfgrey2("  dS%ddR%d @ [% 4.1f % 4.1f % 4.1f] = %v (num: %v)\n", n, i, rst[0], rst[1], rst[2], shape.dSdR[n][i], dSndRi)
				tol2 := tol
				if name == "tri15" && n == 11 && i == 1 {
					tol2 = 1.0e-9
				}
				if math.Abs(shape.dSdR[n][i]-dSndRi) > tol2 {
					tst.Errorf("%s dS%ddR%d failed with err = %g\n", name, n, i, math.Abs(shape.dSdR[n][i]-dSndRi))
					return
				}
				//utl.CheckScalar(tst, fmt.Sprintf("dS%ddR%d", n, i), tol2, dSdR[n][i], dSndRi)
			}
		}

		// check face vertices
		tol = 1e-17
		errS = 0.0
		if name == "tri10" {
			tol = 1e-14
		}

		nfaces := len(shape.FaceLocalV)
		if nfaces == 0 {
			continue
		}
		for k := 0; k < nfaces; k++ {
			for n := range shape.FaceLocalV[k] {
				rst := []float64{0, 0, 0}
				for i := 0; i < shape.Gndim; i++ {
					rst[i] = shape.NatCoords[i][n]
				}
				shape.Func(shape.S, shape.dSdR, rst[0], rst[1], rst[2], false)
				utl.Pforan("S = %v\n", shape.S)
				for m := range shape.FaceLocalV[k] {
					if n == m {
						errS += math.Abs(shape.S[m] - 1.0)
					} else {
						errS += math.Abs(shape.S[m])
					}
				}
			}
		}
		utl.Pforan("%g\n", errS)
		if errS > tol {
			tst.Errorf("%s failed with err = %g\n", name, errS)
			return
		}

		utl.PfGreen("OK\n")
	}
}
