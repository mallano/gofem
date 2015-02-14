// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"math"
	"testing"

	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// Check checks derivatives
func Check(tst *testing.T, mdl Model, pc0, sl0, pcf float64, npts int, tolCc, tolD1a, tolD1b, tolD2a, tolD2b float64, verbose bool, pcSkip []float64, tolSkip float64, doplot bool) {

	// nonrate model
	nr_mdl, is_nonrate := mdl.(Nonrate)
	utl.Pforan("is_nonrate = %v\n", is_nonrate)

	// for all pc stations
	Pc := utl.LinSpace(pc0, pcf, npts)
	Sl := make([]float64, npts)
	Sl[0] = sl0
	var err error
	for i := 1; i < npts; i++ {

		// update and plot
		Sl[i], err = Update(mdl, Pc[i-1], Sl[i-1], Pc[i]-Pc[i-1])
		if err != nil {
			tst.Errorf("Update failed: %v\n", err)
			return
		}
		if doplot {
			plt.PlotOne(Pc[i], Sl[i], "'ko'")
		}

		// skip point on checking of derivatives
		if doskip(Pc[i], pcSkip, tolSkip) {
			continue
		}

		// wetting flag
		wet := Pc[i]-Pc[i-1] < 0

		// check Cc = dsl/dpc
		utl.Pforan("\npc=%g, sl=%g, wetting=%v\n", Pc[i], Sl[i], wet)
		if is_nonrate {

			// analytical Cc
			Cc_ana, err := mdl.Cc(Pc[i], Sl[i], wet)
			if err != nil {
				tst.Errorf("Cc failed: %v\n", err)
				return
			}

			// numerical Cc
			Cc_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
				return nr_mdl.Sl(x)
			}, Pc[i], 1e-3)
			utl.CheckAnaNum(tst, "Cc = ∂sl/∂pc    ", tolCc, Cc_ana, Cc_num, verbose)
		}

		// compute all derivatives
		L, Lx, J, Jx, Jy, err := mdl.Derivs(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("Derivs failed: %v\n", err)
			return
		}
		L_ana_A := L
		L_ana_B, err := mdl.L(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("L failed: %v\n", err)
			return
		}
		Lx_ana := Lx
		Jx_ana := Jx
		Jy_ana := Jy
		J_ana_A := J
		J_ana_B, err := mdl.J(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("J failed: %v\n", err)
			return
		}

		// numerical L = ∂Cc/∂pc
		L_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Cctmp, _ := mdl.Cc(x, Sl[i], wet)
			return Cctmp
		}, Pc[i], 1e-3)
		utl.CheckAnaNum(tst, "L  = ∂Cc/∂pc    ", tolD1a, L_ana_A, L_num, verbose)

		// numerical Lx := ∂²Cc/∂pc²
		Lx_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ltmp, _, _, _, _, _ := mdl.Derivs(x, Sl[i], wet)
			return Ltmp
		}, Pc[i], 1e-3)
		utl.CheckAnaNum(tst, "Lx = ∂²Cc/∂pc²  ", tolD2a, Lx_ana, Lx_num, verbose)

		// numerical J := ∂Cc/∂sl (version A)
		J_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ccval, _ := mdl.Cc(Pc[i], x, wet)
			return Ccval
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "J  = ∂Cc/∂sl    ", tolD1b, J_ana_A, J_num, verbose)

		// numerical Jx := ∂²Cc/(∂pc ∂sl)
		Jx_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ltmp, _, _, _, _, _ := mdl.Derivs(Pc[i], x, wet)
			return Ltmp
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "Jx = ∂²Cc/∂pc∂sl", tolD2b, Jx_ana, Jx_num, verbose)

		// numerical Jy := ∂²Cc/∂sl²
		Jy_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Jtmp, _ := mdl.J(Pc[i], x, wet)
			return Jtmp
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "Jy = ∂²Cc/∂sl²  ", tolD2b, Jy_ana, Jy_num, verbose)

		// check A and B derivatives
		utl.CheckScalar(tst, "L_A == L_B", 1e-17, L_ana_A, L_ana_B)
		utl.CheckScalar(tst, "J_A == J_B", 1e-17, J_ana_A, J_ana_B)
	}
}

// doskip analyse whether a point should be skip or not
func doskip(x float64, xskip []float64, tol float64) bool {
	for _, v := range xskip {
		if math.Abs(x-v) < tol {
			return true
		}
	}
	return false
}
