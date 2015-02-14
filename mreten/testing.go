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
	utl.Pforan("is_nonrate = %v\n", nr_mdl)
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
			utl.CheckAnaNum(tst, "Cc            ", tolCc, Cc_ana, Cc_num, verbose)
		}

		// compute all derivatives
		err = mdl.Derivs(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("Derivs failed: %v\n", err)
			return
		}
		DCcDpc_ana := D.DCcDpc
		D2CcDpc2_ana := D.D2CcDpc2
		D2CcDsl2_ana := D.D2CcDsl2
		D2CcDpcDsl_ana := D.D2CcDpcDsl
		A_DCcDsl_ana := D.DCcDsl
		B_DCcDsl_ana, err := mdl.DCcDsl(Pc[i], Sl[i], wet)
		if err != nil {
			tst.Errorf("DCcDsl failed: %v\n", err)
			return
		}

		// check A and B derivatives
		utl.CheckScalar(tst, "DCcDsl A==B", 1e-17, A_DCcDsl_ana, B_DCcDsl_ana)

		// numerical DCcDpc = ∂Cc/∂pc
		DCcDpc_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ccval, _ := mdl.Cc(x, Sl[i], wet)
			return Ccval
		}, Pc[i], 1e-3)
		utl.CheckAnaNum(tst, "∂Cc/∂pc       ", tolD1a, DCcDpc_ana, DCcDpc_num, verbose)

		// numerical D2CcDpc2 := ∂²Cc/∂pc²
		D2CcDpc2_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			mdl.Derivs(x, Sl[i], wet)
			return D.DCcDpc
		}, Pc[i], 1e-3)
		utl.CheckAnaNum(tst, "∂²Cc/∂pc²     ", tolD2a, D2CcDpc2_ana, D2CcDpc2_num, verbose)

		// numerical DCcDsl := ∂Cc/∂sl
		A_DCcDsl_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ccval, _ := mdl.Cc(Pc[i], x, wet)
			return Ccval
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "∂Cc/∂sl (A)   ", tolD1b, A_DCcDsl_ana, A_DCcDsl_num, verbose)

		// numerical DCcDsl := ∂Cc/∂sl (version B)
		B_DCcDsl_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ccval, _ := mdl.Cc(Pc[i], x, wet)
			return Ccval
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "∂Cc/∂sl (B)   ", tolD1b, B_DCcDsl_ana, B_DCcDsl_num, verbose)

		// numerical D2CcDsl2 := ∂²Cc/∂sl²
		D2CcDsl2_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			mdl.Derivs(Pc[i], x, wet)
			return D.DCcDsl
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "∂²Cc/∂sl²     ", tolD2b, D2CcDsl2_ana, D2CcDsl2_num, verbose)

		// numerical D2CcDpcDsl := ∂²Cc/(∂pc ∂sl)
		D2CcDpcDsl_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			mdl.Derivs(Pc[i], x, wet)
			return D.DCcDpc
		}, Sl[i], 1e-3)
		utl.CheckAnaNum(tst, "∂²Cc/(∂pc ∂sl)", tolD2b, D2CcDpcDsl_ana, D2CcDpcDsl_num, verbose)
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
