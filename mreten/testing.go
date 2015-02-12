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
func Check(tst *testing.T, mdl Model, pc0, sl0, pcf float64, npts int, tolCc, tolD1, tolD2 float64, verbose bool, pcSkip []float64, tolSkip float64, doplot bool) {

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

		// skip point
		if doskip(Pc[i], pcSkip, tolSkip) {
			continue
		}

		// update and plot
		Sl[i], err = Update(mdl, Pc[i-1], Sl[i-1], Pc[i]-Pc[i-1])
		if err != nil {
			tst.Errorf("Update failed: %v\n", err)
			return
		}
		if doplot {
			plt.PlotOne(Pc[i], Sl[i], "'ko'")
		}

		// wetting flag
		wet := Pc[i]-Pc[i-1] < 0

		// check Cc = dsl/dpc
		utl.Pf("\n")
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

			// comparison
			utl.CheckAnaNum(tst, utl.Sf("Cc        @ %+.3f", Pc[i]), tolCc, Cc_ana, Cc_num, verbose)
		}

		// compute all derivatives
		mdl.Derivs(Pc[i], Sl[i], wet)
		DCcDpc_ana := D.DCcDpc
		D2CcDpc2_ana := D.D2CcDpc2

		// numerical DCcDpc = ∂Cc/∂pc
		DCcDpc_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			Ccval, _ := mdl.Cc(x, Sl[i], wet)
			return Ccval
		}, Pc[i], 1e-3)

		// comparison
		utl.CheckAnaNum(tst, utl.Sf("∂Cc/∂pc   @ %+.3f", Pc[i]), tolD1, DCcDpc_ana, DCcDpc_num, verbose)

		// numerical D2CcDpc2 := ∂²Cc/∂pc²
		D2CcDpc2_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			mdl.Derivs(x, Sl[i], wet)
			return D.DCcDpc
		}, Pc[i], 1e-3)

		// comparison
		utl.CheckAnaNum(tst, utl.Sf("∂²Cc/∂pc² @ %+.3f", Pc[i]), tolD2, D2CcDpc2_ana, D2CcDpc2_num, verbose)

		// if is nonrate, skip additional derivatives checks
		if is_nonrate {
			continue
		}

		// analytical derivatives
		DCcDsl_ana := D.DCcDsl

		// numerical DCcDsl := ∂Cc/∂sl
		DCcDsl_num, _ := num.DerivCentral(func(x float64, args ...interface{}) float64 {
			mdl.Derivs(x, Sl[i], wet)
			return D.DCcDsl
		}, Pc[i], 1e-3)

		// comparison
		utl.CheckAnaNum(tst, utl.Sf("∂Cc/∂sl   @ %+.3f", Pc[i]), tolD2, DCcDsl_ana, DCcDsl_num, verbose)
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
