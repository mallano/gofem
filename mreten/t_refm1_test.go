// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"testing"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_refm1a(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("refm1a")

	model := GetModel("testsim", "mat1", "ref-m1", false)
	model.Init(model.GetPrms())
	mdl := model.(RateType)
	utl.Pforan("mdl = %v\n", mdl)

	// plot
	if false {

		var der Derivs
		pc0 := -5.0
		sl0 := 1.0
		Δpc := 25.0
		wetting := false
		npts := 31

		// ode solver
		//   x      = [0.0, 1.0]
		//   pc     = pc0 + x * Δpc
		//   y[0]   = sl
		//   f(x,y) = dy/dx = dsl/dpc * dpc/dx = Cc * Δpc
		//   J(x,y) = df/dy = DCcDsl * Δpc
		fcn := func(f []float64, x float64, y []float64, args ...interface{}) (err error) {
			f[0], err = mdl.Cc(pc0+x*Δpc, y[0], wetting)
			f[0] *= Δpc
			return
		}
		jac := func(dfdy *la.Triplet, x float64, y []float64, args ...interface{}) (err error) {
			if dfdy.Max() == 0 {
				dfdy.Init(1, 1, 1)
			}
			err = mdl.Derivs(&der, pc0+x*Δpc, y[0], wetting)
			dfdy.Start()
			dfdy.Put(0, 0, der.DCcDsl)
			return nil
		}

		var odesol ode.ODE
		odesol.Init("Radau5", 1, fcn, jac, nil, nil, true)
		odesol.SetTol(1e-10, 1e-7)

		Pc := make([]float64, npts)
		Sl := make([]float64, npts)
		Pc[0] = pc0
		Sl[0] = sl0
		X := utl.LinSpace(0, 1, npts)
		for i := 1; i < npts; i++ {
			y := []float64{Sl[i-1]}
			err := odesol.Solve(y, X[i-1], X[i], X[i]-X[i-1], false)
			if err != nil {
				tst.Errorf("ode failed: %v\n", err)
				return
			}
			//utl.Pforan("x, xb = %v, %v\n", X[i-1], X[i])
			Pc[i] = pc0 + X[i]*Δpc
			Sl[i] = y[0]
		}
		//utl.Pforan("Pc = %v\n", Pc)
		//utl.Pforan("Sl = %v\n", Sl)

		plt.Plot(Pc, Sl, "'b.-', label='ref-m1', clip_on=0")
		plt.Gll("$p_c$", "$s_{\\ell}$", "")
		plt.Cross()
		plt.Show()
	}
}
