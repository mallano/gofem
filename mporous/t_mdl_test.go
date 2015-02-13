// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_mdl01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("mdl01")

	var mdl Model
	var cnd mconduct.M1
	var lrm mreten.RefM1
	mdl_prms := mdl.GetPrms(true)
	cnd_prms := cnd.GetPrms(true)
	lrm_prms := lrm.GetPrms(true)
	var prms fun.Prms
	for _, p := range mdl_prms {
		prms = append(prms, p)
	}
	for _, p := range cnd_prms {
		prms = append(prms, p)
	}
	for _, p := range lrm_prms {
		prms = append(prms, p)
	}
	utl.Pfcyan("prms = [\n%v\n", prms)

	mdl.Init(prms, &cnd, &lrm)
	utl.Pforan("mdl = %v\n", mdl)

	var sta StateLG
	sta.Sl = 1.0
	Δpl := -20.0
	err := mdl.Update(&sta, Δpl, 0)
	if err != nil {
		tst.Errorf("Updated failed: %v\n", err)
		return
	}
	utl.Pforan("sl = %v\n", sta.Sl)

	pc0 := -5.0
	sl0 := 1.0
	pcf := 20.0
	nptsA := 41

	doplot := true
	if doplot {
		plt.Reset()
	}
	mreten.Plot(mdl.Lrm, pc0, sl0, pcf, nptsA, "'b.-'", "'r+-'", "ref-m1_drying")
	if doplot {
		mreten.PlotEnd(true)
	}
}
