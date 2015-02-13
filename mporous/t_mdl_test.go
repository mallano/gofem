// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

import (
	"testing"

	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/fun"
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

	var lrm mreten.RefM1
	lrm_prms := lrm.GetPrms()
	cnd_prms := GetLGcndPrms()
	mdl_prms := GetModelPrms()
	utl.Pforan("cnd_prms = [\n%v\n", cnd_prms)
	utl.Pforan("lrm_prms = [\n%v\n", lrm_prms)
	utl.Pforan("mdl_prms = [\n%v\n", mdl_prms)

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
}
