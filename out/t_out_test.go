// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_out01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("out01")

	var styles []string
	var labels []string
	var extra string

	stageIdx := 0
	regionIdx := 0

	ti := 0.0
	tf := -1.0

	defer With("data/p01.sim", "data", stageIdx, regionIdx)()
	Tplot("pl", &At{2.5, 10, 0}, ti, tf, styles, labels, extra)
	Tplot("pl", IdsOrTags{-5}, ti, tf, styles, labels, extra)
	Show()
}
