// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_p01(tst *testing.T) {

	/* this tests simulates seepage flow along a column
	 * by reducing the initial hydrostatic pressure at
	 * at the bottom of the column
	 *
	 *       8 o----o----o 9 (-5)
	 *       |   14    |
	 *       |         |
	 *    21 o    o    o 22 (-6)
	 *       |   26    |
	 *       |         |
	 *     6 o----o----o 7 (-4)
	 *       |   13    |
	 *       |         |
	 *    19 |    o    o 20 (-6)
	 *       |   25    |
	 *       |         |
	 *     4 o----o----o 5 (-3)
	 *       |   12    |
	 *       |         |
	 *    17 o    o    o 18 (-6)
	 *       |   24    |
	 *       |         |
	 *     2 o----o----o 3 (-2)
	 *       |   11    |
	 *       |         |
	 *    15 o    o    o 16 (-6)
	 *       |   23    |
	 *       |         |
	 *     0 o----o----o 1 (-1)
	 *           10
	 */

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	utl.Tsilent = false
	utl.TTitle("p01")

	// run simulation
	if !Start("data/p01.sim", true, !utl.Tsilent) {
		tst.Errorf("test failed\n")
		return
	}

	// make sure to flush log
	defer End()

	// run simulation
	if !Run() {
		tst.Errorf("test failed\n")
		return
	}
}
