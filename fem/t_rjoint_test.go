// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_rjoint01a(tst *testing.T) {

	defer func() {
		if err := recover(); err != nil {
			io.Pf("ERROR: %v\n", err)
		}
	}()

	//verbose()
	chk.PrintTitle("rjoint01a")

	// run simulation
	if !Start("data/rjoint01a.sim", true, chk.Verbose) {
		tst.Errorf("Start failed\n")
		return
	}

	// make sure to flush log
	defer End()

	// run simulation
	if !Run() {
		tst.Errorf("Run failed\n")
		return
	}
}
