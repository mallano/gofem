// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// catch errors
	utl.Tsilent = false
	var tst testing.T
	defer func() {
		if mpi.Rank() == 0 {
			if err := recover(); err != nil {
				utl.PfRed("ERROR: %v\n", err)
			}
			if tst.Failed() {
				utl.PfRed("test failed\n")
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// start global variables and log
	if !fem.Start("data/spo751.sim", true, true) {
		tst.Error("Start failed\n")
		return
	}

	// make sure to flush log
	defer fem.End()

	// run simulation
	if !fem.Run() {
		tst.Error("Run failed\n")
		return
	}

	// check
	skipK := true
	tolK := 1e-17
	tolu := 1e-12
	tols := 1e-14
	verb := false
	fem.TestingCompareResultsU(&tst, "data/spo751.sim", "cmp/spo751.cmp", tolK, tolu, tols, skipK, verb)
}
