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

	utl.Tsilent = false
	defer func() {
		if err := recover(); err != nil {
			if mpi.IsOn() {
				if mpi.Rank() == 0 {
					utl.PfRed("\n\nSome error has happened: %v\n", err)
				}
			} else {
				utl.PfRed("\n\nSome error has happened: %v\n", err)
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// run simulation
	fem.Start("data/bh16.sim", true, true)
	defer fem.End()
	fem.Run()

	// check
	skipK := true
	tolK := 1e-12
	tolu := 1e-15
	tols := 1e-12
	verb := true
	var tst testing.T
	fem.TestingCompareResultsU(&tst, "data/bh16.sim", "cmp/bh16.cmp", tolK, tolu, tols, skipK, verb)
}
