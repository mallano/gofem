// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
)

func main() {

	// catch errors
	var tst testing.T
	defer func() {
		if mpi.Rank() == 0 {
			if err := recover(); err != nil {
				io.PfRed("ERROR: %v\n", err)
			}
			if tst.Failed() {
				io.PfRed("test failed\n")
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// start global variables and log
	if !fem.Start("data/p01.sim", true, true) {
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
}
