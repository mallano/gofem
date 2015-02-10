// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/utl"
)

func main() {

	// error handler
	utl.Tsilent = false
	defer func() {
		if err := recover(); err != nil {
			utl.PfRed("Some error has happened: %v\n", err)
		}
	}()

	// input data
	matOld := "matOld.mat"
	matNew := "matNew.mat"
	convSymb := true
	flag.Parse()
	if len(flag.Args()) > 0 {
		matOld = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		matNew = flag.Arg(1)
	}
	if len(flag.Args()) > 2 {
		convSymb = utl.Atob(flag.Arg(2))
	}

	// print input data
	utl.Pf("\nInput data\n")
	utl.Pf("==========\n")
	utl.Pf("  matOld   = %20s // old material filename\n", matOld)
	utl.Pf("  matNew   = %20s // new material filename\n", matNew)
	utl.Pf("  convSymb = %20v // do convert symbols\n", convSymb)
	utl.Pf("\n")

	// convert old => new
	inp.MatfileOld2New("", matNew, matOld, convSymb)
	utl.Pfblue2("conversion successful\n")
}
