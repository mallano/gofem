// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"

	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/io"
)

func main() {

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
		convSymb = io.Atob(flag.Arg(2))
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  matOld   = %20s // old material filename\n", matOld)
	io.Pf("  matNew   = %20s // new material filename\n", matNew)
	io.Pf("  convSymb = %20v // do convert symbols\n", convSymb)
	io.Pf("\n")

	// convert old => new
	inp.MatfileOld2New("", matNew, matOld, convSymb)
	io.Pfblue2("conversion successful\n")
}
