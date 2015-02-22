// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gosl/chk"
)

func Test_plot01(tst *testing.T) {

	// finalise analysis process and catch errors
	defer End()

	// test title
	verbose()
	chk.PrintTitle("plot01")

	// start analysis process
	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"
	Start(datadir+"p02.sim", 0, 0)

	// define entities
	Define("A", N{1})
	Define("B", At{0, 1})
	Define("left", Along{{0, 0}, {0, 10}})

	// load results
	LoadResults(nil)

	plA := GetRes("pl", "A", 0)

	Splot("liquid pressure")
	Plt("t", "pl", "B", "'b.-'", -1)
	Plt("t", plA, "A", "'r.-'", -1)

	Splot("")
	Plt("pl", "pl", "A", "'ko-'", -1)

	Splot("")
	Plt("pl", "y", "left", "'bo-'", 0)
	Plt("pl", "y", "left", "'bo-'", -1)

	Splot("")
	Plt("y", "pl", "left", "'bo-'", 0)
	Plt("y", "pl", "left", "'bo-'", -1)

	Draw("", "", true)
}
