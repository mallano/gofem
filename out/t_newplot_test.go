// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gosl/chk"
)

func Test_plot(tst *testing.T) {

	verbose()

	chk.PrintTitle("plot")
	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"

	Read(datadir+"p02.sim", 0, 0)

	SetPoint(1, "A")
	SetPointAt([]float64{0, 1}, "B")
	SetAlong([]float64{0, 0}, []float64{0, 10}, "left")

	Apply()

	plA := Get("pl", "A")
	//io.Pforan("plA = %v\n", plA)

	//os.Exit(0)

	SetTitle("My fancy plot")
	Plot("t", "pl", "B", "'b.-'", "")
	Plot("t", plA, "A", "'r.-'", "")

	Subplot()
	Plot("pl", "pl", "A", "'ko-'", "")

	Subplot()
	Plot("pl", "y", "left", "'bo-'", "")

	Subplot()
	Plot("y", "pl", "left", "'bo-'", "")

	Show()

}
