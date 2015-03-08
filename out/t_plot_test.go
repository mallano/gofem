// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func Test_plot01(tst *testing.T) {

	// finalise analysis process and catch errors
	defer func() {
		if err := recover(); err != nil {
			tst.Fail()
			io.PfRed("ERROR: %v\n", err)
		} else {
			fem.End()
		}
	}()

	// test title
	//verbose()
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

	// linewidth and markersize; -1 => use default
	lw, ms := -1.0, -1.0

	Splot("liquid pressure")
	Plot("t", "pl", "B", plt.Fmt{"b", ".", "-", lw, ms, ""}, -1)
	Plot("t", plA, "A", plt.Fmt{"r", ".", "-", lw, ms, ""}, -1)

	Splot("")
	Plot("pl", "pl", "A", plt.Fmt{"k", "o", "-", lw, ms, ""}, -1)

	Splot("")
	Plot("pl", "y", "left", plt.Fmt{"b", "o", "-", lw, ms, ""}, 0)
	Plot("pl", "y", "left", plt.Fmt{"g", "o", "-", lw, ms, ""}, -1)

	Splot("")
	io.Pforan("T = %v\n", T)
	last := len(T) - 1
	Plot("y", "pl", "left", plt.Fmt{"b", "o", "-", lw, ms, io.Sf("t=%g", T[0])}, 0)
	Plot("y", "pl", "left", plt.Fmt{"m", "*", "-", 2, ms, io.Sf("t=%g", T[last])}, -1)

	//Draw("", "", true)
	Draw("", "", false)
}
