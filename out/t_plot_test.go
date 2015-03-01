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

	Splot("liquid pressure")
	Plt("t", "pl", "B", plt.FmtS{"b.-"}, -1)
	Plt("t", plA, "A", plt.FmtS{"r.", "-"}, -1)

	Splot("")
	Plt("pl", "pl", "A", plt.FmtS{"k", "o", "-"}, -1)

	Splot("")
	Plt("pl", "y", "left", plt.FmtS{"bo-"}, 0)
	Plt("pl", "y", "left", plt.FmtS{"go-"}, -1)

	Splot("")
	Plt("y", "pl", "left", plt.FmtS{"b", "o", "-"}, 0)
	Plt("y", "pl", "left", plt.FmtL{C: "m", M: "*", Ls: "-", Lw: 2}, -1)

	//Draw("", "", true)
	Draw("", "", false)
}
