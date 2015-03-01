// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// start analysis process
	out.Start("onepulse-qua9co.sim", 0, 0)

	// define entities
	out.Define("A B C D E", out.N{-5, -4, -3, -2, -1})
	//out.Define("A", out.N{-5})

	// load results
	out.LoadResults(nil)

	out.Splot("liquid pressure")
	out.Plt("t", "pl", "A", plt.FmtS{"b*-"}, -1)
	out.Plt("t", "pl", "B", plt.FmtS{"go-"}, -1)
	out.Plt("t", "pl", "C", plt.FmtS{"mx-"}, -1)
	out.Plt("t", "pl", "D", plt.FmtL{M: "+", C: "orange", Ls: "-"}, -1)
	out.Plt("t", "pl", "E", plt.FmtS{"r^-"}, -1)

	// show
	out.Draw("", "", true)
}
