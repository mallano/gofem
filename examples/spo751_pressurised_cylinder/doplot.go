// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/ana"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// finalise analysis process and catch errors
	defer out.End()
	var sol ana.Hill
	sol.Init(fun.Prms{
		&fun.Prm{N: "a", V: 100},
		&fun.Prm{N: "b", V: 200},
		&fun.Prm{N: "E", V: 210},
		&fun.Prm{N: "ν", V: 0.3},
		&fun.Prm{N: "σy", V: 0.24},
		&fun.Prm{N: "P", V: 0.2},
	})

	// start analysis process
	datadir := "$GOPATH/src/github.com/cpmech/gofem/examples/spo751_pressurised_cylinder/"
	out.Start(datadir+"spo751.sim", 0, 0)

	// define entities
	out.Define("A", out.N{0})
	out.Define("B", out.N{20})
	out.Define("bottom", out.Along{{100, 0}, {101, 0}})

	// load results
	out.LoadResults(nil)

	out.Splot("Pressure at inner and outer face")
	out.Plt("ux", "t", "A", plt.FmtS{"b.-"}, -1)
	out.Plt("ux", "t", "B", plt.FmtS{"b.-"}, -1)

	//out.Splot("radial stress")
	//X, _, _ := out.GetXYZ("ux", "bottom")
	//SR := make([]float64, len(X))
	//ST := make([]float64, len(X))
	//c := sol.Getc()
	//sol.Sig(c, X, SR, ST)

	//out.Plt(X, ST, "bottom", plt.FmtS{"b.-"}, -1)

	// show
	out.Draw("", "", true)
}
