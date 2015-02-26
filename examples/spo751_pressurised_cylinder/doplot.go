// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import "github.com/cpmech/gofem/out"

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// start analysis process
	datadir := "$GOPATH/src/github.com/cpmech/gofem/examples/spo751_pressurised_cylinder/"
	out.Start(datadir+"spo751.sim", 0, 0)

	// define entities
	out.Define("A", out.N{20})

	// load results
	out.LoadResults(nil)

	out.Splot("Internal pressurised cylinder")
	out.Plt("ux", "t", "A", "'b.-'", -1)

	out.Draw("", "", false)
}
