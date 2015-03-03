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
	out.Start("coarse-elast-d2-q9.sim", 0, 0)

	// define entities
	out.Define("a", out.P{{18, 8}})

	// load results
	out.LoadResults(nil)

	pg_a := out.GetRes("pg", "a", -1)
	pl_a := out.GetRes("pl", "a", -1)
	pc_a := make([]float64, len(pg_a))
	for i, _ := range pg_a {
		pc_a[i] = pg_a[i] - pl_a[i]
	}

	//out.Plt("t", "sl", "a", plt.FmtS{""}, -1)
	out.Plt(pc_a, "sl", "a", plt.FmtS{""}, -1)

	// show
	out.Draw("", "", true)
}
