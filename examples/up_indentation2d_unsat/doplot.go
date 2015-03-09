// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// input data
	simfn := "a-coarse-elast-d2-q9"
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}

	// start analysis process
	out.Start(simfn, 0, 0)

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

	out.Plot(pc_a, "sl", "a", plt.Fmt{M: "o"}, -1)
	out.Csplot.Xlbl = "$p_c$"

	// show
	out.Draw("", "", true)
}
