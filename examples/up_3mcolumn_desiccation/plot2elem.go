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
	simfn := "onepulse-qua9co.sim"
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
	out.Define("A C E", out.N{-5, -3, -1})
	out.Define("a b", out.P{{1, 8}, {0, 8}})

	// load results
	out.LoadResults(nil)

	// styles
	me := 10
	Sn := []plt.Fmt{
		plt.Fmt{C: "b", M: "*", Me: me},
		plt.Fmt{C: "m", M: "x", Me: me},
		plt.Fmt{C: "r", M: "^", Me: me},
	}
	Se := []plt.Fmt{
		plt.Fmt{C: "b", M: "*", Me: me},
		plt.Fmt{C: "g", M: "o", Me: me},
	}

	// pl
	out.Splot("liquid pressure")
	for i, l := range []string{"A", "C", "E"} {
		out.Plot("t", "pl", l, Sn[i], -1)
	}

	// uy
	out.Splot("displacements")
	for i, l := range []string{"A", "C", "E"} {
		out.Plot("t", "uy", l, Sn[i], -1)
	}

	out.Splot("liquid saturation")
	for i, l := range []string{"a", "b"} {
		out.Plot("t", "sl", l, Se[i], -1)
	}

	out.Splot("stresses")
	for i, l := range []string{"a", "b"} {
		out.Plot("t", "sy", l, Se[i], -1)
	}

	// show
	out.Draw("", "", true)
}
