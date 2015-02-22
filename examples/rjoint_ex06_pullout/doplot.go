// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/plt"
)

func main() {

	// finalise analysis process and catch errors upon exit
	defer out.End()

	// start analysis process
	out.Start("o2elast.sim", 0, 0)

	// define entities
	out.Define("rod", out.Along{{0.05, 0.2, 0.05}, {0.05, 0.6, 0.05}})

	// load results
	out.LoadResults([]float64{0.2, 0.4, 0.6, 0.8, 0.9, 0.98, 1})

	// distances along rod
	dist := out.GetDist("uy", "rod")

	// get results
	for i, _ := range out.I {
		uy := out.GetX("uy", "rod", i)
		plt.Plot(dist, uy, "")
	}

	// plot
	plt.Show()
}
