// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import "github.com/cpmech/gofem/out"

func main() {

	// finalise analysis process and catch errors upon exit
	defer out.End()

	// start analysis process
	out.Start("spo754.sim", 0, 0)

	// all nodes
	out.Define("all nodes", out.AllNodes())

	// load results
	out.LoadResults([]float64{-1})

	// plot
	out.Contour("uy", -1)
}
