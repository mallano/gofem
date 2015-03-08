// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"encoding/json"

	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
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

	// read comparison results
	rcmp_nod := read_pyfem_rod_data("cmp/pyfem_o2_rod_nod.dat")

	// plot uy along y for selected times
	out.Splot("rod displacements")
	for i, _ := range out.I {
		out.Plot("y", "uy", "rod", plt.Fmt{}, i)
	}
	for _, d := range rcmp_nod {
		plt.Plot(d.Y, d.Uy, "'k+', ms=5")
	}

	// show
	out.Draw("", "", true)
}

// auxiliary /////////////////////////////////////////////////////////////////////////////////////////

type PyfemRodData struct {
	// all
	Time float64
	X    []float64
	Y    []float64
	Z    []float64

	// rod
	Ea []float64 // len(Ea) == Nips
	Sa []float64
	Uy []float64 // len(Uy) == Nnodes

	// joint
	Tau  []float64 // len == Nips
	Sigc []float64
	Ur   []float64 // ω
	W_pa []float64 // ωp
}

func read_pyfem_rod_data(fn string) []PyfemRodData {

	// read file
	b, err := io.ReadFile(fn)
	if err != nil {
		chk.Panic("%v\n", err.Error())
	}

	// decode
	var dat []PyfemRodData
	err = json.Unmarshal(b, &dat)
	if err != nil {
		chk.Panic("%v\n", err.Error())
	}
	return dat
}
