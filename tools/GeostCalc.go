// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"flag"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/io"
)

func main() {

	// input data
	simfile := "simfile.sim"
	zmin := 0.0
	zmax := 3.0
	npts := 11

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfile = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		zmin = io.Atof(flag.Arg(1))
	}
	if len(flag.Args()) > 2 {
		zmax = io.Atob(flag.Arg(2))
	}
	if len(flag.Args()) > 3 {
		npts = io.Atoi(flag.Arg(3))
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfile  = %30s // simulation filename\n", simfile)
	io.Pf("  zmin     = %30s // min elevation\n", zmin)
	io.Pf("  zmax     = %30v // max elevation\n", zmax)
	io.Pf("  npts     = %30v // number of points\n", npts)
	io.Pf("\n")

	// sim file
	sim := inp.ReadSim("", simfile, false)
	if sim == nil {
		io.PfRed("cannot read sim file\n")
		return
	}

	// layer
	var lay fem.GeoLayer
	lay.Zmin = zmin
	lay.Zmax = zmax
	lay.Cl = sim.WaterRho0 / sim.WaterBulk
	//if !lay.ReadPorousParameters(sim.Regions[0],
	// TODO

}
