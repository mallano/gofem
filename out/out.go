// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package out implements FE simulation output handling for analyses and plotting
package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/io"
)

// Global variables
var (

	// constants
	TolC     float64 // tolerance to compare x-y-z coordinates
	TolT     float64 // tolerance to compare times
	Ndiv     int     // bins n-division
	SubpNrow int     // override subplot configuration parameters
	SubpNcol int     // override subplot configuration parameters

	// data
	Sum     *fem.Summary     // summary of results
	Dom     *fem.Domain      // FE domain
	Ipoints []*fem.OutIpData // all integration points
	Cid2ips [][]int          // [ncells][nip] maps cell id to index in Ipoints
	NodBins gm.Bins          // bins for nodes
	IpsBins gm.Bins          // bins for integration points

	// results
	R ResultsMap // maps labels => points
	I []int      // selected output indices
	T []float64  // selected output times

	// suplot data
	Spd map[string][]int // [nkeys] subplot data
)

// End must be called and the end to flush log file
func End() {
	if err := recover(); err != nil {
		io.PfRed("ERROR: %v\n", err)
	} else {
		fem.End()
	}
}

// Start starts handling of results given a simulation input file
func Start(simfnpath string, stageIdx, regionIdx int) (startisok bool) {

	// constants
	TolC = 1e-8
	TolT = 1e-3
	Ndiv = 20
	SubpNrow = 0
	SubpNcol = 0

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		chk.Panic("cannot start analysis process with simfnpath=%q\n", simfnpath)
	}

	// read summary
	Sum = fem.ReadSum()
	if Sum == nil {
		chk.Panic("cannot read summary file for simulation=%q\n", simfnpath)
	}

	// allocate domain
	Dom = fem.NewDomain(fem.Global.Sim.Regions[regionIdx])
	if !Dom.SetStage(stageIdx, fem.Global.Sim.Stages[stageIdx]) {
		chk.Panic("cannot allocate domain\n")
	}

	// clear previous data
	Ipoints = make([]*fem.OutIpData, 0)
	Cid2ips = make([][]int, len(Dom.Msh.Cells))
	R = make(map[string]Points)
	I = make([]int, 0)
	T = make([]float64, 0)
	Spd = make(map[string][]int)

	// bins
	m := Dom.Msh
	xi := []float64{m.Xmin, m.Ymin}
	xf := []float64{m.Xmax, m.Ymax}
	if m.Ndim == 3 {
		xi = append(xi, m.Zmin)
		xf = append(xf, m.Zmax)
	}
	NodBins.Init(xi, xf, Ndiv)
	IpsBins.Init(xi, xf, Ndiv)

	// add nodes to bins
	for activeId, nod := range Dom.Nodes {
		err := NodBins.Append(nod.Vert.C, activeId)
		if err != nil {
			return
		}
	}

	// add integration points to slice of ips and to bins
	for cid, ele := range Dom.Cid2elem {
		if ele == nil {
			continue
		}
		dat := ele.OutIpsData()
		nip := len(dat)
		ids := make([]int, nip)
		for i, d := range dat {
			id := len(Ipoints)
			ids[i] = id
			Ipoints = append(Ipoints, d)
			IpsBins.Append(d.X, id)
		}
		Cid2ips[cid] = ids
	}

	// success
	return true
}
