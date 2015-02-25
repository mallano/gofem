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
	"github.com/cpmech/gosl/utl"
)

// constants
var (
	TolC = 1e-8 // tolerance to compare x-y-z coordinates
	TolT = 1e-3 // tolerance to compare times
	Ndiv = 20   // bins n-division
)

// Global variables
var (

	// data set by Start
	//Sum       *fem.Summary     // summary of results
	Dom       *fem.Domain      // FE domain
	Ipoints   []*fem.OutIpData // all integration points
	Cid2ips   [][]int          // [ncells][nip] maps cell id to index in Ipoints
	Ipkey2ips map[string][]int // maps ip keys to indices in Ipoints
	Ipkeys    map[string]bool  // all ip keys
	NodBins   gm.Bins          // bins for nodes
	IpsBins   gm.Bins          // bins for integration points

	// auxiliary

	// results loaded by LoadResults
	R ResultsMap // maps labels => points
	I []int      // selected output indices
	T []float64  // selected output times

	// subplots
	Splots []*SplotDat // all subplots
	Csplot *SplotDat   // current subplot
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

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		chk.Panic("cannot start analysis process with simfnpath=%q\n", simfnpath)
	}

	// read summary
	if !fem.Global.Sum.Read() {
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
	Ipkey2ips = make(map[string][]int)
	Ipkeys = make(map[string]bool)
	R = make(map[string]Points)
	I = make([]int, 0)
	T = make([]float64, 0)
	Splots = make([]*SplotDat, 0)

	// bins
	m := Dom.Msh
	xi := []float64{m.Xmin, m.Ymin}
	xf := []float64{m.Xmax, m.Ymax}
	if m.Ndim == 3 {
		xi = append(xi, m.Zmin)
		xf = append(xf, m.Zmax)
	}
	err := NodBins.Init(xi, xf, Ndiv)
	if err != nil {
		chk.Panic("cannot initialise bins for nodes: %v", err)
	}
	err = IpsBins.Init(xi, xf, Ndiv)
	if err != nil {
		chk.Panic("cannot initialise bins for integration points: %v", err)
	}

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
			err = IpsBins.Append(d.X, id)
			if err != nil {
				chk.Panic("cannot append to bins of integration points: %v", err)
			}
			for key, _ := range d.V {
				utl.StrIntsMapAppend(&Ipkey2ips, key, id)
				Ipkeys[key] = true
			}
		}
		Cid2ips[cid] = ids
	}

	// success
	return true
}
