// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package out implements FE simulation output handling for analyses and plotting
//  The main structures containing results are:
//   TseriesTimes -- slice of time values; e.g. T = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0}
//   TseriesRres  -- all results for keys, quantities (aka points), and times.
//        TseriesRes[nkeys][npts..?][ntimes]. Example:
//        TseriesRes = [][][]float64{
//          {
//            {100, 99, 98}, // pl @ point A (bottom of column) for 3 time outputs
//            {  0,  0,  0}, // pl @ point B (top of column) for 3 time outputs
//          },
//          {
//            {0,      0,      0}, // uy @ point A for 3 time outputs
//            {0, -0.001, -0.002}, // uy @ point B for 3 time outputs
//          },
//        }
//   Spd -- contains subplot data for further configurations. Example:
//          spd := {"pl" : []int{1, 1, 2}} == "pl" => nrow,ncol,idx
package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// IpDat holds integration point data and serves to aid locating points in space
type IpDat struct {
	Cid int            // id of parent element
	X   []float64      // ip's coordinates
	P   *mporous.State // state @ p-element's ip
	U   *msolid.State  // state @ u-element's ip
}

// Global variables
var (

	// constants
	TolC     float64 // tolerance to compare x-y-z coordinates
	Ndiv     int     // bins n-division
	SubpNrow int     // override subplot configuration parameters
	SubpNcol int     // override subplot configuration parameters

	// data
	Sum     *fem.Summary // summary of results
	Dom     *fem.Domain  // FE domain
	Ipoints []*IpDat     // all integration points
	Cid2ips [][]int      // [ncells][nip] maps cell id to index in Ipoints
	NodBins gm.Bins      // bins for nodes
	IpsBins gm.Bins      // bins for integration points

	// time-series data
	TseriesKeys    []string       // [nkeys] all keys
	TseriesData    []*TseriesDat  // [nkeys] all items
	TseriesTimes   []float64      // [ntimes] time series
	TseriesRes     [][][]float64  // [nkeys][nqts..?][ntimes] all time-series results
	TseriesKey2idx map[string]int // [nkeys] key => index in TseriesKeys or TseriesData

	// suplot data
	Spd map[string][]int // [nkeys] subplot data
)

// End must be called and the end to flush log file
func End() {
	fem.End()
}

// Start starts handling of results given a simulation input file
func Start(simfnpath string, stageIdx, regionIdx int) (startisok bool) {

	// constants
	TolC = 1e-8
	Ndiv = 20
	SubpNrow = 0
	SubpNcol = 0

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		return
	}

	// read summary
	Sum = fem.ReadSum()
	if Sum == nil {
		return
	}

	// allocate domain
	Dom = fem.NewDomain(fem.Global.Sim.Regions[regionIdx])
	if !Dom.SetStage(stageIdx, fem.Global.Sim.Stages[stageIdx]) {
		return
	}

	// clear previous data
	TseriesClear()
	Ipoints = make([]*IpDat, 0)
	Cid2ips = make([][]int, len(Dom.Msh.Cells))
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
	for _, ele := range Dom.Elems {
		var cell *inp.Cell
		var ips []*shp.Ipoint
		var xmat [][]float64
		var pstates []*mporous.State
		var ustates []*msolid.State
		switch e := ele.(type) {
		case *fem.ElemP:
			cell = e.Cell
			ips = e.IpsElem
			xmat = e.X
			pstates = e.States
		case *fem.ElemU:
			cell = e.Cell
			ips = e.IpsElem
			xmat = e.X
			ustates = e.States
		}
		if cell != nil {
			ids := make([]int, len(ips))
			for idx, ip := range ips {
				x := cell.Shp.IpRealCoords(xmat, ip)
				id := len(Ipoints)
				ids[idx] = id
				var pstate *mporous.State
				var ustate *msolid.State
				if pstates != nil {
					pstate = pstates[idx]
				}
				if ustates != nil {
					ustate = ustates[idx]
				}
				Ipoints = append(Ipoints, &IpDat{cell.Id, x, pstate, ustate})
				IpsBins.Append(x, id)
			}
			Cid2ips[cell.Id] = ids
		}
	}

	// success
	return true
}

// Apply applies commands to generate T and R.
func Apply() (err error) {
	TseriesStart()
	TseriesTimes = make([]float64, Sum.NumTidx)
	for tidx := 0; tidx < Sum.NumTidx; tidx++ {
		if !Dom.In(tidx) {
			return chk.Err("Domain.In failed. See log files\n")
		}
		TseriesTimes[tidx] = Dom.Sol.T
		for i, dat := range TseriesData {
			for j, q := range dat.Qts {
				TseriesRes[i][j][tidx] = *q.Value
			}
		}
	}
	return
}

// Show shows plot
//  extra -- is a function to carry out extra configurations
func Show(extra func()) {
	plot_all()
	if extra != nil {
		extra()
	}
	plt.Show()
	return
}

// Save saves plot
func Save(dirout, filename string, extra func()) (err error) {
	plot_all()
	if extra != nil {
		extra()
	}
	plt.SaveD(dirout, filename)
	return
}

// plot_all plots all results
func plot_all() {
	plt.Reset()
	nplots := len(TseriesRes)
	nrow, ncol := utl.BestSquare(nplots)
	if SubpNrow > 0 {
		nrow = SubpNrow
	}
	if SubpNcol > 0 {
		ncol = SubpNcol
	}
	Spd = make(map[string][]int)
	for i := 0; i < nplots; i++ {
		key := TseriesKeys[i]
		Spd[key] = []int{nrow, ncol, i + 1}
		plt.Subplot(nrow, ncol, i+1)
		for j, Y := range TseriesRes[i] {
			sty := TseriesData[i].Sty[j]
			args := sty.GetArgs("clip_on=0")
			plt.Plot(TseriesTimes, Y, args)
		}
		plt.Gll(GetTexLabel("time", ""), GetTexLabel(key, ""), "")
	}
	return
}
