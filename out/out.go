// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// IpDat holds integration point data
type IpDat struct {
	Eid int            // id of parent element
	X   []float64      // ip's coordinates
	P   *mporous.State // state @ p-element's ip
	U   *msolid.State  // state @ u-element's ip
}

// SubpDat maps keys to subplot configuration parameters
//  Example: "pl" => 1,1,2 == row,col,idx
type SubpDat map[string][]int

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
	NodBins gm.Bins      // bins for nodes
	IpsBins gm.Bins      // bins for integration points
)

// End must be called and the end to flush log file
func End() {
	fem.End()
}

// Start starts handling and plotting of results given a simulation input file
func Start(simfnpath string, stageIdx, regionIdx int) (startisok bool) {

	// constants
	TolC = 1e-8
	Ndiv = 20

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
	for activeId, n := range Dom.Nodes {
		err := NodBins.Append(n.Vert.C, activeId)
		if err != nil {
			return
		}
	}

	// add integration points to slice of ips and to bins
	for _, ele := range Dom.Elems {
		switch e := ele.(type) {
		case *fem.ElemP:
			for idx, ip := range e.IpsElem {
				C := e.Cell.Shp.IpRealCoords(e.X, ip)
				id := len(Ipoints)
				Ipoints = append(Ipoints, &IpDat{e.Cell.Id, C, e.States[idx], nil})
				IpsBins.Append(C, id)
			}
		}
	}

	// success
	return true
}

func Splot(key string, loc LineLocator, times []float64, styles []*plt.LineData) {
}

func Plot(keyx, keyy string, loc PointLocator, styles []*plt.LineData) {
}

// Show shows plot
//  extra -- is a function to carry out extra configurations
func Show(extra func(spd SubpDat)) (err error) {
	spd, err := plot_all()
	if err != nil {
		return
	}
	if extra != nil {
		extra(spd)
	}
	plt.Show()
	return
}

// Save saves plot
func Save(dirout, filename string, extra func(spd SubpDat)) (err error) {
	spd, err := plot_all()
	if err != nil {
		return
	}
	if extra != nil {
		extra(spd)
	}
	plt.SaveD(dirout, filename)
	return
}

// read_results reads all results
//  Note: R[nkeys][nqts...?][ntimes]
func read_results() (T []float64, R [][][]float64, err error) {
	R = TplotStart()
	T = make([]float64, Sum.NumTidx)
	for tidx := 0; tidx < Sum.NumTidx; tidx++ {
		if !Dom.ReadSol(tidx) {
			return nil, nil, utl.Err("ReadSol failed. See log files\n")
		}
		T[tidx] = Dom.Sol.T
		for i, dat := range TplotData {
			for j, q := range dat.Qts {
				R[i][j][tidx] = *q.Value
			}
		}
	}
	return
}

// plot_all plots all results. It returns a map of keys and subplot indices. Ex:
//  spd: "pl" => 1,1,2 == row,col,idx
func plot_all() (spd SubpDat, err error) {
	T, R, err := read_results()
	if err != nil {
		return
	}
	nplots := len(R)
	nrow, ncol := utl.BestSquare(nplots)
	if SubpNrow > 0 {
		nrow = SubpNrow
	}
	if SubpNcol > 0 {
		ncol = SubpNcol
	}
	spd = make(map[string][]int)
	for i := 0; i < nplots; i++ {
		key := TplotKeys[i]
		spd[key] = []int{nrow, ncol, i + 1}
		plt.Subplot(nrow, ncol, i+1)
		for j, Y := range R[i] {
			sty := TplotData[i].Sty[j]
			args := sty.GetArgs("clip_on=0")
			plt.Plot(T, Y, args)
		}
	}
	return
}
