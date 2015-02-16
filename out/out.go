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

// Global variables
var (

	// constants
	TolC float64 // tolerance to compare x-y-z coordinates
	Ndiv int     // bins n-division

	// data
	Sum     *fem.Summary // summary of results
	Dom     *fem.Domain  // FE domain
	Ipoints []*IpDat     // all integration points
	NodBins gm.Bins      // bins for nodes
	IpsBins gm.Bins      // bins for integration points
)

// With starts handling and plotting of results given a simulation input file
// It returns a callback function that must be called in order to release resources and flush files
func With(simfnpath string, stageIdx, regionIdx int, commands func()) (ok bool) {

	// constants
	TolC = 1e-8
	Ndiv = 20

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		return
	}
	defer fem.End()

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

	// do plot
	commands()
	return true
}

func Splot(key string, loc LineLocator, times []float64, styles []*plt.LineData) {
}

func Plot(keyx, keyy string, loc PointLocator, styles []*plt.LineData) {
}

// Plotall plots all results. It returns a map of keys and subplot indices. Ex:
//  splots: "pl" => 1,1,2 == row,col,idx
func PlotAll() (splots map[string][]int, err error) {
	T, R, err := read_results()
	if err != nil {
		return
	}
	nplots := len(R)
	nrow, ncol := utl.BestSquare(nplots)
	splots = make(map[string][]int)
	for i := 0; i < nplots; i++ {
		key := TplotKeys[i]
		splots[key] = []int{nrow, ncol, i + 1}
		plt.Subplot(nrow, ncol, i+1)
		for j, Y := range R[i] {
			sty := TplotData[i].Sty[j]
			args := sty.GetArgs("clip_on=0")
			plt.Plot(T, Y, args)
		}
	}
	return
}

// Show shows plot
func Show() {
	plt.Show()
}

// Save saves plot
func Save(dirout, filename string) {
	plt.SaveD(dirout, filename)
}

// read_results reads all results
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
				R[i][j] = append(R[i][j], *q.Value)
			}
		}
	}
	return
}
