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

	// for time-plots

)

// With starts handling and plotting of results given a simulation input file
// It returs a callback function that must be called in order to release resources and flush files
func With(simfnpath string, stageIdx, regionIdx int) (err error) {

	// constants
	TolC = 1e-8
	Ndiv = 20

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		return utl.Err("cannot load sim file %q\n", simfnpath)
	}

	// read summary
	Sum = fem.ReadSum()
	if Sum == nil {
		return utl.Err("cannot load summary file %q\n", simfnpath)
	}

	// allocate domain
	Dom = fem.NewDomain(fem.Global.Sim.Regions[regionIdx])
	if !Dom.SetStage(stageIdx, fem.Global.Sim.Stages[stageIdx]) {
		return utl.Err("SetStage failed\n")
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
		err = NodBins.Append(n.Vert.C, activeId)
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
	return
}

func Splot(key string, loc LineLocator, times []float64, styles []*plt.LineData) {
}

func Plot(keyx, keyy string, loc PointLocator, styles []*plt.LineData) {
}

func Show() (err error) {

	// read results
	R := TplotStart()
	T := make([]float64, Sum.NumTidx)
	for tidx := 0; tidx < Sum.NumTidx; tidx++ {
		if !Dom.ReadSol(tidx) {
			return utl.Err("ReadSol failed. See log files\n")
		}
		T[tidx] = Dom.Sol.T
		for i, dat := range TplotData {
			for j, q := range dat.Qts {
				R[i][j] = append(R[i][j], *q.Value)
			}
		}
	}

	// plot
	nplots := len(R)
	nrow, ncol := utl.BestSquare(nplots)
	utl.Pforan("nplots = %v\n", nplots)
	for k := 0; k < nplots; k++ {
		key := TplotKeys[k]
		plt.Subplot(nrow, ncol, k+1)
		utl.Pforan("key = %v\n", key)
		for _, Y := range R[k] {
			plt.Plot(T, Y, "clip_on=0")
		}
	}
	utl.Pforan("nrow,ncol = %v, %v\n", nrow, ncol)
	plt.Show()
	return
}

func Save(eps bool) {
}
