// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/gm"
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
	Dom     *fem.Domain // FE domain
	Ipoints []*IpDat    // all integration points
	NodBins gm.Bins     // bins for nodes
	IpsBins gm.Bins     // bins for integration points
	SelNodT []*fem.Node // selected nodes to plot along time
	SelIpsT []*IpDat    // selected ips to plot along time
)

// With starts handling and plotting of results given a simulation input file
// It returs a callback function that must be called in order to release resources and flush files
func With(simfnpath, resdir string, stageIdx, regionIdx int) func() {

	// constants
	TolC = 1e-8
	Ndiv = 20

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		return errFunc(utl.Err("cannot load sim file %q\n", simfnpath))
	}

	// allocate domain
	Dom = fem.NewDomain(fem.Global.Sim.Regions[regionIdx])
	if !Dom.SetStage(stageIdx, fem.Global.Sim.Stages[stageIdx]) {
		return errFunc(utl.Err("SetStage failed\n"))
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
			return errFunc(err)
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

	// return callback function
	return fem.End
}

func Tplot(key string, loc Locator, t0, tf float64, styles, labels []string, extra string) {
}

func Splot(key string, times []float64, styles, labels, extraeach []string, extraall string) {
}

func Show() {

	//nplots := 4
	//nrow, ncol := utl.BestSquare(nplots)
}

func Save(eps bool) {
}
