// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"strings"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/gm"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// global data
var (

	// FEM data
	Ips    []*IpDat // data at ips
	IpBins gm.Bins  // bins for nodes

	// containers
	TrkPts     []*PointTracker  // tracking individual points
	TrkListPts []*PointsTracker // tracking collection of points

	// subplots
	Subplots    []*SubPlotData // all subplots
	ConfigFunc  func()         // configuration function
	SubplotRows int            // number of rows in figure
	SubplotCols int            // number of cols in figure
	Tserie      []float64
	CSubplot    *SubPlotData // current subplot
)

// IpDat saves data from ips
type IpDat struct {
	Labels []string   // labels
	X      []float64  // coordinates
	Data   []*float64 // element's state connections
}

// PointTracker saves information for one tracked point along time
type PointTracker struct {

	// essential
	Tag    string // point tag
	NodeId int    // point id
	IpId   int    // integration point id

	// derived
	Keys []string
	Data map[string][]float64 // [nkey][ntimes]
	X    []float64            // coordinates
}

// AlongLine defines a selection along line
type AlongLine struct {
	A []float64 // first point on line
	B []float64 // second point on line
}

// PointsTracker saves information for a group of tracked points along time
type PointsTracker struct {
	Tag      string                 // point id
	IdOrTags []int                  // id or tags
	Along    AlongLine              // along line definition
	NodePts  []int                  // container for nodes
	IpPts    []int                  // container for integration points
	Keys     []string               // keys
	Data     []map[string][]float64 // [npoints][nkey][ntimes]
	X        [][]float64            // [npoints][dir] slice of coordinates
}

// Series stores all data for a plot series (X vs Y)
type Series struct {
	Xkey   string
	Ykey   string
	Tag    string
	Label  string
	Format string
	X      []float64
	Y      []float64
	Tidx   int
}

// SubPlotData stores all data for one subplot
type SubPlotData struct {
	Title     string
	Xlbl      string
	Ylbl      string
	Xunt      string
	Yunt      string
	Xscale    float64
	Yscale    float64
	AllSeries []*Series
}

// SubPlot resets settings to start a new subplot
func Subplot() {
	// new subplot
	CSubplot = new(SubPlotData)

	// add the new subplot
	Subplots = append(Subplots, CSubplot)
}

func SetSubplots(nrows, ncols int) {
	SubplotRows = nrows
	SubplotCols = ncols
}

func SetTitle(title string) {
	if CSubplot == nil {
		Subplot()
	}
	CSubplot.Title = title
}

func SetAxisLabels(xlbl, ylbl string) {
	if CSubplot == nil {
		Subplot()
	}
	CSubplot.Xlbl = xlbl
	CSubplot.Ylbl = ylbl
}

func SetDataScale(xs, ys float64) {
	if CSubplot == nil {
		Subplot()
	}
	CSubplot.Xscale = xs
	CSubplot.Yscale = ys
}

// SetPoint select a point for tracking
func SetPoint(id int, ptag string) {
	// add node
	x := Dom.Nodes[id].Vert.C
	pt := &PointTracker{Tag: ptag, NodeId: id, IpId: -1, X: x}
	pt.Data = make(map[string][]float64, 0)
	TrkPts = append(TrkPts, pt)
}

// SetPointAt select a point for tracking
func SetPointAt(point []float64, ptag string) {
	//pt := &PointTracker{Tag: ptag, NodeId: -1, IpId: -1, X: point}

	// add node
	id := NodBins.Find(point)
	if id > 0 {
		pt := &PointTracker{Tag: ptag, NodeId: id, IpId: -1, X: point}
		pt.Data = make(map[string][]float64, 0)
		TrkPts = append(TrkPts, pt)
	}

	// add ip
	id = IpBins.Find(point)
	if id > 0 {
		pt := &PointTracker{Tag: ptag, NodeId: -1, IpId: id, X: point}
		pt.Data = make(map[string][]float64, 0)
		TrkPts = append(TrkPts, pt)
	}
}

// SetPointAt select a group of points for tracking
func SetPoints(ids_tags []int, ptag string) {
	//ptst := &PointsTracker{Tag: ptag, IdOrTags: ids_tags}

	ptst := &PointsTracker{}
	ptst.Tag = ptag

	for _, idortag := range ids_tags {
		if idortag < 0 {
			//add tags
			tag := idortag
			verts := Dom.Msh.VertTag2verts[tag]
			for _, v := range verts {
				id := v.Id
				ptst.X = append(ptst.X, v.C)
				ptst.NodePts = append(ptst.NodePts, id)
			}
		} else {
			//add nodes
			id := idortag
			x := Dom.Nodes[id].Vert.C
			ptst.X = append(ptst.X, x)
			ptst.NodePts = append(ptst.NodePts, id)
		}
	}
	TrkListPts = append(TrkListPts, ptst)
}

// SetPointAt select a group of points for tracking
func SetAlong(A, B []float64, ptag string) {
	//ptst := &PointsTracker{Tag: ptag, Along: Along{A, B}}

	ptst := &PointsTracker{}
	ptst.Tag = ptag

	// add nodes
	ids := NodBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		x := Dom.Nodes[id].Vert.C
		ptst.X = append(ptst.X, x)
		ptst.NodePts = append(ptst.NodePts, id)
	}

	// add ips
	ids = IpBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		// x=??
		ptst.NodePts = append(ptst.NodePts, id)
	}

	TrkListPts = append(TrkListPts, ptst)
}

func GetCoordByKey(dir, key string, tag string) []float64 {
	cidx := 0
	if dir == "y" {
		cidx = 1
	}
	if dir == "z" {
		cidx = 2
	}

	crds := []float64{}
	for _, tps := range TrkListPts {
		if tps.Tag == tag {
			for i, dmap := range tps.Data {
				_, ok := dmap[key]
				if ok {
					crds = append(crds, tps.X[i][cidx])
				}
			}
			break
		}
	}
	return crds
}

func GetPair(dir, key string, tag string) (a, b []float64) {
	cidx := 0
	if dir == "y" {
		cidx = 1
	}
	if dir == "z" {
		cidx = 2
	}

	tidx := 0
	crds := []float64{}
	vals := []float64{}

	for _, tps := range TrkListPts {
		if tps.Tag == tag {
			for i, dmap := range tps.Data {
				vec, ok := dmap[key]
				if ok {
					crds = append(crds, tps.X[i][cidx])
					vals = append(vals, vec[tidx])
				}
			}
			break
		}
	}
	return crds, vals
}

func Get(key, tag string) []float64 {

	for _, tp := range TrkPts {
		if tp.Tag == tag {
			return tp.Data[key]
		}
	}

	tidx := 0
	vals := []float64{}

	for _, tps := range TrkListPts {
		if tps.Tag == tag {
			for _, dmap := range tps.Data {
				vec, ok := dmap[key]
				if ok {
					vals = append(vals, vec[tidx])
				}
			}
			return vals
		}
	}

	chk.Panic("Error: no data found for key %v with tag %v", key, tag)
	return []float64{}
}

type PlotOpts struct {
	LineWidth  float64
	MarkerSize float64
	Times      []float64
}

func Plot(xhdl, yhdl interface{}, ptag, format, label string, opt_times ...[]float64) {
	// create subplot if none
	if CSubplot == nil {
		Subplot()
	}

	// check ptag
	if ptag == "" {
		chk.Panic("Error: no tag provided for plot")
	}

	// set default label
	if label == "" {
		label = ptag
	}

	// check if plot is for various times
	single_time := true
	var times []float64
	var labels []string
	if len(opt_times) > 0 {
		single_time = false
		times = opt_times[0]
		labels = strings.Split(label, ",")
	}
	_ = times
	_ = labels

	if single_time {

		// fill a new Series
		srs := new(Series)
		// get handlers for X and Y data: string or []float64
		srs.Xkey, _ = xhdl.(string)
		srs.Ykey, _ = yhdl.(string)
		srs.X, _ = xhdl.([]float64)
		srs.Y, _ = yhdl.([]float64)
		srs.Tag = ptag
		srs.Format = format
		srs.Label = label
		CSubplot.AllSeries = append(CSubplot.AllSeries, srs)

		// check
		if (srs.Xkey == "" && len(srs.X) == 0) || (srs.Ykey == "" && len(srs.Y) == 0) {
			chk.Err("wrong handlers for %v in Plot", ptag)
		}

		// check subplot labels
		if CSubplot.Xlbl == "" {
			CSubplot.Xlbl = srs.Xkey
		}
		if CSubplot.Ylbl == "" {
			CSubplot.Ylbl = srs.Ykey
		}
	}

	//if strings.Contains("x y z", srs.Xkey)
}

func Show() {
	// complete all data in all subplots and series
	for _, subplot := range Subplots {
		for _, srs := range subplot.AllSeries {
			if srs.Xkey == "t" {
				srs.X = Tserie
			}

			if srs.Ykey == "t" {
				srs.Y = Tserie
			}

			if strings.Contains("x y z", srs.Xkey) && len(srs.X) == 0 {
				srs.X = GetCoordByKey(srs.Xkey, srs.Ykey, srs.Tag)
				srs.Y = Get(srs.Ykey, srs.Tag)
			}

			if strings.Contains("x y z", srs.Ykey) && len(srs.Y) == 0 {
				srs.X = Get(srs.Xkey, srs.Tag)
				srs.Y = GetCoordByKey(srs.Ykey, srs.Xkey, srs.Tag)
			}

			if len(srs.X) == 0 {
				srs.X = Get(srs.Xkey, srs.Tag)
			}

			if len(srs.Y) == 0 {
				srs.Y = Get(srs.Ykey, srs.Tag)
			}
		}
	}

	// Show
	nplots := len(Subplots)
	nrow, ncol := utl.BestSquare(nplots)

	nrow = imax(nrow, SubplotRows)
	ncol = imax(ncol, SubplotCols)
	for i, subplot := range Subplots {
		plt.Subplot(nrow, ncol, i+1)
		for _, srs := range subplot.AllSeries {
			plt.Plot(srs.X, srs.Y, srs.Format)
		}
		plt.Gll(subplot.Xlbl, subplot.Ylbl, "")
		plt.Title(subplot.Title, "")
	}
	plt.Show()
}

// Start starts handling of results given a simulation input file
func Read(simfnpath string, stageIdx, regionIdx int) {

	// constants
	Ndiv := 20

	// start FE global structure
	erasefiles := false
	verbose := false
	if !fem.Start(simfnpath, erasefiles, verbose) {
		chk.Panic("Error in fem.Start")
	}

	// read summary
	Sum = fem.ReadSum()
	if Sum == nil {
		chk.Panic("\nError reading fem summary")
	}

	// allocate domain
	Dom = fem.NewDomain(fem.Global.Sim.Regions[regionIdx])
	if !Dom.SetStage(stageIdx, fem.Global.Sim.Stages[stageIdx]) {
		chk.Panic("\nError setting stage %v", stageIdx)
		return
	}

	// clear previous data
	Ips = make([]*IpDat, 0)

	// bins
	m := Dom.Msh
	xi := []float64{m.Xmin, m.Ymin}
	xf := []float64{m.Xmax, m.Ymax}
	if m.Ndim == 3 {
		xi = append(xi, m.Zmin)
		xf = append(xf, m.Zmax)
	}
	NodBins.Init(xi, xf, Ndiv)
	IpBins.Init(xi, xf, Ndiv)

	// add nodes to bins
	for activeId, nod := range Dom.Nodes {
		err := NodBins.Append(nod.Vert.C, activeId)
		if err != nil {
			return
		}
	}

	// add integration points to slice of ips and to bins
	for _, ele := range Dom.Elems {
		data := ele.OutIpsData()
		for _, ipd := range data {
			id := len(Ips)
			var labels []string
			var values []*float64
			for key, val := range ipd.V {
				labels = append(labels, key)
				values = append(values, val)
			}
			Ips = append(Ips, &IpDat{labels, ipd.X, values})
			IpBins.Append(ipd.X, id)
		}
	}
}

func Apply() {
	ntidx := len(Sum.Times)
	for tidx := 0; tidx < ntidx; tidx++ {
		Tserie = append(Tserie, float64(tidx))
		if !Dom.In(tidx) {
			chk.Panic("Domain.In failed. See log files\n")
		}

		// collect data to Point Trackers
		for _, tp := range TrkPts {

			// data from nodes
			if tp.NodeId > 0 {
				id := tp.NodeId
				node := Dom.Nodes[id]
				keys := node.GetKeys()
				for _, key := range keys {
					dof := node.GetDof(key)
					utl.StrDblsMapAppend(&tp.Data, key, Dom.Sol.Y[dof.Eq])
				}
			}

			// data from ips
			if tp.IpId > 0 {
				id := tp.NodeId
				ip := Ips[id]
				keys := ip.Labels
				for i, key := range keys {
					utl.StrDblsMapAppend(&tp.Data, key, *ip.Data[i])
				}
			}
		}

		// collect data to 'Lists of Points' Trackers
		for _, tps := range TrkListPts {

			// data from nodes
			if len(tps.NodePts) > 0 {
				for j, id := range tps.NodePts {
					node := Dom.Nodes[id]
					keys := node.GetKeys()
					for _, key := range keys {
						dof := node.GetDof(key)
						tps.Data = append(tps.Data, make(map[string][]float64)) // important
						utl.StrDblsMapAppend(&tps.Data[j], key, Dom.Sol.Y[dof.Eq])
					}

				}

			}

			// data from ips
			if len(tps.IpPts) > 0 {
				for j, id := range tps.IpPts {
					ip := Ips[id]
					keys := ip.Labels
					for i, key := range keys {
						tps.Data = append(tps.Data, make(map[string][]float64)) // important
						utl.StrDblsMapAppend(&tps.Data[j], key, *ip.Data[i])
					}
				}
			}
		}
	}
}

func imax(a, b int) int {
	if a < b {
		return b
	}
	return a
}
