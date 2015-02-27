// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gosl/io"
)

func geost_compute_rho(z, pl, pg float64, mdl *mporous.Model) (ρ float64) {

	// compute updated saturation
	divus := 0.0
	s, err := mdl.NewState(pl, pg, divus)
	if err != nil {
		io.PfRed("geost_compute_rho failed: %v\n", err)
		return
	}

	// saturations
	sl := s.Sl
	sg := 1.0 - sl

	// n variables
	nf := mdl.Nf0
	ns := s.Ns0
	nl := nf * sl
	ng := nf * sg

	// ρ variables
	ρl := nl * s.RhoL
	ρg := ng * s.RhoG
	ρs := ns * mdl.RhoS0
	ρ = ρl + ρg + ρs
	return
}

func geost_compute_p(z, pl, pg float64, mdl *mporous.Model) (p float64) {

	// compute updated saturation
	divus := 0.0
	s, err := mdl.NewState(pl, pg, divus)
	if err != nil {
		io.PfRed("geost_compute_p failed: %v\n", err)
		return
	}

	// saturations and voids' pressure
	sl := s.Sl
	sg := 1.0 - sl
	p = sl*pl + sg*pg
	return
}

type geost_fcn func(z float64) float64

// Layer holds information of one soil layer
type Layer struct {
	Tag   int       // element tag
	Zmin  float64   // min z-coordinate
	Zmax  float64   // max z-coordinate
	DsigV geost_fcn // ΔσV=f(z) : increment of total σV with elevation along this layer
	P     geost_fcn // p=f(x) : computes voids' pressure: p := sl * pl + sg * pg
}

// Layers is a set of Layer
type Layers []*Layer

// Len the length of Layers
func (o Layers) Len() int {
	return len(o)
}

// Swap swaps two Layers
func (o Layers) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Layers considering Dist
func (o Layers) Less(i, j int) bool {
	return o[i].Zmin < o[j].Zmin
}

// SetGeoSt sets the initial state to a hydrostatic condition
func (o *Domain) SetGeoSt(stg *inp.Stage) (ok bool) {

	// check for geost data
	geo := stg.GeoSt
	if geo == nil {
		return true
	}

	// check for homogeneous flag
	if !geo.Hom {
		LogErrCond(true, "geost: can only handle homogeneous domains for now")
		return
	}

	// max elevation and water level
	ndim := o.Msh.Ndim
	zmax := o.Msh.Ymax
	if ndim == 3 {
		zmax = o.Msh.Zmax
	}
	zwater := zmax
	if geo.Zwater > zmax {
		zwater = geo.Zwater // e.g. ponding
	}
	if geo.Unsat {
		zwater = geo.Zwater // e.g. internal water level specified for unsaturated condition
	}

	// set Sol
	for _, n := range o.Nodes {
		z := n.Vert.C[ndim-1]
		dof := n.GetDof("pl")
		if dof != nil {
			pl := (zwater - z) * geo.GamW
			o.Sol.Y[dof.Eq] = pl
		}
	}

	// get region
	if LogErrCond(len(Global.Sim.Regions) != 1, "geost: can only handle one domain for now") {
		return
	}
	reg := Global.Sim.Regions[0]

	// number of element tags
	ntags := len(reg.ElemsData)

	// find layers
	var layers Layers
	tag2lay := make(map[int]*Layer)
	for _, d := range reg.ElemsData {

		// get gravity
		if LogErrCond(len(Global.Sim.Stages) < 1, "geost: the first stage must be defined") {
			return
		}
		stg := Global.Sim.Stages[0]
		if LogErrCond(len(stg.EleConds) != ntags, "geost: the number of element conditions in stage definition must be equal to the number of element tags.") {
			return
		}
		var grav float64
		var foundgrav bool
		for _, econd := range stg.EleConds {
			if econd.Tag == d.Tag {
				for j, key := range econd.Keys {
					if key == "g" {
						fcn := Global.Sim.Functions.Get(econd.Funcs[j])
						if LogErrCond(fcn == nil, "geost: cannot find gravity function in functions database") {
							return
						}
						grav = fcn.F(0, nil)
						foundgrav = true
						break
					}
				}
				break
			}
		}
		if LogErrCond(!foundgrav, "geost: cannot find gravity function/value corresponding to element tag = %d", d.Tag) {
			return
		}

		// find min and max z-coordinates
		cells := o.Msh.CellTag2cells[d.Tag]
		if LogErrCond(len(cells) < 1, "geost: there are no cells with tag = %d", d.Tag) {
			return
		}
		layerZmin := o.Msh.Verts[cells[0].Verts[0]].C[ndim-1]
		layerZmax := layerZmin
		for _, c := range cells {
			for _, v := range c.Verts {
				layerZmin = min(layerZmin, o.Msh.Verts[v].C[ndim-1])
				layerZmax = max(layerZmax, o.Msh.Verts[v].C[ndim-1])
			}
		}

		// get model
		cell := cells[0]
		elem := o.Cid2elem[cell.Id]
		var mdl *mporous.Model
		switch e := elem.(type) {
		case *ElemUP:
			mdl = e.P.Mdl
		default:
			LogErrCond(true, "geost: cannot handle element type %v. only \"up\" is available now", cell.Type)
			return
		}

		// DsigV function
		DsigV := func(z float64) float64 {
			if z < layerZmin || z > layerZmax {
				io.PfRed("geost: wrong elevation in layer detected: z=%g is outside [%g, %g]\n", z, layerZmin, layerZmax)
				return 0
			}
			pl := (zwater - z) * geo.GamW
			pg := 0.0
			ρ := geost_compute_rho(z, pl, pg, mdl)
			ΔσV := (layerZmax - z) * ρ * grav
			return ΔσV
		}

		// compute increment of total vertical stress
		σtop := DsigV(layerZmax)
		σbot := DsigV(layerZmin)
		ΔσV := σbot - σtop
		io.Pforan("zmin=%g zmax=%g σtop=%g σbot=%g ΔσV=%g\n", layerZmin, layerZmax, σtop, σbot, ΔσV)

		// new layer
		//lay := &Layer{d.Tag, zmin, zmax, ΔσV, rhofcn, pfcn}
		lay := &Layer{}
		n := len(layers)
		layers = append(layers, lay)
		tag2lay[d.Tag] = layers[n]
	}

	// debug
	io.Pforan("layers =\n%+v\n", layers)
	//for _, e := range o.Elems {
	//io.Pforan("e = %v\n", e)
	//}
	return

	// set elements
	for _, e := range o.ElemIntvars {

		// get element's integration points data
		ele := e.(Elem)
		d := ele.OutIpsData()
		nip := len(d)

		// build map with pressures @ ips
		pl := make([]float64, nip)
		for i, ip := range d {
			z := ip.X[ndim-1]
			pl[i] = (zwater - z) * geo.GamW
		}
		io.Pfcyan("pl @ ip = %v\n", pl)
		ivs := map[string][]float64{"pl": pl}

		// set element's states
		if LogErrCond(!e.SetIvs(ivs), "hydrostatic: element's internal values setting failed") {
			return
		}
	}

	// success
	log.Printf("dom: initial hydrostatic state set with zmax=%g zwater=%g γw=%g", zmax, zwater, geo.GamW)
	return true
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////

func (o Layers) String() string {
	if len(o) == 0 {
		return "[]"
	}
	l := "[\n"
	for i, lay := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("  { \"Tag\":%d, \"Zmin\":%g, \"Zmax\":%g, \"DsigV\":%g }", lay.Tag, lay.Zmin, lay.Zmax, lay.DsigV)
	}
	l += "\n]"
	return l
}
