// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/io"
)

// SetHydroSt sets the initial state to a hydrostatic condition
func (o *Domain) SetGeoSt(stg *inp.Stage) (ok bool) {

	// check for hydrost data
	geo := stg.GeoSt
	if geo == nil {
		return true
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
	io.PfYel("zwater = %v\n", zwater)
	io.PfYel("zmax   = %v\n", zmax)

	// set Sol
	for _, n := range o.Nodes {
		z := n.Vert.C[ndim-1]
		dof := n.GetDof("pl")
		if dof != nil {
			pl := (zwater - z) * geo.GamW
			o.Sol.Y[dof.Eq] = pl
			io.Pforan("z=%g pl=%v\n", z, pl)
		}
	}

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
