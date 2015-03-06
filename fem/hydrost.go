// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"
)

// SetHydroSt sets the initial state to a hydrostatic condition
func (o *Domain) SetHydroSt(stg *inp.Stage) (ok bool) {

	// check for hydrost data
	hst := stg.HydroSt
	if hst == nil {
		return true
	}

	// max elevation and water level
	ndim := o.Msh.Ndim
	zmax := o.Msh.Ymax
	if ndim == 3 {
		zmax = o.Msh.Zmax
	}
	zwater := zmax
	if hst.Zwater > zmax {
		zwater = hst.Zwater // e.g. ponding
	}

	// set Sol
	for _, n := range o.Nodes {
		z := n.Vert.C[ndim-1]
		dof := n.GetDof("pl")
		if dof != nil {
			pl := (zwater - z) * hst.GamW
			o.Sol.Y[dof.Eq] = pl
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
			pl[i] = (zwater - z) * hst.GamW
		}
		ivs := map[string][]float64{"pl": pl}

		// set element's states
		if LogErrCond(!e.SetIvs(ivs), "hydrostatic: element's internal values setting failed") {
			return
		}
	}

	// success
	log.Printf("dom: initial hydrostatic state set with zmax=%g zwater=%g γw=%g", zmax, zwater, hst.GamW)
	return true
}
