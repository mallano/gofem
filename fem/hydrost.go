// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import "github.com/cpmech/gofem/inp"

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
	if hst.Unsat {
		zwater = hst.Zwater // e.g. internal water level specified for unsaturated condition
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
	return true
}
