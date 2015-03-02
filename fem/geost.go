// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"
	"sort"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

// Layer holds information of one soil layer
type Layer struct {
	Tags   []int          // element's tags (elements in layer)
	Ndim   int            // space dimension
	GamW   float64        // unit weight of water
	Zwater float64        // water-table
	Zmin   float64        // min z-coordinate of layer
	Zmax   float64        // max z-coordinate of layer
	Grav   float64        // gravity constant
	Mdl    *mporous.Model // material model
	K0     float64        // earth-pressure @ rest
	Nu     float64        // to compute out-of-plane stress in plane strain case
	DsigV  float64        // ΔσV : total σV added by this layer. negative (-) means compression
}

/* State computes the stress values at elevation z
 *
 *  Input:
 *   σ0abs  -- stress added by layers above this layer (absolute value) or surcharge load (positive)
 *   z      -- elevation
 *
 *  Output:
 *   σV       -- total vertical stress
 *   sx sy sz -- effective stresses
 *   pl pg    -- liquid and gas pressures
 *
 *  Note that:
 *   σh' = K0・σv'
 *   σh' = ν・σv' / (1 - ν)
 *   σh' - ν・σh' = ν・σv'
 *   σh' = ν・(σv' + σh')
 *  Thus:, with σr' representing the out-of-plane effective stress in plane-strain case:
 *   σh' = σr'
 */
func (o Layer) State(σ0abs, z float64) (σV, sx, sy, sz, pl, pg float64, err error) {

	// check
	if z < o.Zmin || z > o.Zmax {
		err = chk.Err("geost: wrong elevation in layer detected: z=%g is outside [%g, %g]\n", z, o.Zmin, o.Zmax)
		return
	}

	// pressures and mixture density
	pl = (o.Zwater - z) * o.GamW
	pg = 0.0
	ρ, p := geost_compute_rho_and_p(z, pl, pg, o.Mdl)

	// stresses
	ΔσVabs := (o.Zmax - z) * ρ * o.Grav
	σVabs := σ0abs + ΔσVabs
	σV = -σVabs
	σVe := σV + p
	σHe := o.K0 * σVe

	// y == vertical
	if o.Ndim == 2 {
		sx, sy, sz = σHe, σVe, σHe
		return
	}

	// z == vertical
	sx, sy, sz = σHe, σHe, σVe
	return
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

// Less compares Layers: sort from top to bottom
func (o Layers) Less(i, j int) bool {
	return o[i].Zmin > o[j].Zmin
}

// SetGeoSt sets the initial state to a hydrostatic condition
func (o *Domain) SetGeoSt(stg *inp.Stage) (ok bool) {

	// check for geost data
	geo := stg.GeoSt
	if geo == nil {
		return true
	}

	// max elevation and water level
	ndim := o.Msh.Ndim
	zmax, zwater := get_zmax_zwater(geo, o.Msh)

	// set Sol.Y
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

	// initialise layers
	var layers Layers
	layers = make([]*Layer, len(geo.Layers))
	tag2lay := make(map[int]*Layer)
	for idx, laytags := range geo.Layers {

		// new layer
		var lay Layer
		tag0 := laytags[0]
		lay.Tags = laytags
		lay.Ndim = ndim
		lay.GamW = geo.GamW
		lay.Zwater = zwater
		lay.Zmin = zmax
		lay.Zmax = 0

		// gravity
		lay.Grav, ok = get_gravity(stg, tag0)
		if !ok {
			return
		}

		// get model
		edat := reg.Etag2data(tag0)
		lay.Mdl = GetAndInitPorousModel(edat.Mat)

		// parameters
		if geo.UseK0 {
			lay.K0 = geo.K0
			lay.Nu = geo.K0 / (1.0 + geo.K0)
		} else {
			lay.K0 = geo.Nu / (1.0 - geo.Nu)
			lay.Nu = geo.Nu
		}
		if LogErrCond(lay.K0 < 1e-7 || lay.Nu < 0 || lay.Nu > 0.4999, "geost: K0 or Nu is incorect: K0=%g, Nu=%g", lay.K0, lay.Nu) {
			return
		}

		// for each tag of cells in this layer
		for _, tag := range laytags {

			// check tags
			cells := o.Msh.CellTag2cells[tag]
			if LogErrCond(len(cells) < 1, "geost: there are no cells with tag = %d", tag) {
				return
			}

			// find min and max z-coordinates
			for _, c := range cells {
				for _, v := range c.Verts {
					lay.Zmin = min(lay.Zmin, o.Msh.Verts[v].C[ndim-1])
					lay.Zmax = max(lay.Zmax, o.Msh.Verts[v].C[ndim-1])
				}
			}

			// set mat
			tag2lay[tag] = layers[idx]
		}

		// increment of vertical stress added by this layer
		lay.DsigV = num.TrapzRange(lay.Zmin, lay.Zmax, 10, func(z float64) float64 {
			pl := (lay.Zwater - z) * lay.GamW
			pg := 0.0
			ρ, _ := geost_compute_rho_and_p(z, pl, pg, lay.Mdl)
			return ρ * lay.Grav
		})

		// append layer
		layers[idx] = &lay
	}

	// sort layers from top to bottom
	sort.Sort(layers)
	//io.Pforan("layers = %v\n", layers)

	// loop over layers, from top to bottom
	var σ0abs float64
	var err error
	for i, lay := range layers {
		if i > 0 {
			σ0abs += layers[i-1].DsigV
		}
		for _, tag := range lay.Tags {
			cells := o.Msh.CellTag2cells[tag]
			for _, c := range cells {
				elem := o.Cid2elem[c.Id]
				switch ele := elem.(type) {
				case ElemIntvars:

					// get element's integration points data
					e := ele.(Elem)
					d := e.OutIpsData()
					nip := len(d)

					// build maps of pressures and effective stresses
					sx := make([]float64, nip)
					sy := make([]float64, nip)
					sz := make([]float64, nip)
					pl := make([]float64, nip)
					pg := make([]float64, nip)
					for i, ip := range d {
						z := ip.X[ndim-1]
						_, sx[i], sy[i], sz[i], pl[i], pg[i], err = lay.State(σ0abs, z)
						if LogErr(err, "geost: cannot compute State") {
							return
						}

					}
					ivs := map[string][]float64{"sx": sx, "sy": sy, "sz": sz, "pl": pl, "pg": pg}

					// set element's states
					if LogErrCond(!ele.SetIvs(ivs), "geost: element's internal values setting failed") {
						return
					}
				}
			}
		}
	}

	// success
	log.Printf("dom: initial hydrostatic state set with zmax=%g zwater=%g γw=%g", zmax, zwater, geo.GamW)
	return true
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////

func geost_compute_rho_and_p(z, pl, pg float64, mdl *mporous.Model) (ρ, p float64) {

	// compute updated saturation
	divus := 0.0
	s, err := mdl.NewState(pl, pg, divus)
	if err != nil {
		io.PfRed("geost_compute_rho failed: %v\n", err)
		return
	}

	// saturations and voids' pressure
	sl := s.Sl
	sg := 1.0 - sl
	p = sl*pl + sg*pg

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

// get_zmax_zwater returns the max elevation and water level
func get_zmax_zwater(geo *inp.GeoStData, msh *inp.Mesh) (zmax, zwater float64) {
	ndim := msh.Ndim
	zmax = msh.Ymax
	if ndim == 3 {
		zmax = msh.Zmax
	}
	zwater = zmax
	if geo.Zwater > zmax {
		zwater = geo.Zwater // e.g. ponding
	}
	if geo.Unsat {
		zwater = geo.Zwater // e.g. internal water level specified for unsaturated condition
	}
	return
}

// get_gravity returns the initial gravity value defined in element's data structure
func get_gravity(stg *inp.Stage, elemTag int) (grav float64, ok bool) {
	if LogErrCond(len(Global.Sim.Stages) < 1, "geost: the first stage must be defined") {
		return
	}
	var foundgrav bool
	for _, econd := range stg.EleConds {
		if econd.Tag == elemTag {
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
	if LogErrCond(!foundgrav, "geost: cannot find gravity function/value corresponding to element tag = %d", elemTag) {
		return
	}
	return grav, true
}

// String prints a json formatted string with Layers' content
func (o Layers) String() string {
	if len(o) == 0 {
		return "[]"
	}
	l := "[\n"
	for i, lay := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("  { \"Tags\":%v, \"Zmin\":%g, \"Zmax\":%g, \"DsigV\":%g }", lay.Tags, lay.Zmin, lay.Zmax, lay.DsigV)
	}
	l += "\n]"
	return l
}
