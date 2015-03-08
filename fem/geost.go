// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"sort"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/ode"
)

// geostate holds state @ top of layer
type geostate struct{ pl, ρL, sl, ρ, σV float64 }

// GeoLayer holds information of one soil layer. It computes pressures (σVabs, pl) and
// densities (ρL, ρ) based on the following model
//
//    ρL  = ρL0 + Cl・pl   thus   dρL/dpl = Cl
//    ρ   = nf・sl・ρL + (1 - nf)・ρS
//    ns  = 1 - nf
//
//    Z(z) = zmax + T・(z - zmax)   with 0 ≤ T ≤ 1
//    dZ   = (z - zmax)・dT
//    dpl  = ρL(pl)・g・(-dZ)
//    dpl  = ρL(pl)・g・(zmax - z)・dT
//    dsl  = -Cc・dpl
//    dσV  = ρ(pl,sl)・g・(zmax - z)・dT
//    Δz   = zmax - z
//
//            / dpl/dT \   / ρL(pl)・g・Δz                   \
//            | dρL/dT |   | Cl・dpl/dT                      |
//    dY/dT = | dsl/dT | = | -Cc・dpl/dT                     |
//            | dρ/dT  |   | nf・dsl/dT・ρL + nf・sl・dρL/dT |
//            \ dσV/dT /   \ ρ(pl,sl)・g・Δz                 /
//
type GeoLayer struct {
	Tags  []int      // tags of cells within this layer
	Zmin  float64    // coordinate (elevation) at bottom of layer
	Zmax  float64    // coordinate (elevation) at top of layer
	Nodes []*Node    // nodes in layer
	Elems []Elem     // elements in layer
	Cl    float64    // liquid compressibility
	RhoS0 float64    // initial density of solids
	nf0   float64    // initial (constant) porosity
	K0    float64    // earth-pressure at rest
	Dpl   float64    // liquid pressure added by this layer
	DsigV float64    // absolute value of vertical stress increment added by this layer
	top   *geostate  // state @ top of layer
	fcn   ode.Cb_fcn // function for ode solver
	Jac   ode.Cb_jac // Jacobian for ode solver
	sol   ode.ODE    // ode solver
}

// Start starts ODE solver for computing state variables in Calc
//  prev -- previous state @ top of this layer
func (o *GeoLayer) Start(prev *geostate) {

	// set state @ top
	o.top = prev

	// x := {pl, ρL, sl, ρ, σV} == geostate
	g := Global.Sim.Gfcn.F(0, nil)
	nf := o.nf0
	o.fcn = func(f []float64, x float64, y []float64, args ...interface{}) error {
		Δz := args[0].(float64)
		ρL := y[1]
		sl := y[2]
		ρ := y[3]
		Cc := 0.0
		f[0] = ρL * g * Δz             // dpl/dT
		f[1] = o.Cl * f[0]             // dρL/dT
		f[2] = -Cc * f[0]              // dsl/dT
		f[3] = nf*f[2]*ρL + nf*sl*f[1] // dρ/dT
		f[4] = ρ * g * Δz              // dσV/dT
		return nil
	}

	// set ODE (using numerical Jacobian)
	silent := true
	o.sol.Init("Radau5", 5, o.fcn, nil, nil, nil, silent)
	o.sol.Distr = false // must be sure to disable this; otherwise it causes problems in parallel runs
}

// Calc computes state @ level z
func (o GeoLayer) Calc(z float64) (*geostate, error) {
	y := []float64{o.top.pl, o.top.ρL, o.top.sl, o.top.ρ, o.top.σV}
	Δz := o.Zmax - z
	err := o.sol.Solve(y, 0, 1, 1, false, Δz)
	if err != nil {
		err = chk.Err("geost: failed when calculating state sing ODE solver: %v", err)
		return nil, err
	}
	return &geostate{y[0], y[1], y[2], y[3], y[4]}, nil
}

// GeoLayers is a set of Layer
type GeoLayers []*GeoLayer

// Len the length of Layers
func (o GeoLayers) Len() int {
	return len(o)
}

// Swap swaps two Layers
func (o GeoLayers) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Layers: sort from top to bottom
func (o GeoLayers) Less(i, j int) bool {
	return o[i].Zmin > o[j].Zmin
}

// SetGeoSt sets the initial state to a hydrostatic condition
func (o *Domain) SetGeoSt(stg *inp.Stage) (ok bool) {

	// check layers definition
	geo := stg.GeoSt
	if LogErrCond(len(geo.Layers) < 1, "geost: layers must be defined by stating what tags belong to which layer") {
		return
	}

	// get region
	if LogErrCond(len(Global.Sim.Regions) != 1, "geost: can only handle one domain for now") {
		return
	}
	reg := Global.Sim.Regions[0]

	// fix UseK0
	nlayers := len(geo.Layers)
	if len(geo.UseK0) != nlayers {
		geo.UseK0 = make([]bool, nlayers)
	}

	// initialise layers
	var L GeoLayers
	L = make([]*GeoLayer, nlayers)
	ndim := Global.Ndim
	nodehandled := make(map[int]bool)
	for i, tags := range geo.Layers {

		// new layer
		L[i] = new(GeoLayer)
		L[i].Tags = tags
		L[i].Zmin = Global.Sim.MaxElev
		L[i].Zmax = 0
		L[i].Cl = Global.Sim.WaterRho0 / Global.Sim.WaterBulk

		// get porous parameters
		if !L[i].get_porous_parameters(reg, tags[0]) {
			return
		}

		// parameters
		if geo.UseK0[i] {
			L[i].K0 = geo.K0[i]
		} else {
			L[i].K0 = geo.Nu[i] / (1.0 - geo.Nu[i])
		}
		if LogErrCond(L[i].K0 < 1e-7, "geost: K0 or Nu is incorect: K0=%g, Nu=%g", L[i].K0, geo.Nu) {
			return
		}

		// for each tag of cells in this layer
		for _, tag := range tags {

			// check tags
			cells := o.Msh.CellTag2cells[tag]
			if LogErrCond(len(cells) < 1, "geost: there are no cells with tag = %d", tag) {
				return
			}

			// set nodes and elements and find min and max z-coordinates
			for _, c := range cells {
				L[i].Elems = append(L[i].Elems, o.Cid2elem[c.Id])
				for _, v := range c.Verts {
					if !nodehandled[v] {
						L[i].Nodes = append(L[i].Nodes, o.Vid2node[v])
					}
					L[i].Zmin = min(L[i].Zmin, o.Msh.Verts[v].C[ndim-1])
					L[i].Zmax = max(L[i].Zmax, o.Msh.Verts[v].C[ndim-1])
					nodehandled[v] = true
				}
			}
		}
	}

	// sort layers from top to bottom
	sort.Sort(L)
	//io.Pfyel("layers = %v\n", L)

	// set previous/top states in layers and compute Sol.Y
	var err error
	for i, lay := range L {

		// previous state
		var top *geostate
		ρS := lay.RhoS0
		nf := lay.nf0
		if i == 0 {
			pl := 0.0
			ρL := Global.Sim.WaterRho0
			sl := 1.0
			ρ := nf*sl*ρL + (1.0-nf)*ρS
			σV := 0.0
			top = &geostate{pl, ρL, sl, ρ, σV}
		} else {
			top, err = L[i-1].Calc(L[i-1].Zmin)
			if LogErr(err, "cannot compute state @ bottom of layer") {
				return
			}
			ρL := top.ρL
			sl := top.sl
			top.ρ = nf*sl*ρL + (1.0-nf)*ρS
		}

		// start layer
		//io.PfYel("top = %v\n", top)
		lay.Start(top)

		// set nodes
		for _, nod := range lay.Nodes {
			z := nod.Vert.C[ndim-1]
			s, err := lay.Calc(z)
			if LogErr(err, io.Sf("cannot compute state @ node z = %g", z)) {
				return
			}
			dof := nod.GetDof("pl")
			if dof != nil {
				o.Sol.Y[dof.Eq] = s.pl
			}
		}

		// set elements
		for _, elem := range lay.Elems {
			if ele, okk := elem.(ElemIntvars); okk {

				// build slices
				coords := ele.Ipoints()
				nip := len(coords)
				pl := make([]float64, nip)
				ρL := make([]float64, nip)
				sx := make([]float64, nip)
				sy := make([]float64, nip)
				sz := make([]float64, nip)
				for i := 0; i < nip; i++ {
					z := coords[i][ndim-1]
					s, err := lay.Calc(z)
					if LogErr(err, io.Sf("cannot compute state @ ip z = %g", z)) {
						return
					}
					pl[i], ρL[i] = s.pl, s.ρL
					p := s.pl * s.sl
					σVe := -s.σV + p
					σHe := lay.K0 * σVe
					sx[i], sy[i], sz[i] = σHe, σVe, σHe
					if ndim == 3 {
						sx[i], sy[i], sz[i] = σHe, σHe, σVe
					}
				}
				ivs := map[string][]float64{"pl": pl, "ρL": ρL, "sx": sx, "sy": sy, "sz": sz}

				// set element's states
				if LogErrCond(!ele.SetIniIvs(o.Sol, ivs), "geost: element's internal values setting failed") {
					return
				}
			}
		}
	}
	return true
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////////

// get_porous_parameters extracts parameters based on region data
func (o *GeoLayer) get_porous_parameters(reg *inp.Region, ctag int) (ok bool) {
	edat := reg.Etag2data(ctag)
	mat := Global.Sim.Mdb.Get(edat.Mat)
	if mat.Model != "group" {
		LogErrCond(true, "geost: material type describing layer must be 'group' with porous data")
		return
	}
	//var RhoL0, BulkL float64
	if matname, found := io.Keycode(mat.Extra, "p"); found {
		m := Global.Sim.Mdb.Get(matname)
		for _, p := range m.Prms {
			switch p.N {
			case "RhoS0":
				o.RhoS0 = p.V
			case "nf0":
				o.nf0 = p.V
			}
		}
	}
	if LogErrCond(o.RhoS0 < 1e-7, "geost: initial density of solids RhoS0=%g is incorrect", o.RhoS0) {
		return
	}
	if LogErrCond(o.nf0 < 1e-7, "geost: initial porosity nf0=%g is incorrect", o.nf0) {
		return
	}
	return true
}

// String prints geostate
func (o *geostate) String() string {
	return io.Sf("pl=%g ρL=%g sl=%g ρ=%g σV=%g\n", o.pl, o.ρL, o.sl, o.ρ, o.σV)
}

// String prints a json formatted string with GeoLayers' content
func (o GeoLayers) String() string {
	if len(o) == 0 {
		return "[]"
	}
	l := "[\n"
	for i, lay := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("  { \"Tags\":%v, \"Zmin\":%g, \"Zmax\":%g, \"nf0\":%g, \"RhoS0\":%g, \"Cl\":%g\n", lay.Tags, lay.Zmin, lay.Zmax, lay.nf0, lay.RhoS0, lay.Cl)
		l += "    \"Nodes\":["
		for j, nod := range lay.Nodes {
			if j > 0 {
				l += ","
			}
			l += io.Sf("%d", nod.Vert.Id)
		}
		l += "],\n    \"Elems\":["
		for j, ele := range lay.Elems {
			if j > 0 {
				l += ","
			}
			l += io.Sf("%d", ele.Id())
		}
		l += "] }"
	}
	l += "\n]"
	return l
}
