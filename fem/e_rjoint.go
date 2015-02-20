// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

// Rjoint implements the rod-joint (interface/link) element for reinforced solids
type Rjoint struct {

	// basic data
	Cell *inp.Cell   // cell
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Ndim int         // space dimension
	Nu   int         // total number of unknowns

	// parameters
	K1      float64
	K2      float64
	H       float64
	Coulomb bool

	// supporting elements
	Rod *Rod
	Sld *ElemU

	// optional data
	UseB      bool    // use B matrix
	Thickness float64 // thickness (for plane-stress)
	Debug     bool    // debugging flag

	// integration points
	IpsElem []*shp.Ipoint // integration points of element
	IpsFace []*shp.Ipoint // integration points corresponding to faces

	// material model and internal variables
	Model    msolid.Solid // material model
	MdlSmall msolid.Small // model specialisation for small strains
	MdlLarge msolid.Large // model specialisation for large deformations

	// internal variables
	States    []*msolid.State
	StatesBkp []*msolid.State

	// natural boundary conditions
	NatBcs []*NaturalBc

	// local starred variables
	ζs    [][]float64 // [nip][ndim] t2 star vars: ζ* = α1.u + α2.v + α3.a
	χs    [][]float64 // [nip][ndim] t2 star vars: χ* = α4.u + α5.v + α6.a
	divχs []float64   // [nip] divergent of χs (for coupled sims)

	// scratchpad. computed @ each ip
	grav []float64   // [ndim] gravity vector
	us   []float64   // [ndim] displacements @ ip
	fi   []float64   // [nu] internal forces
	K    [][]float64 // [nu][nu] consistent tangent (stiffness) matrix
	B    [][]float64 // [nsig][nu] B matrix for axisymetric case
	D    [][]float64 // [nsig][nsig] constitutive consistent tangent matrix

	// strains
	ε  []float64 // total (updated) strains
	Δε []float64 // incremental strains leading to updated strains

	// for debugging
	fex []float64 // x-components of external surface forces
	fey []float64 // y-components of external syrface forces
	fez []float64 // z-components of external syrface forces
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	iallocators["rjoint"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info {
		return &Info{}
	}

	// element allocator
	eallocators["rjoint"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {

		// basic data
		var o Rjoint
		o.Cell = msh.Cells[cid]

		// material model name
		matname := edat.Mat
		matdata := Global.Mdb.Get(matname)
		if LogErrCond(matdata == nil, "Mdb.Get failed\n") {
			return nil
		}
		mdlname := matdata.Model

		// model
		// TODO
		if false {
			o.Model = msolid.GetModel(Global.Sim.Data.FnameKey, matname, mdlname, false)
			if LogErrCond(o.Model == nil, "cannot find model named %s\n", mdlname) {
				return nil
			}
			err := o.Model.Init(o.Ndim, Global.Sim.Data.Pstress, matdata.Prms)
			if LogErr(err, "Model.Init failed") {
				return nil
			}
		}

		// parameters
		for _, p := range matdata.Prms {
			switch p.N {
			case "k1":
				o.K1 = p.V
			case "k2":
				o.K2 = p.V
			case "h":
				o.H = p.V
			case "mu":
				if p.V > 0.0 {
					o.Coulomb = true
				}
			}
		}

		// return new element
		return &o
	}
}

// connect rod and solid ////////////////////////////////////////////////////////////////////////////

// Connect connects rod/solid elements in this Rjoint
func (o *Rjoint) Connect(cid2elem []Elem) (nnzK int, ok bool) {

	// get rod and solid elements
	rodId := o.Cell.JlinId
	sldId := o.Cell.JsldId
	o.Rod = cid2elem[rodId].(*Rod)
	o.Sld = cid2elem[sldId].(*ElemU)
	if LogErrCond(o.Rod == nil, "cannot find joint's rod cell with id == %d", rodId) {
		return
	}
	if LogErrCond(o.Sld == nil, "cannot find joint's solid cell with id == %d", sldId) {
		return
	}

	io.Pforan("my rod = %+v\n", o.Rod)
	io.Pforan("my sld = %+v\n", o.Sld)

	// integration points
	// TODO
	//nip = len(o.IpsElem)

	// scratchpad. computed @ each ip

	return 0, true
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *Rjoint) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	return true
}

// SetEleConds set element conditions
func (o *Rjoint) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	return true
}

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *Rjoint) SetNatBcs(key string, idxface int, f fun.Func, extra string) (ok bool) {
	o.NatBcs = append(o.NatBcs, &NaturalBc{key, idxface, f, extra})
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *Rjoint) InterpStarVars(sol *Solution) (ok bool) {
	return true
}

// adds -R to global residual vector fb
func (o *Rjoint) AddToRhs(fb []float64, sol *Solution) (ok bool) {
	return true
}

// adds element K to global Jacobian matrix Kb
func (o *Rjoint) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {
	return true
}

// Update perform (tangent) update
func (o *Rjoint) Update(sol *Solution) (ok bool) {
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *Rjoint) InitIvs(sol *Solution) (ok bool) {
	return true
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *Rjoint) SetIvs(zvars map[string][]float64) (ok bool) {
	return true
}

// BackupIvs create copy of internal variables
func (o *Rjoint) BackupIvs() (ok bool) {
	return true
}

// RestoreIvs restore internal variables from copies
func (o *Rjoint) RestoreIvs() (ok bool) {
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o Rjoint) Encode(enc Encoder) (ok bool) {
	return true
}

// Decode decodes internal variables
func (o Rjoint) Decode(dec Decoder) (ok bool) {
	return true
}

// OutIpsData returns data from all integration points for output
func (o Rjoint) OutIpsData() (labels []string, data []*OutIpData) {
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *Rjoint) ipvars(idx int, sol *Solution) (ok bool) {
	return true
}
