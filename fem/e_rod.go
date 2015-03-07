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

// Rod represents a structural rod element (for only axial loads)
type Rod struct {

	// basic data
	Cid int         // cell/element id
	X   [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Shp *shp.Shape  // shape structure
	Nu  int         // total number of unknowns == 2 * nsn

	// parameters
	A float64 // cross-sectional area

	// variables for dynamics
	Rho  float64  // density of solids
	Gfcn fun.Func // gravity function

	// integration points
	IpsElem []*shp.Ipoint // integration points of element

	// vectors and matrices
	K   [][]float64 // global K matrix
	M   [][]float64 // global M matrices
	Rus []float64   // residual: Rus = fi - fx

	// problem variables
	Umap []int // assembly map (location array/element equations)

	// material model and internal variables
	Model     msolid.OnedSolid
	States    []*msolid.OnedState
	StatesBkp []*msolid.OnedState

	// scratchpad. computed @ each ip
	grav []float64 // [ndim] gravity vector
	us   []float64 // [ndim] displacements @ ip
	fi   []float64 // [nu] internal forces
	ue   []float64 // local u vector
}

// register element
func init() {

	// information allocator
	infogetters["rod"] = func(cellType string, faceConds []*FaceCond) *Info {

		// new info
		var info Info

		// number of nodes in element
		nverts := shp.GetNverts(cellType)

		// solution variables
		ykeys := []string{"ux", "uy"}
		if Global.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz"}
		}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// maps
		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz"}

		// t1 and t2 variables
		info.T2vars = ykeys
		return &info
	}

	// element allocator
	eallocators["rod"] = func(cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o Rod
		o.Cid = cid
		o.X = x
		o.Shp = shp.Get(cellType)
		ndim := Global.Ndim
		o.Nu = ndim * o.Shp.Nverts

		var err error

		// material model name
		matname := edat.Mat
		matdata := Global.Mdb.Get(matname)
		if LogErrCond(matdata == nil, "materials database failed on getting %q material\n", matname) {
			return nil
		}
		mdlname := matdata.Model
		o.Model = msolid.GetOnedSolid(Global.Sim.Data.FnameKey, matname, mdlname, false)
		if LogErrCond(o.Model == nil, "cannot find model named %s\n", mdlname) {
			return nil
		}
		err = o.Model.Init(ndim, matdata.Prms)
		if LogErr(err, "Model.Init failed") {
			return nil
		}

		// parameters
		for _, p := range matdata.Prms {
			switch p.N {
			case "A":
				o.A = p.V
			case "rho":
				o.Rho = p.V
			}
		}

		// integration points
		var nip int
		if s_nip, found := io.Keycode(edat.Extra, "nip"); found {
			nip = io.Atoi(s_nip)
		}
		o.IpsElem, err = shp.GetIps(o.Shp.Type, nip)
		if LogErr(err, "GetIps failed") {
			return nil
		}
		nip = len(o.IpsElem)

		// scratchpad. computed @ each ip
		o.K = la.MatAlloc(o.Nu, o.Nu)
		o.M = la.MatAlloc(o.Nu, o.Nu)
		o.ue = make([]float64, o.Nu)
		o.Rus = make([]float64, o.Nu)

		// scratchpad. computed @ each ip
		o.grav = make([]float64, ndim)
		o.us = make([]float64, ndim)
		o.fi = make([]float64, o.Nu)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o Rod) Id() int { return o.Cid }

// SetEqs set equations
func (o *Rod) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	ndim := Global.Ndim
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Shp.Nverts; m++ {
		for i := 0; i < ndim; i++ {
			r := i + m*ndim
			o.Umap[r] = eqs[m][i]
		}
	}
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *Rod) InterpStarVars(sol *Solution) (ok bool) {

	// skip steady cases
	if Global.Sim.Data.Steady {
		return true
	}

	// TODO
	return true
}

// SetEleConds set element conditions
func (o *Rod) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	if key == "g" {
		o.Gfcn = f
	}
	return true
}

// adds -R to global residual vector fb
func (o Rod) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// for each integration point
	nverts := o.Shp.Nverts
	ndim := Global.Ndim
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// auxiliary
		coef := ip.W
		Jvec := o.Shp.Jvec3d
		G := o.Shp.Gvec
		σ := o.States[idx].Sig

		// update fb with internal forces
		for m := 0; m < nverts; m++ {
			for i := 0; i < ndim; i++ {
				r := o.Umap[i+m*ndim]
				fb[r] -= coef * o.A * σ * G[m] * Jvec[i] // -fi
			}
		}
	}
	return true
}

// adds element K to global Jacobian matrix Kb
func (o Rod) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// zero K matrix
	la.MatFill(o.K, 0)
	la.MatFill(o.M, 0)

	// for each integration point
	nverts := o.Shp.Nverts
	ndim := Global.Ndim
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// auxiliary
		coef := ip.W
		Jvec := o.Shp.Jvec3d
		G := o.Shp.Gvec
		J := o.Shp.J

		// add contribution to consistent tangent matrix
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				for i := 0; i < ndim; i++ {
					for j := 0; j < ndim; j++ {
						r := i + m*ndim
						c := j + n*ndim
						E, err := o.Model.CalcD(o.States[idx], firstIt)
						if LogErr(err, "AddToKb") {
							return
						}
						o.K[r][c] += coef * o.A * E * G[m] * G[n] * Jvec[i] * Jvec[j] / J
					}
				}
			}
		}
	}

	// add K to sparse matrix Kb
	for i, I := range o.Umap {
		for j, J := range o.Umap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return true
}

// Update perform (tangent) update
func (o *Rod) Update(sol *Solution) (ok bool) {

	// for each integration point
	nverts := o.Shp.Nverts
	ndim := Global.Ndim
	for idx, _ := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// auxiliary
		Jvec := o.Shp.Jvec3d
		G := o.Shp.Gvec
		J := o.Shp.J

		// compute strains
		Δε := 0.0
		for m := 0; m < nverts; m++ {
			for i := 0; i < ndim; i++ {
				r := o.Umap[i+m*ndim]
				Δε += G[m] * Jvec[i] * sol.ΔY[r] / J
			}
		}

		// call model update => update stresses
		if LogErr(o.Model.Update(o.States[idx], 0.0, Δε), "Update") {
			return
		}
	}
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *Rod) InitIvs(sol *Solution) (ok bool) {
	nip := len(o.IpsElem)
	o.States = make([]*msolid.OnedState, nip)
	o.StatesBkp = make([]*msolid.OnedState, nip)
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Model.InitIntVars()
		o.StatesBkp[i] = o.States[i].GetCopy()
	}
	return true
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *Rod) SetIvs(zvars map[string][]float64) (ok bool) {
	return true
}

// BackupIvs create copy of internal variables
func (o *Rod) BackupIvs() (ok bool) {
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return true
}

// RestoreIvs restore internal variables from copies
func (o *Rod) RestoreIvs() (ok bool) {
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o Rod) Encode(enc Encoder) (ok bool) {
	return !LogErr(enc.Encode(o.States), "Encode")
}

// Decode decodes internal variables
func (o Rod) Decode(dec Decoder) (ok bool) {
	return !LogErr(dec.Decode(&o.States), "Decode")
}

// OutIpsData returns data from all integration points for output
func (o Rod) OutIpsData() (data []*OutIpData) {
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *Rod) ipvars(idx int, sol *Solution) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.Shp.CalcAtIp(o.X, o.IpsElem[idx], true), "ipvars") {
		return
	}

	// skip if steady (this must be after CalcAtIp, because callers will need S and G)
	if Global.Sim.Data.Steady {
		return true
	}

	// clear variables
	ndim := Global.Ndim
	for i := 0; i < ndim; i++ {
		o.us[i] = 0
	}

	// recover u-variables @ ip
	for m := 0; m < o.Shp.Nverts; m++ {
		for i := 0; i < ndim; i++ {
			r := o.Umap[i+m*ndim]
			o.us[i] += o.Shp.S[m] * sol.Y[r]
		}
	}
	return true
}
