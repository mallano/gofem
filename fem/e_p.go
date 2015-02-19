// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// ElemP implements an element for transient seepage analyses [1]
//  References:
//   [1] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816,
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type ElemP struct {

	// basic data
	Cell *inp.Cell   // cell
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Ndim int         // space dimension
	Np   int         // total number of unknowns == number of vertices

	// integration points
	IpsElem []*shp.Ipoint // integration points of element
	IpsFace []*shp.Ipoint // integration points corresponding to faces

	// material model and internal variables
	Mdl *mporous.Model

	// problem variables
	Pmap []int // assembly map (location array/element equations)

	// internal variables
	States    []*mporous.State
	StatesBkp []*mporous.State

	// gravity
	Gfcn fun.Func // gravity function

	// natural boundary conditions
	NatBcs []*NaturalBc // natural boundary conditions

	// seepage face
	HasSeep    bool    // indicates if this element has seepage faces
	Vid2seepId []int   // [nverts] maps local vertex id to index in Fmap
	SeepId2vid []int   // [nseep] maps seepage face variable id to local vertex id
	Fmap       []int   // [nseep] map of "fl" variables (seepage face)
	Macaulay   bool    // use discrete ramp function instead of smooth ramp
	βrmp       float64 // coefficient for Sramp

	// local starred variables
	ψl []float64 // [nip] ψl* = β1.p + β2.dpdt

	// scratchpad. computed @ each ip
	g   []float64   // [ndim] gravity vector
	pl  float64     // pl: liquid pressure
	gpl []float64   // [ndim] ∇pl: gradient of liquid pressure
	ρwl []float64   // [ndim] ρl*wl: weighted liquid relative velocity
	tmp []float64   // [ndim] temporary (auxiliary) vector
	K   [][]float64 // [np][np] K := dRpl/dpl consistent tangent matrix
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	iallocators["p"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info {

		// new info
		var info Info

		// number of nodes in element
		cell := msh.Cells[cid]
		nverts := cell.Shp.Nverts

		// solution variables
		ykeys := []string{"pl"}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// vertices on seepage faces
		if len(cell.SeepVerts) > 0 {
			ykeys = append(ykeys, "fl")
		}
		for m, _ := range cell.SeepVerts {
			info.Dofs[m] = append(info.Dofs[m], "fl")
		}

		// maps
		info.Y2F = map[string]string{"pl": "ql"}

		// t1 and t2 variables
		info.T1vars = ykeys
		return &info
	}

	// element allocator
	eallocators["p"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {

		// basic data
		var o ElemP
		o.Cell = msh.Cells[cid]
		o.X = BuildCoordsMatrix(o.Cell, msh)
		o.Ndim = msh.Ndim
		o.Np = o.Cell.Shp.Nverts

		// integration points
		var nip, nipf int
		if s_nip, found := io.Keycode(edat.Extra, "nip"); found {
			nip = io.Atoi(s_nip)
		}
		if s_nipf, found := io.Keycode(edat.Extra, "nipf"); found {
			nipf = io.Atoi(s_nipf)
		}
		var err error
		o.IpsElem, err = shp.GetIps(o.Cell.Shp.Type, nip)
		if LogErr(err, "GetIps failed for solid element") {
			return nil
		}
		o.IpsFace, err = shp.GetIps(o.Cell.Shp.FaceType, nipf)
		if LogErr(err, "GetIps failed for face") {
			return nil
		}
		nip = len(o.IpsElem)
		nipf = len(o.IpsFace)

		// materials
		matname := edat.Mat
		cndmat, lrmmat, pormat, err := Global.Mdb.GroupGet3(matname, "c", "l", "p")
		if LogErr(err, "Mdb.GroupGet3 failed") {
			return nil
		}

		// models
		simfnk := Global.Sim.Data.FnameKey
		getnew := false
		cnd := mconduct.GetModel(simfnk, cndmat.Name, cndmat.Model, getnew)
		if LogErrCond(cnd == nil, "cannot allocate conductivity models with name=%q", cndmat.Model) {
			return nil
		}
		lrm := mreten.GetModel(simfnk, lrmmat.Name, lrmmat.Model, getnew)
		if LogErrCond(lrm == nil, "cannot allocate liquid retention model with name=%q", lrmmat.Model) {
			return nil
		}
		o.Mdl = mporous.GetModel(simfnk, pormat.Name, getnew)
		if LogErrCond(o.Mdl == nil, "cannot allocate model for porous medium with name=%q", pormat.Name) {
			return nil
		}
		if LogErr(cnd.Init(cndmat.Prms), "cannot initialise conductivity model") {
			return nil
		}
		if LogErr(lrm.Init(lrmmat.Prms), "cannot initialise liquid retention model") {
			return nil
		}
		if LogErr(o.Mdl.Init(pormat.Prms, cnd, lrm), "cannot initialise porous model") {
			return nil
		}

		// local starred variables
		o.ψl = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.g = make([]float64, o.Ndim)
		o.gpl = make([]float64, o.Ndim)
		o.ρwl = make([]float64, o.Ndim)
		o.tmp = make([]float64, o.Ndim)
		o.K = la.MatAlloc(o.Np, o.Np)

		// seepage face
		nseep := len(o.Cell.SeepVerts)
		if nseep > 0 {
			o.HasSeep = true
			o.SeepId2vid = utl.IntBoolMapSort(o.Cell.SeepVerts)
			o.Vid2seepId = utl.IntVals(o.Cell.Shp.Nverts, -1)
			o.Fmap = make([]int, nseep)
			for i, m := range o.SeepId2vid {
				o.Vid2seepId[m] = i
			}

			// use macaulay function ?
			if s_mac, found := io.Keycode(edat.Extra, "mac"); found {
				o.Macaulay = io.Atob(s_mac)
			}

			// coefficient for smooth ramp function
			o.βrmp = 70.0
			if s_bet, found := io.Keycode(edat.Extra, "bet"); found {
				o.βrmp = io.Atof(s_bet)
			}
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs sets equations
func (o *ElemP) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	o.Pmap = make([]int, o.Np)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		o.Pmap[m] = eqs[m][0]
	}
	if o.HasSeep {
		for i, m := range o.SeepId2vid {
			o.Fmap[i] = eqs[m][1]
		}
	}
	return true
}

// SetEleConds sets element conditions
func (o *ElemP) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return true
}

// SetSurfLoads sets surface loads (natural boundary conditions)
func (o *ElemP) SetNatBcs(key string, idxface int, f fun.Func, extra string) (ok bool) {
	o.NatBcs = append(o.NatBcs, &NaturalBc{key, idxface, f, extra})
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemP) InterpStarVars(sol *Solution) (ok bool) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Cell.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		// interpolate starred variables
		o.ψl[idx] = 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			o.ψl[idx] += o.Cell.Shp.S[m] * sol.Psi[o.Pmap[m]]
		}
	}
	return true
}

// AddToRhs adds -R to global residual vector fb
func (o ElemP) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// for each integration point
	β1 := Global.DynCoefs.β1
	nverts := o.Cell.Shp.Nverts
	var coef, plt, klr, RhoL, Cpl float64
	var err error
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef = o.Cell.Shp.J * ip.W
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		plt = β1*o.pl - o.ψl[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].Sl)
		RhoL = o.States[idx].RhoL
		Cpl, err = o.States[idx].Lvars(o.Mdl)
		if LogErr(err, "calc of tpm variables failed") {
			return
		}

		// compute ρwl
		// see Eq. (6) of [1]
		for i := 0; i < o.Ndim; i++ {
			o.ρwl[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.ρwl[i] += klr * o.Mdl.Klsat[i][j] * (RhoL*o.g[i] - o.gpl[i])
			}
		}

		// add negative of residual term to fb
		// see Eqs. (12) and (17) of [1]
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			fb[r] -= coef * S[m] * Cpl * plt
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.ρwl[i] // += coef * div(ρl*wl)
			}
		}
	}

	// external 'fluxes'
	return o.add_fluxloads_to_rhs(fb, sol)
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o ElemP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// clear matrices
	la.MatFill(o.K, 0)

	// for each integration point
	Cl := o.Mdl.Cl
	β1 := Global.DynCoefs.β1
	nverts := o.Cell.Shp.Nverts
	var coef, plt, klr, RhoL, Cpl, dCpldpl, dklrdpl float64
	var err error
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef = o.Cell.Shp.J * ip.W
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G
		plt = β1*o.pl - o.ψl[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].Sl)
		RhoL = o.States[idx].RhoL
		Cpl, dCpldpl, dklrdpl, err = o.States[idx].Lderivs(o.Mdl)
		if LogErr(err, "calc of tpm derivatives failed") {
			return
		}

		// K := dRpl/dpl
		// see Eqs. (18), (A.2) and (A.3) of [1]
		for n := 0; n < nverts; n++ {
			for j := 0; j < o.Ndim; j++ {
				o.tmp[j] = S[n]*dklrdpl*(RhoL*o.g[j]-o.gpl[j]) + klr*(S[n]*Cl*o.g[j]-G[n][j])
			}
			for m := 0; m < nverts; m++ {
				o.K[m][n] += coef * S[m] * S[n] * (dCpldpl*plt + β1*Cpl)
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.K[m][n] -= coef * G[m][i] * o.Mdl.Klsat[i][j] * o.tmp[j]
					}
				}
			}
		}
	}

	// add to sparse matrix Kb
	for i, I := range o.Pmap {
		for j, J := range o.Pmap {
			Kb.Put(I, J, o.K[i][j])
		}
	}
	return true
}

// Update performs (tangent) update
func (o *ElemP) Update(sol *Solution) (ok bool) {

	// for each integration point
	var Δpl float64
	for idx, _ := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], false), "Update") {
			return
		}

		// compute Δpl @ ip by means of interpolating from nodes
		Δpl = 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			r := o.Pmap[m]
			Δpl += o.Cell.Shp.S[m] * sol.ΔY[r]
		}

		// update state
		if LogErr(o.Mdl.Update(o.States[idx], Δpl, 0, 0), "update failed") {
			return
		}
	}
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs resets (and fix) internal variables after primary variables have been changed
func (o *ElemP) InitIvs(sol *Solution) (ok bool) {

	// for each integration point
	nip := len(o.IpsElem)
	o.States = make([]*mporous.State, nip)
	o.StatesBkp = make([]*mporous.State, nip)
	var err error
	for idx, _ := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], false), "InitIvs") {
			return
		}

		// compute pl @ ip by means of interpolating from nodes
		o.pl = 0
		for m := 0; m < o.Cell.Shp.Nverts; m++ {
			r := o.Pmap[m]
			o.pl += o.Cell.Shp.S[m] * sol.Y[r]
		}

		// state initialisation
		o.States[idx], err = o.Mdl.NewState(o.pl, 0, 0)
		if LogErr(err, "state initialisation failed") {
			return
		}

		// backup copy
		o.StatesBkp[idx] = o.States[idx].GetCopy()
	}
	return true
}

// SetIvs sets secondary variables; e.g. during initialisation via files
func (o *ElemP) SetIvs(zvars map[string][]float64) (ok bool) {
	return true
}

// BackupIvs creates copy of internal variables
func (o *ElemP) BackupIvs() (ok bool) {
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return true
}

// RestoreIvs restores internal variables from copies
func (o *ElemP) RestoreIvs() (ok bool) {
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElemP) Encode(enc Encoder) (ok bool) {
	return !LogErr(enc.Encode(o.States), "Encode")
}

// Decode decodes internal variables
func (o ElemP) Decode(dec Decoder) (ok bool) {
	return !LogErr(dec.Decode(&o.States), "Decode")
}

// OutIpsData returns data from all integration points for output
func (o ElemP) OutIpsData() (labels []string, data []*OutIpData) {
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemP) ipvars(idx int, sol *Solution) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true), "ipvars") {
		return
	}

	// gravity
	o.g[o.Ndim-1] = 0
	if o.Gfcn != nil {
		o.g[o.Ndim-1] = -o.Gfcn.F(sol.T, nil)
	}

	// clear pl and its gradient @ ip
	o.pl = 0
	for i := 0; i < o.Ndim; i++ {
		o.gpl[i] = 0
	}

	// compute pl and its gradient @ ip by means of interpolating from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		r := o.Pmap[m]
		o.pl += o.Cell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i] += o.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}
	return true
}

// add_fluxloads_to_rhs adds surfaces loads to rhs
func (o ElemP) add_fluxloads_to_rhs(fb []float64, sol *Solution) (ok bool) {
	return true
}

// ramp implements the ramp function
func (o *ElemP) ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.βrmp)
}

// rampderiv returns the ramp function first derivative
func (o *ElemP) rampderiv(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.βrmp)
}
