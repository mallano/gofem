// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// ElemP represents a solid element with p (e.g. pressure) as primary variable
type ElemP struct {

	// basic data
	Cell *inp.Cell   // cell
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Ndim int         // space dimension
	Np   int         // total number of unknowns == number of vertices

	// gravity
	Gfcn fun.Func // gravity function

	// integration points
	IpsElem []*shp.Ipoint // integration points of element
	IpsFace []*shp.Ipoint // integration points corresponding to faces

	// material model and internal variables
	Model *mporous.Model

	// problem variables
	Pmap []int // assembly map (location array/element equations)

	// internal variables
	States    []*mporous.StateLG
	StatesBkp []*mporous.StateLG

	// natural boundary conditions
	NatBcs []*NaturalBc //natural boundary conditions

	// local starred variables
	ψl []float64 // [nip] ψl* = β1.p + β2.dpdt

	// scratchpad. computed @ each ip
	grav    []float64   // [ndim] gravity vector
	dpldt   float64     // dpl/dt: liquid pressure rate
	pl      float64     // pl: liquid pressure
	gpl     []float64   // [ndim] ∇pl: gradient of liquid pressure
	ρwl     []float64   // [ndim] ρl*wl: weighted liquid relative velocity
	dRplDpl [][]float64 // [np][np] consistent tangent matrix
	hl      []float64   // [ndim] auxiliary vector
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
		if s_nip, found := utl.Keycode(edat.Extra, "nip"); found {
			nip = utl.Atoi(s_nip)
		}
		if s_nipf, found := utl.Keycode(edat.Extra, "nipf"); found {
			nipf = utl.Atoi(s_nipf)
		}
		o.IpsElem = shp.GetIps(o.Cell.Shp.Type, nip)
		o.IpsFace = shp.GetIps(o.Cell.Shp.FaceType, nipf)
		nip = len(o.IpsElem)
		nipf = len(o.IpsFace)

		// material model
		matname := edat.Mat
		matdata := global.Mdb.Get(matname)
		if LogErrCond(matdata == nil, "Mdb.Get failed\n") {
			return nil
		}
		o.Model = new(mporous.Model)
		o.Model.Init(matdata.Prms)

		// allocate states
		o.States = make([]*mporous.StateLG, nip)
		o.StatesBkp = make([]*mporous.StateLG, nip)
		for i := 0; i < nip; i++ {
			o.States[i] = new(mporous.StateLG)
			o.StatesBkp[i] = new(mporous.StateLG)
		}

		// local starred variables
		o.ψl = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.grav = make([]float64, o.Ndim)
		o.gpl = make([]float64, o.Ndim)
		o.ρwl = make([]float64, o.Ndim)
		o.dRplDpl = la.MatAlloc(o.Np, o.Np)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *ElemP) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	o.Pmap = make([]int, o.Np)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		o.Pmap[m] = eqs[m][0]
	}
	return true
}

// SetEleConds set element conditions
func (o *ElemP) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return true
}

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *ElemP) SetSurfLoads(key string, idxface int, f fun.Func, extra string) (ok bool) {
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

// adds -R to global residual vector fb
func (o ElemP) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// for each integration point
	tpm := mporous.TPM
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol, false) {
			return
		}
		coef := o.Cell.Shp.J * ip.W
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// add residual term to fb
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			fb[r] -= coef * S[m] * tpm.Cpl * o.dpldt
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.ρwl[i]
			}
		}
	}

	// external 'fluxes'
	return o.add_fluxloads_to_rhs(fb, sol)
}

// adds element K to global Jacobian matrix Kb
func (o ElemP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// clear matrices
	la.MatFill(o.dRplDpl, 0)

	// for each integration point
	dc := global.DynCoefs
	tpm := mporous.TPM
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol, true) {
			return
		}
		coef := o.Cell.Shp.J * ip.W
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// dRplDpl
		for m := 0; m < nverts; m++ {
			for n := 0; n < nverts; n++ {
				o.dRplDpl[m][n] += coef * S[m] * S[n] * (tpm.DCplDpl*o.dpldt + dc.β1*tpm.Cpl)
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.dRplDpl[m][n] -= coef * G[m][i] * o.Model.Klsat[i][j] * (S[n]*tpm.DklrDpl*o.hl[j] - tpm.Klr*(S[n]*tpm.Cl*(-o.grav[j])+G[n][j]))
					}
				}
			}
		}
	}
	return true
}

// Update perform (tangent) update
func (o *ElemP) Update(sol *Solution) (ok bool) {
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *ElemP) InitIvs(sol *Solution) (ok bool) {
	return true
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *ElemP) SetIvs(zvars map[string][]float64) (ok bool) {
	return true
}

// BackupIvs create copy of internal variables
func (o *ElemP) BackupIvs() (ok bool) {
	return true
}

// RestoreIvs restore internal variables from copies
func (o *ElemP) RestoreIvs() (ok bool) {
	return true
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// add_fluxloads_to_rhs adds surfaces loads to rhs
func (o ElemP) add_fluxloads_to_rhs(fb []float64, sol *Solution) (ok bool) {
	return true
}

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemP) ipvars(idx int, sol *Solution, tpm_derivs bool) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true), "ipvars") {
		return
	}

	// gravity
	o.grav[o.Ndim-1] = 0
	if o.Gfcn != nil {
		o.grav[o.Ndim-1] = -o.Gfcn.F(sol.T, nil)
	}

	// clear pl and its gradient @ ip
	o.pl = 0
	for i := 0; i < o.Ndim; i++ {
		o.gpl[i] = 0
	}

	// compute pl and gradient of pl @ ip by means of interpolation from nodes
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		r := o.Pmap[m]
		o.pl += o.Cell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i] += o.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}

	// TPM variables
	if LogErr(mporous.CalcL(o.pl, o.States[idx], o.Model, tpm_derivs), "ipvars") {
		return
	}

	// local variables
	o.dpldt = global.DynCoefs.β1*o.pl - o.ψl[idx]

	// auxiliary vectors
	for i := 0; i < o.Ndim; i++ {
		o.hl[i] = o.States[idx].RhoL*o.grav[i] - o.gpl[i]
		o.ρwl[i] = 0
		for j := 0; j < o.Ndim; j++ {
			o.ρwl[i] += mporous.TPM.Klr * o.Model.Klsat[i][j] * o.hl[j]
		}
	}
	return true
}
