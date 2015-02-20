// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

// ElemUP represents an element for porous media based on the u-p formulation
type ElemUP struct {

	// basic data
	Ndim int // space dimension

	// flags
	Lbb bool // Ladyženskaja-Babuška-Brezzi element

	// basis elemsnts
	u *ElemU // u-element
	p *ElemP // p-element

	// scratchpad. computed @ each ip
	divvs float64   // divergence of velocity of solids
	bs    []float64 // auxiliary vector
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	iallocators["up"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info {

		// new info
		var info Info

		// cell
		cell := msh.Cells[cid]

		// lbb flag
		lbb := true
		if val, found := io.Keycode(edat.Extra, "nolbb"); found {
			lbb = !io.Atob(val)
		}

		// number of u and p nodes
		nne_u := cell.Shp.Nverts
		nne_p := nne_u
		if lbb {
			shp_p := shp.Get(cell.Shp.BasicType)
			nne_p = shp_p.Nverts
		}

		// solution variables
		ykeys := []string{"ux", "uy", "pl"}
		if msh.Ndim == 3 {
			ykeys = []string{"ux", "uy", "uz", "pl"}
		}
		ukeys := ykeys[:msh.Ndim]
		info.Dofs = make([][]string, nne_u)
		for m := 0; m < nne_p; m++ {
			info.Dofs[m] = ykeys
		}
		for m := nne_p; m < nne_u; m++ {
			info.Dofs[m] = ukeys
		}

		// maps
		info.Y2F = map[string]string{"ux": "fx", "uy": "fy", "uz": "fz", "pl": "ql"}

		// t1 and t2 variables
		info.T1vars = []string{"pl"}
		info.T2vars = ukeys
		return &info
	}

	// element allocator
	eallocators["up"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {

		// basic data
		var o ElemUP
		o.Ndim = msh.Ndim

		// flags
		o.Lbb = true
		if val, found := io.Keycode(edat.Extra, "nolbb"); found {
			o.Lbb = !io.Atob(val)
		}

		// basis elements
		alloc_u := eallocators["u"]
		alloc_p := eallocators["p"]
		o.u = alloc_u(edat, cid, msh).(*ElemU)
		o.p = alloc_p(edat, cid, msh).(*ElemP)

		// scratchpad. computed @ each ip
		o.bs = make([]float64, o.Ndim)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *ElemUP) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	eqs_u := make([][]int, o.u.Cell.Shp.Nverts)
	eqs_p := make([][]int, o.p.Cell.Shp.Nverts)
	for m := 0; m < o.u.Cell.Shp.Nverts; m++ {
		eqs_u[m] = eqs[m][:o.Ndim]
	}
	idxp := o.Ndim
	for m := 0; m < o.p.Cell.Shp.Nverts; m++ {
		eqs_p[m] = []int{eqs[m][idxp]}
	}
	if !o.u.SetEqs(eqs_u, mixedform_eqs) {
		return
	}
	return o.p.SetEqs(eqs_p, nil)
}

// SetEleConds set element conditions
func (o *ElemUP) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	if !o.u.SetEleConds(key, f, extra) {
		return
	}
	return o.p.SetEleConds(key, f, extra)
}

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *ElemUP) SetNatBcs(key string, idxface int, f fun.Func, extra string) (ok bool) {
	u_surfkeys := o.u.surfloads_keys()
	if u_surfkeys[key] {
		if !o.u.SetNatBcs(key, idxface, f, extra) {
			return
		}
		return true
	}
	return o.p.SetNatBcs(key, idxface, f, extra)
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemUP) InterpStarVars(sol *Solution) (ok bool) {
	if !o.u.InterpStarVars(sol) {
		return
	}
	return o.p.InterpStarVars(sol)
}

// adds -R to global residual vector fb
func (o ElemUP) AddToRhs(fb []float64, sol *Solution) (ok bool) {
	return true
}

// adds element K to global Jacobian matrix Kb
func (o ElemUP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {
	return true
}

// Update perform (tangent) update
func (o *ElemUP) Update(sol *Solution) (ok bool) {
	if !o.u.Update(sol) {
		return
	}
	return o.p.Update(sol)
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *ElemUP) InitIvs(sol *Solution) (ok bool) {
	if !o.u.InitIvs(sol) {
		return
	}
	return o.p.InitIvs(sol)
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *ElemUP) SetIvs(zvars map[string][]float64) (ok bool) {
	if !o.u.SetIvs(zvars) {
		return
	}
	return o.p.SetIvs(zvars)
}

// BackupIvs create copy of internal variables
func (o *ElemUP) BackupIvs() (ok bool) {
	if !o.u.BackupIvs() {
		return
	}
	return o.p.BackupIvs()
}

// RestoreIvs restore internal variables from copies
func (o *ElemUP) RestoreIvs() (ok bool) {
	if !o.u.RestoreIvs() {
		return
	}
	return o.p.RestoreIvs()
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElemUP) Encode(enc Encoder) (ok bool) {
	if !o.u.Encode(enc) {
		return
	}
	return o.p.Encode(enc)
}

// Decode decodes internal variables
func (o ElemUP) Decode(dec Decoder) (ok bool) {
	if !o.u.Decode(dec) {
		return
	}
	return o.p.Decode(dec)
}

// OutIpsData returns data from all integration points for output
func (o ElemUP) OutIpsData() (labels []string, data []*OutIpData) {
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemUP) ipvars(idx int, sol *Solution, tpm_derivs bool) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.u.Cell.Shp.CalcAtIp(o.u.X, o.u.IpsElem[idx], true), "ipvars") {
		return
	}
	if LogErr(o.p.Cell.Shp.CalcAtIp(o.p.X, o.u.IpsElem[idx], true), "ipvars") {
		return
	}

	// gravity
	o.u.grav[o.Ndim-1] = 0
	if o.u.Gfcn != nil {
		o.u.grav[o.Ndim-1] = -o.u.Gfcn.F(sol.T, nil)
	}

	// clear variables
	o.p.pl = 0
	for i := 0; i < o.Ndim; i++ {
		o.p.gpl[i] = 0
		o.u.us[i] = 0
	}

	// recover p-variables @ ip
	for m := 0; m < o.p.Cell.Shp.Nverts; m++ {
		r := o.p.Pmap[m]
		o.p.pl += o.p.Cell.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.p.gpl[i] += o.p.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}

	// recover u-variables @ ip
	var divus float64
	for m := 0; m < o.u.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := o.u.Umap[i+m*o.Ndim]
			o.u.us[i] += o.u.Cell.Shp.S[m] * sol.Y[r]
			divus += o.u.Cell.Shp.G[m][i] * sol.Y[r]
		}
	}
	return true
}
