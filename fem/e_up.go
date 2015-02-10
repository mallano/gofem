// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
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
	iallocators["up"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) (info Info) {

		// cell
		cell := msh.Cells[cid]

		// lbb flag
		lbb := true
		if val, found := utl.Keycode(edat.Extra, "nolbb"); found {
			lbb = !utl.Atob(val)
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
		return
	}

	// element allocator
	eallocators["up"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {

		// basic data
		var o ElemUP
		o.Ndim = msh.Ndim

		// flags
		o.Lbb = true
		if val, found := utl.Keycode(edat.Extra, "nolbb"); found {
			o.Lbb = !utl.Atob(val)
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
func (o *ElemUP) SetEqs(eqs [][]int, mixedform_eqs []int) {
	eqs_u := make([][]int, o.u.Cell.Shp.Nverts)
	eqs_p := make([][]int, o.p.Cell.Shp.Nverts)
	for m := 0; m < o.u.Cell.Shp.Nverts; m++ {
		eqs_u[m] = eqs[m][:o.Ndim]
	}
	idxp := o.Ndim
	for m := 0; m < o.p.Cell.Shp.Nverts; m++ {
		eqs_p[m] = []int{eqs[m][idxp]}
	}
	o.u.SetEqs(eqs_u, mixedform_eqs)
	o.p.SetEqs(eqs_p, nil)
}

// SetEleConds set element conditions
func (o *ElemUP) SetEleConds(key string, f fun.Func, extra string) {
	o.u.SetEleConds(key, f, extra)
	o.p.SetEleConds(key, f, extra)
}

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *ElemUP) SetSurfLoads(key string, idxface int, f fun.Func, extra string) {
	u_surfkeys := o.u.surfloads_keys()
	if u_surfkeys[key] {
		o.u.SetSurfLoads(key, idxface, f, extra)
		return
	}
	o.p.SetSurfLoads(key, idxface, f, extra)
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemUP) InterpStarVars(sol *Solution) (err error) {
	err = o.u.InterpStarVars(sol)
	if err != nil {
		return
	}
	return o.p.InterpStarVars(sol)
}

// adds -R to global residual vector fb
func (o ElemUP) AddToRhs(fb []float64, sol *Solution) (err error) {
	return
}

// adds element K to global Jacobian matrix Kb
func (o ElemUP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (err error) {
	// success
	return
}

// Update perform (tangent) update
func (o *ElemUP) Update(sol *Solution) (err error) {
	err = o.u.Update(sol)
	if err != nil {
		return
	}
	return o.p.Update(sol)
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *ElemUP) InitIvs(sol *Solution) {
	o.u.InitIvs(sol)
	o.p.InitIvs(sol)
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *ElemUP) SetIvs(zvars map[string][]float64) {
	o.u.SetIvs(zvars)
	o.p.SetIvs(zvars)
}

// BackupIvs create copy of internal variables
func (o *ElemUP) BackupIvs() {
	o.u.BackupIvs()
	o.p.BackupIvs()
}

// RestoreIvs restore internal variables from copies
func (o *ElemUP) RestoreIvs() {
	o.u.RestoreIvs()
	o.p.RestoreIvs()
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// log logs errors
func (o *ElemUP) log(err error, msg string) error {
	if err != nil {
		log.Printf("ElemUP: eid=%d %s failed with %v\n", o.u.Cell.Id, msg, err)
	}
	return err
}

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemUP) ipvars(idx int, sol *Solution, tpm_derivs bool) (err error) {

	// interpolation functions and gradients
	err = o.log(o.u.Cell.Shp.CalcAtIp(o.u.X, o.u.IpsElem[idx], true), "ipvars")
	if err != nil {
		return
	}
	err = o.log(o.p.Cell.Shp.CalcAtIp(o.p.X, o.u.IpsElem[idx], true), "ipvars")
	if err != nil {
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

	// TPM variables
	err = o.log(mporous.CalcLS(o.p.pl, divus, o.p.States[idx], o.p.Model, tpm_derivs), "ipvars")
	if err != nil {
		return
	}

	// local variables
	dc := global.DynCoefs
	o.p.dpldt = dc.β1*o.p.pl - o.p.ψl[idx]
	o.divvs = dc.α4*divus - o.u.divχs[idx]

	// auxiliary vectors
	for i := 0; i < o.Ndim; i++ {
		o.bs[i] = dc.α1*o.u.us[i] - o.u.ζs[idx][i] - o.u.grav[i]
		o.p.hl[i] = -o.p.States[idx].RhoL*o.bs[i] - o.p.gpl[i]
		o.p.ρwl[i] = 0
		for j := 0; j < o.Ndim; j++ {
			o.p.ρwl[i] += mporous.TPM.Klr * o.p.Model.Klsat[i][j] * o.p.hl[j]
		}
	}

	// success
	return
}
