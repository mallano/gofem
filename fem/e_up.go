// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// ElemUP represents an element for porous media based on the u-p formulation [1]
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
type ElemUP struct {

	// underlying elements
	U *ElemU // u-element
	P *ElemP // p-element

	// scratchpad. computed @ each ip
	divvs float64   // divergence of velocity of solids
	bs    []float64 // auxiliary vector
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["up"] = func(ndim int, cellType string, faceConds []*FaceCond) *Info {

		// new info
		var info Info

		// p-element cell type
		p_cellType := cellType
		lbb := !Global.Sim.Data.NoLBB
		if lbb {
			p_cellType = shp.GetBasicType(cellType)
		}

		// underlying cells info
		u_info := infogetters["u"](ndim, cellType, faceConds)
		p_info := infogetters["p"](ndim, p_cellType, faceConds)

		// solution variables
		nverts := shp.GetNverts(cellType)
		info.Dofs = make([][]string, nverts)
		for i, dofs := range u_info.Dofs {
			info.Dofs[i] = append(info.Dofs[i], dofs...)
		}
		for i, dofs := range p_info.Dofs {
			info.Dofs[i] = append(info.Dofs[i], dofs...)
		}

		// maps
		info.Y2F = u_info.Y2F
		for key, val := range p_info.Y2F {
			info.Y2F[key] = val
		}

		// t1 and t2 variables
		info.T1vars = p_info.T1vars
		info.T2vars = u_info.T2vars
		return &info
	}

	// element allocator
	eallocators["up"] = func(ndim int, cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemUP

		// p-element cell type
		p_cellType := cellType
		lbb := !Global.Sim.Data.NoLBB
		if lbb {
			p_cellType = shp.GetBasicType(cellType)
		}

		// underlying elements
		u_allocator := eallocators["u"]
		p_allocator := eallocators["p"]
		u_elem := u_allocator(ndim, cellType, faceConds, cid, edat, x)
		p_elem := p_allocator(ndim, p_cellType, faceConds, cid, edat, x)
		if LogErrCond(u_elem == nil, "cannot allocate underlying u-element") {
			return nil
		}
		if LogErrCond(p_elem == nil, "cannot allocate underlying p-element") {
			return nil
		}
		o.U = u_elem.(*ElemU)
		o.P = p_elem.(*ElemP)

		// scratchpad. computed @ each ip
		o.bs = make([]float64, ndim)

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o ElemUP) Id() int { return o.U.Id() }

// SetEqs set equations
func (o *ElemUP) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	ndim := o.U.Ndim
	eqs_u := make([][]int, o.U.Shp.Nverts)
	eqs_p := make([][]int, o.P.Shp.Nverts)
	for m := 0; m < o.U.Shp.Nverts; m++ {
		eqs_u[m] = eqs[m][:ndim]
	}
	idxp := ndim
	for m := 0; m < o.P.Shp.Nverts; m++ {
		eqs_p[m] = []int{eqs[m][idxp]}
	}
	if !o.U.SetEqs(eqs_u, mixedform_eqs) {
		return
	}
	return o.P.SetEqs(eqs_p, nil)
}

// SetEleConds set element conditions
func (o *ElemUP) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	if !o.U.SetEleConds(key, f, extra) {
		return
	}
	return o.P.SetEleConds(key, f, extra)
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemUP) InterpStarVars(sol *Solution) (ok bool) {
	if !o.U.InterpStarVars(sol) {
		return
	}
	return o.P.InterpStarVars(sol)
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
	if !o.U.Update(sol) {
		return
	}
	return o.P.Update(sol)
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *ElemUP) InitIvs(sol *Solution) (ok bool) {
	if !o.U.InitIvs(sol) {
		return
	}
	return o.P.InitIvs(sol)
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *ElemUP) SetIvs(zvars map[string][]float64) (ok bool) {
	if !o.U.SetIvs(zvars) {
		return
	}
	return o.P.SetIvs(zvars)
}

// BackupIvs create copy of internal variables
func (o *ElemUP) BackupIvs() (ok bool) {
	if !o.U.BackupIvs() {
		return
	}
	return o.P.BackupIvs()
}

// RestoreIvs restore internal variables from copies
func (o *ElemUP) RestoreIvs() (ok bool) {
	if !o.U.RestoreIvs() {
		return
	}
	return o.P.RestoreIvs()
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElemUP) Encode(enc Encoder) (ok bool) {
	if !o.U.Encode(enc) {
		return
	}
	return o.P.Encode(enc)
}

// Decode decodes internal variables
func (o ElemUP) Decode(dec Decoder) (ok bool) {
	if !o.U.Decode(dec) {
		return
	}
	return o.P.Decode(dec)
}

// OutIpsData returns data from all integration points for output
func (o ElemUP) OutIpsData() (data []*OutIpData) {
	u_dat := o.U.OutIpsData()
	p_dat := o.P.OutIpsData()
	nip := len(o.U.IpsElem)
	chk.IntAssert(len(u_dat), nip)
	chk.IntAssert(len(u_dat), len(p_dat))
	data = make([]*OutIpData, nip)
	for i, d := range u_dat {
		for key, val := range p_dat[i].V {
			d.V[key] = val
		}
		data[i] = &OutIpData{d.Eid, d.X, d.V}
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemUP) ipvars(idx int, sol *Solution, tpm_derivs bool) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.U.Shp.CalcAtIp(o.U.X, o.U.IpsElem[idx], true), "ipvars") {
		return
	}
	if LogErr(o.P.Shp.CalcAtIp(o.P.X, o.U.IpsElem[idx], true), "ipvars") {
		return
	}

	// auxiliary
	ndim := o.U.Ndim

	// gravity
	o.U.grav[ndim-1] = 0
	if o.U.Gfcn != nil {
		o.U.grav[ndim-1] = -o.U.Gfcn.F(sol.T, nil)
	}

	// clear variables
	o.P.pl = 0
	for i := 0; i < ndim; i++ {
		o.P.gpl[i] = 0
		o.U.us[i] = 0
	}

	// recover p-variables @ ip
	for m := 0; m < o.P.Shp.Nverts; m++ {
		r := o.P.Pmap[m]
		o.P.pl += o.P.Shp.S[m] * sol.Y[r]
		for i := 0; i < ndim; i++ {
			o.P.gpl[i] += o.P.Shp.G[m][i] * sol.Y[r]
		}
	}

	// recover u-variables @ ip
	var divus float64
	for m := 0; m < o.U.Shp.Nverts; m++ {
		for i := 0; i < ndim; i++ {
			r := o.U.Umap[i+m*ndim]
			o.U.us[i] += o.U.Shp.S[m] * sol.Y[r]
			divus += o.U.Shp.G[m][i] * sol.Y[r]
		}
	}
	return true
}
