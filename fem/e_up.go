// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
)

// ElemUP represents an element for porous media based on the u-p formulation [1]
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816,
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
type ElemUP struct {

	// auxiliary
	Ndim   int         // space dimension
	Fconds []*FaceCond // face conditions; e.g. seepage faces
	CtypeU string      // u: cell type
	CtypeP string      // p: cell type

	// underlying elements
	U *ElemU // u-element
	P *ElemP // p-element

	// scratchpad. computed @ each ip
	divvs float64     // divvs = div(α4・us - χs) = α4・div(us) - div(χs); (see Eq. 35a [1]) divergence of velocity of solids
	bs    []float64   // bs = as - g = α1・u - ζs - g; (Eqs 35b and A.1 [1]) with 'as' being the acceleration of solids and g, gravity
	hl    []float64   // hl = -ρL・bs - ∇pl; Eq (A.1) of [1]
	Kup   [][]float64 // [nu][np] Kup := dRus/dpl consistent tangent matrix
	Kpu   [][]float64 // [np][nu] Kpu := dRpl/dus consistent tangent matrix

	// for seepage face derivatives
	dρldus_ex [][]float64 // [nverts][nverts*ndim] ∂ρl/∂us extrapolted to nodes => if has qb (flux)
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
		o.Ndim = ndim
		o.Fconds = faceConds

		// p-element cell type
		p_cellType := cellType
		lbb := !Global.Sim.Data.NoLBB
		if lbb {
			p_cellType = shp.GetBasicType(cellType)
		}

		// cell types
		o.CtypeU = cellType
		o.CtypeP = p_cellType

		// allocate u element
		u_allocator := eallocators["u"]
		u_elem := u_allocator(ndim, cellType, faceConds, cid, edat, x)
		if LogErrCond(u_elem == nil, "cannot allocate underlying u-element") {
			return nil
		}
		o.U = u_elem.(*ElemU)

		// make sure p-element uses the same nubmer of integration points than u-element
		edat.Nip = len(o.U.IpsElem)
		//edat.Nipf = len(o.U.IpsFace) // TODO: check if this is necessary

		// allocate p-element
		p_allocator := eallocators["p"]
		p_elem := p_allocator(ndim, p_cellType, faceConds, cid, edat, x)
		if LogErrCond(p_elem == nil, "cannot allocate underlying p-element") {
			return nil
		}
		o.P = p_elem.(*ElemP)

		// scratchpad. computed @ each ip
		o.bs = make([]float64, ndim)
		o.hl = make([]float64, ndim)
		o.Kup = la.MatAlloc(o.U.Nu, o.P.Np)
		o.Kpu = la.MatAlloc(o.P.Np, o.U.Nu)

		// seepage terms
		if o.P.DoExtrap {
			p_nverts := o.P.Shp.Nverts
			u_nverts := o.U.Shp.Nverts
			o.dρldus_ex = la.MatAlloc(p_nverts, u_nverts*ndim)
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o ElemUP) Id() int { return o.U.Id() }

// SetEqs set equations
func (o *ElemUP) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {

	// u: equations
	u_getter := infogetters["u"]
	u_info := u_getter(o.Ndim, o.CtypeU, o.Fconds)
	u_nverts := len(u_info.Dofs)
	u_eqs := make([][]int, u_nverts)
	for i := 0; i < u_nverts; i++ {
		nkeys := len(u_info.Dofs[i])
		u_eqs[i] = make([]int, nkeys)
		for j := 0; j < nkeys; j++ {
			u_eqs[i][j] = eqs[i][j]
		}
	}

	// p: equations
	p_getter := infogetters["p"]
	p_info := p_getter(o.Ndim, o.CtypeP, o.Fconds)
	p_nverts := len(p_info.Dofs)
	p_eqs := make([][]int, p_nverts)
	for i := 0; i < p_nverts; i++ {
		start := len(u_info.Dofs[i])
		nkeys := len(p_info.Dofs[i])
		p_eqs[i] = make([]int, nkeys)
		for j := 0; j < nkeys; j++ {
			p_eqs[i][j] = eqs[i][start+j]
		}
	}

	// set equations
	if !o.U.SetEqs(u_eqs, mixedform_eqs) {
		return
	}
	return o.P.SetEqs(p_eqs, nil)
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

	// for each integration point
	ndim := o.Ndim
	u_nverts := o.U.Shp.Nverts
	p_nverts := o.P.Shp.Nverts
	var r int
	for idx, ip := range o.U.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.P.Shp.CalcAtIp(o.P.X, ip, true), "InterpStarVars") {
			return
		}
		if LogErr(o.U.Shp.CalcAtIp(o.U.X, ip, true), "InterpStarVars") {
			return
		}
		S := o.U.Shp.S
		G := o.U.Shp.G
		Sb := o.P.Shp.S

		// clear local variables
		o.P.ψl[idx], o.U.divχs[idx] = 0, 0
		for i := 0; i < ndim; i++ {
			o.U.ζs[idx][i], o.U.χs[idx][i] = 0, 0
		}

		// p-variables
		for m := 0; m < p_nverts; m++ {
			r = o.P.Pmap[m]
			o.P.ψl[idx] += Sb[m] * sol.Psi[r]
		}

		// u-variables
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < ndim; i++ {
				r = o.U.Umap[i+m*ndim]
				o.U.ζs[idx][i] += S[m] * sol.Zet[r]
				o.U.χs[idx][i] += S[m] * sol.Chi[r]
				o.U.divχs[idx] += G[m][i] * sol.Chi[r]
			}
		}
	}
	return true
}

// adds -R to global residual vector fb
func (o ElemUP) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// clear variables
	if o.P.DoExtrap {
		la.VecFill(o.P.ρl_ex, 0)
	}

	// for each integration point
	dc := Global.DynCoefs
	ndim := o.Ndim
	u_nverts := o.U.Shp.Nverts
	p_nverts := o.P.Shp.Nverts
	var coef, plt, klr, ρl, ρ, p, Cpl, Cvs, divus, divvs float64
	var err error
	var r int
	for idx, ip := range o.U.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef = o.U.Shp.J * ip.W
		S := o.U.Shp.S
		G := o.U.Shp.G
		Sb := o.P.Shp.S
		Gb := o.P.Shp.G

		// auxiliary
		σe := o.U.States[idx].Sig
		divus = o.P.States[idx].Divus
		divvs = dc.α4*divus - o.U.divχs[idx] // divergence of Eq. (35a) [1]

		// tpm variables
		plt = dc.β1*o.P.pl - o.P.ψl[idx] // Eq. (35c) [1]
		klr = o.P.Mdl.Cnd.Klr(o.P.States[idx].Sl)
		ρl, ρ, p, Cpl, Cvs, err = o.P.States[idx].LSvars(o.P.Mdl)
		if LogErr(err, "calc of tpm variables failed") {
			return
		}

		// debug
		if false {
			if math.Abs(plt) < 1e-14 {
				plt = 0
			}
			io.Pf("pl=%13.10f plt=%13.10f klr=%13.10f ρl=%13.10f sl=%13.10f ρ=%13.10f Cpl=%13.10f Cvs=%13.10f\n", o.P.pl, plt, klr, ρl, o.P.States[idx].Sl, ρ, Cpl, Cvs)
		}

		// compute ρwl. see Eq (34b) and (35) of [1]
		for i := 0; i < ndim; i++ {
			o.P.ρwl[i] = 0
			for j := 0; j < ndim; j++ {
				o.P.ρwl[i] += klr * o.P.Mdl.Klsat[i][j] * o.hl[j]
			}
		}

		// p: add negative of residual term to fb; see Eqs. (38a) and (45a) of [1]
		for m := 0; m < p_nverts; m++ {
			r = o.P.Pmap[m]
			fb[r] -= coef * Sb[m] * (Cpl*plt + Cvs*divvs)
			for i := 0; i < ndim; i++ {
				fb[r] += coef * Gb[m][i] * o.P.ρwl[i] // += coef * div(ρl*wl)
			}
			if o.P.DoExtrap { // Eq. (19) of [2]
				o.P.ρl_ex[m] += o.P.Emat[m][idx] * ρl
			}
		}

		// u: add negative of residual term to fb; see Eqs. (38b) and (45b) [1]
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < ndim; i++ {
				r = o.U.Umap[i+m*ndim]
				fb[r] -= coef * S[m] * ρ * o.bs[i]
				for j := 0; j < ndim; j++ {
					fb[r] -= coef * tsr.M2T(σe, i, j) * G[m][j]
				}
				fb[r] += coef * p * G[m][i]
			}
		}
	}

	// external forces
	if len(o.U.NatBcs) > 0 {
		if !o.U.add_surfloads_to_rhs(fb, sol) {
			return
		}
	}

	// contribution from natural boundary conditions
	if len(o.P.NatBcs) > 0 {
		return o.P.add_natbcs_to_rhs(fb, sol)
	}
	return true
}

// adds element K to global Jacobian matrix Kb
func (o ElemUP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// clear matrices
	ndim := o.Ndim
	u_nverts := o.U.Shp.Nverts
	p_nverts := o.P.Shp.Nverts
	la.MatFill(o.P.Kpp, 0)
	for i := 0; i < o.U.Nu; i++ {
		for j := 0; j < o.P.Np; j++ {
			o.Kup[i][j] = 0
			o.Kpu[j][i] = 0
		}
		for j := 0; j < o.U.Nu; j++ {
			o.U.K[i][j] = 0
		}
	}
	if o.P.DoExtrap {
		for i := 0; i < p_nverts; i++ {
			o.P.ρl_ex[i] = 0
			for j := 0; j < p_nverts; j++ {
				o.P.dρldpl_ex[i][j] = 0
			}
			for j := 0; j < o.U.Nu; j++ {
				o.dρldus_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	dc := Global.DynCoefs
	var coef, plt, klr, ρL, Cl, divus, divvs float64
	var ρl, ρ, Cpl, Cvs, dρdpl, dpdpl, dCpldpl, dCvsdpl, dklrdpl, dCpldusM, dρldusM, dρdusM float64
	var err error
	var r, c int
	for idx, ip := range o.U.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef = o.U.Shp.J * ip.W
		S := o.U.Shp.S
		G := o.U.Shp.G
		Sb := o.P.Shp.S
		Gb := o.P.Shp.G

		// auxiliary
		divus = o.P.States[idx].Divus
		divvs = dc.α4*divus - o.U.divχs[idx] // divergence of Eq (35a) [1]

		// tpm variables
		plt = dc.β1*o.P.pl - o.P.ψl[idx] // Eq (35c) [1]
		klr = o.P.Mdl.Cnd.Klr(o.P.States[idx].Sl)
		ρL = o.P.States[idx].RhoL
		Cl = o.P.Mdl.Cl
		ρl, ρ, Cpl, Cvs, dρdpl, dpdpl, dCpldpl, dCvsdpl, dklrdpl, dCpldusM, dρldusM, dρdusM, err = o.P.States[idx].LSderivs(o.P.Mdl)
		if LogErr(err, "calc of tpm derivatives failed") {
			return
		}

		// debug
		if false {
			io.Pf("coef=%23.10f Gb=%v Klsat=%v dklrdpl = %23.10f\n", coef, Gb, o.P.Mdl.Klsat, dklrdpl)
			io.Pf("Sb=%v hl=%v klr=%v Cl=%v bs=%v\n", Sb, o.hl, klr, Cl, o.bs)
			io.Pf("dCpldpl=%v plt=%v dCvsdpl=%v divvs=%v β1=%v Cpl=%v\n", dCpldpl, plt, dCvsdpl, divvs, dc.β1, Cpl)
		}

		// Kpu, Kup and Kpp
		for n := 0; n < p_nverts; n++ {
			for j := 0; j < ndim; j++ {

				// Kpu := ∂Rl^n/∂us^m and Kup := ∂Rus^m/∂pl^n; see Eq (47) of [1]
				for m := 0; m < u_nverts; m++ {
					c = j + m*ndim

					// add ∂rlb/∂us^m: Eqs (A.3) and (A.6) of [1]
					o.Kpu[n][c] += coef * Sb[n] * (dCpldusM*plt + dc.α4*Cvs) * G[m][j]

					// add ∂(ρl.wl)/∂us^m: Eq (A.8) of [1]
					for i := 0; i < ndim; i++ {
						o.Kpu[n][c] += coef * Gb[n][i] * S[m] * dc.α1 * ρL * klr * o.P.Mdl.Klsat[i][j]
					}

					// add ∂rl/∂pl^n and ∂p/∂pl^n: Eqs (A.9) and (A.11) of [1]
					o.Kup[c][n] += coef * (S[m]*Sb[n]*dρdpl*o.bs[j] - G[m][j]*Sb[n]*dpdpl)

					// for seepage face
					if o.P.DoExtrap {
						o.dρldus_ex[n][c] += o.P.Emat[n][idx] * dρldusM * G[m][j]
					}
				}

				// term in brackets in Eq (A.7) of [1]
				o.P.tmp[j] = Sb[n]*dklrdpl*o.hl[j] - klr*(Sb[n]*Cl*o.bs[j]+Gb[n][j])
			}

			// Kpp := ∂Rl^m/∂pl^n; see Eq (47) of [1]
			for m := 0; m < p_nverts; m++ {

				// add ∂rlb/dpl^n: Eq (A.5) of [1]
				o.P.Kpp[m][n] += coef * Sb[m] * Sb[n] * (dCpldpl*plt + dCvsdpl*divvs + dc.β1*Cpl)

				// add ∂(ρl.wl)/∂us^m: Eq (A.7) of [1]
				for i := 0; i < ndim; i++ {
					for j := 0; j < ndim; j++ {
						o.P.Kpp[m][n] -= coef * Gb[m][i] * o.P.Mdl.Klsat[i][j] * o.P.tmp[j]
					}
				}

				// inner summation term in Eq (22) of [2]
				if o.P.DoExtrap {
					o.P.dρldpl_ex[m][n] += o.P.Emat[m][idx] * Cpl * Sb[n]
				}
			}

			// Eq. (19) of [2]
			if o.P.DoExtrap {
				o.P.ρl_ex[n] += o.P.Emat[n][idx] * ρl
			}
		}

		// Kuu: add ∂rub^m/∂us^n; see Eqs (47) and (A.10) of [1]
		for m := 0; m < u_nverts; m++ {
			for i := 0; i < ndim; i++ {
				r = i + m*ndim
				for n := 0; n < u_nverts; n++ {
					for j := 0; j < ndim; j++ {
						c = j + n*ndim
						o.U.K[r][c] += coef * S[m] * (S[n]*dc.α1*ρ*tsr.It[i][j] + dρdusM*o.bs[i]*G[n][j])
					}
				}
			}
		}

		// Kuu: add stiffness term ∂(σe・G^m)/∂us^n
		if LogErr(o.U.MdlSmall.CalcD(o.U.D, o.U.States[idx], firstIt), "AddToKb") {
			return
		}
		IpAddToKt(o.U.K, u_nverts, ndim, coef, G, o.U.D)
	}

	// contribution from natural boundary conditions
	if o.P.HasSeep {
		if !o.P.add_natbcs_to_jac(sol) {
			return
		}
		if !o.add_natbcs_to_jac(sol) {
			return
		}
	}

	// debug
	//if true {
	if false {
		if o.Id() == 69 {
			o.debug_print_K()
			panic("stop")
		}
	}

	// add K to sparse matrix Kb
	//    _             _
	//   |  Kuu Kup  0   |
	//   |  Kpu Kpp Kpf  |
	//   |_ Kfu Kfp Kff _|
	//
	for i, I := range o.P.Pmap {
		for j, J := range o.P.Pmap {
			Kb.Put(I, J, o.P.Kpp[i][j])
		}
		for j, J := range o.P.Fmap {
			Kb.Put(I, J, o.P.Kpf[i][j])
			Kb.Put(J, I, o.P.Kfp[j][i])
		}
		for j, J := range o.U.Umap {
			Kb.Put(I, J, o.Kpu[i][j])
			Kb.Put(J, I, o.Kup[j][i])
		}
	}
	for i, I := range o.P.Fmap {
		for j, J := range o.P.Fmap {
			Kb.Put(I, J, o.P.Kff[i][j])
		}
	}
	for i, I := range o.U.Umap {
		for j, J := range o.U.Umap {
			Kb.Put(I, J, o.U.K[i][j])
		}
		// TODO:
		//for j, J := range o.P.Fmap {
		//Kb.Put(J, I, o.Kfu[j][i])
		//}
	}
	return true
}

// Update perform (tangent) update
func (o *ElemUP) Update(sol *Solution) (ok bool) {

	// auxiliary
	ndim := o.Ndim

	// for each integration point
	var Δpl, divusNew float64
	var r int
	for idx, ip := range o.U.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.P.Shp.CalcAtIp(o.P.X, ip, false), "Update") {
			return
		}
		if LogErr(o.U.Shp.CalcAtIp(o.U.X, ip, true), "Update") {
			return
		}

		// compute Δpl @ ip
		Δpl = 0
		for m := 0; m < o.P.Shp.Nverts; m++ {
			r = o.P.Pmap[m]
			Δpl += o.P.Shp.S[m] * sol.ΔY[r]
		}

		// compute divus @ ip
		divusNew = 0
		for m := 0; m < o.U.Shp.Nverts; m++ {
			for i := 0; i < ndim; i++ {
				r := o.U.Umap[i+m*ndim]
				divusNew += o.U.Shp.G[m][i] * sol.Y[r]
			}
		}

		// p: update internal state
		if LogErr(o.P.Mdl.Update(o.P.States[idx], Δpl, 0, divusNew), "p: update failed") {
			return
		}

		// u: update internal state
		if !o.U.ipupdate(idx, o.U.Shp.S, o.U.Shp.G, sol) {
			return
		}

		// debug
		//io.Pf("%3d : Δpl=%13.10f pcNew=%13.10f sl=%13.10f RhoL=%13.10f Wet=%v\n", o.Id(), Δpl, o.P.States[idx].Pg-o.P.States[idx].Pl, o.P.States[idx].Sl, o.P.States[idx].RhoL, o.P.States[idx].Wet)
	}
	return true
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
func (o *ElemUP) ipvars(idx int, sol *Solution) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.P.Shp.CalcAtIp(o.P.X, o.U.IpsElem[idx], true), "ipvars") {
		return
	}
	if LogErr(o.U.Shp.CalcAtIp(o.U.X, o.U.IpsElem[idx], true), "ipvars") {
		return
	}

	// auxiliary
	ndim := o.Ndim
	dc := Global.DynCoefs
	ρL := o.P.States[idx].RhoL

	// gravity
	o.P.g[ndim-1] = 0
	if o.P.Gfcn != nil {
		o.P.g[ndim-1] = -o.P.Gfcn.F(sol.T, nil)
	}

	// clear gpl and recover u-variables @ ip
	for i := 0; i < ndim; i++ {
		o.P.gpl[i] = 0 // clear gpl here
		o.U.us[i] = 0
		for m := 0; m < o.U.Shp.Nverts; m++ {
			r := o.U.Umap[i+m*ndim]
			o.U.us[i] += o.U.Shp.S[m] * sol.Y[r]
		}
	}

	// recover p-variables @ ip
	o.P.pl = 0
	for m := 0; m < o.P.Shp.Nverts; m++ {
		r := o.P.Pmap[m]
		o.P.pl += o.P.Shp.S[m] * sol.Y[r]
		for i := 0; i < ndim; i++ {
			o.P.gpl[i] += o.P.Shp.G[m][i] * sol.Y[r]
		}
	}

	// compute bs and hl. see Eqs (A.1) of [1]
	for i := 0; i < ndim; i++ {
		o.bs[i] = dc.α1*o.U.us[i] - o.U.ζs[idx][i] - o.P.g[i]
		o.hl[i] = -ρL*o.bs[i] - o.P.gpl[i]
	}
	return true
}

// add_natbcs_to_jac adds contribution from natural boundary conditions to Jacobian
func (o ElemUP) add_natbcs_to_jac(sol *Solution) (ok bool) {

	// compute surface integral
	ndim := o.Ndim
	u_nverts := o.U.Shp.Nverts
	var tmp float64
	var z, pl, fl, plmax, g, rmp float64
	for _, nbc := range o.P.NatBcs {

		// temporary value == function evaluation
		tmp = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for _, ipf := range o.P.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			if LogErr(o.P.Shp.CalcAtFaceIp(o.P.X, ipf, iface), "up: add_natbcs_to_jac") {
				return
			}
			Sf := o.P.Shp.Sf
			Jf := la.VecNorm(o.P.Shp.Fnvec)
			coef := ipf.W * Jf

			// select natural boundary condition type
			switch nbc.Key {
			case "seepP", "seepH":

				// variables extrapolated to face
				_, z, pl, fl = o.P.fipvars(iface, sol)

				// specified condition
				plmax = tmp
				if nbc.Key == "seepH" {
					plmax = max(o.P.γl*(tmp-z), 0.0)
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.P.ramp(fl + o.P.κ*g)
				for i, m := range o.P.Shp.FaceLocalV[iface] {
					for n := 0; n < u_nverts; n++ {
						for j := 0; j < ndim; j++ {
							c := j + n*ndim
							for l, r := range o.P.Shp.FaceLocalV[iface] {
								o.Kpu[m][c] += coef * Sf[i] * Sf[l] * o.dρldus_ex[r][c] * rmp
							}
						}
					}
				}
			}
		}
	}
	return true
}

func (o ElemUP) debug_print_K() {
	la.PrintMat("Kpp", o.P.Kpp, "%20.10f", false)
	la.PrintMat("Kpf", o.P.Kpf, "%20.10f", false)
	la.PrintMat("Kfp", o.P.Kfp, "%20.10f", false)
	la.PrintMat("Kff", o.P.Kff, "%20.10f", false)
	la.PrintMat("Kpu", o.Kpu, "%20.10f", false)
	la.PrintMat("Kup", o.Kup, "%20.10f", false)
	la.PrintMat("Kuu", o.U.K, "%20.10f", false)
}
