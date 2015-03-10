// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
)

// ElemU represents a solid element with displacements u as primary variables
type ElemU struct {

	// basic data
	Cid int         // cell/element id
	X   [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Shp *shp.Shape  // shape structure
	Nu  int         // total number of unknowns

	// variables for dynamics
	Rho  float64  // density of solids
	Cdam float64  // coefficient for damping
	Gfcn fun.Func // gravity function

	// optional data
	UseB      bool    // use B matrix
	Thickness float64 // thickness (for plane-stress)
	Debug     bool    // debugging flag

	// integration points
	IpsElem []*shp.Ipoint // integration points of element
	IpsFace []*shp.Ipoint // integration points corresponding to faces

	// material model and internal variables
	Model    msolid.Model // material model
	MdlSmall msolid.Small // model specialisation for small strains
	MdlLarge msolid.Large // model specialisation for large deformations

	// internal variables
	States    []*msolid.State // [nip] states
	StatesBkp []*msolid.State // [nip] backup states

	// problem variables
	Umap []int // assembly map (location array/element equations)

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
	infogetters["u"] = func(cellType string, faceConds []*FaceCond) *Info {

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
	eallocators["u"] = func(cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemU
		o.Cid = cid
		o.X = x
		o.Shp = shp.Get(cellType)
		ndim := Global.Ndim
		o.Nu = ndim * o.Shp.Nverts

		// parse flags
		o.UseB, o.Debug, o.Thickness = GetSolidFlags(edat.Extra)

		// integration points
		o.IpsElem, o.IpsFace = GetIntegrationPoints(edat.Nip, edat.Nipf, cellType)
		if o.IpsElem == nil || o.IpsFace == nil {
			return nil
		}
		nip := len(o.IpsElem)

		// model
		var prms fun.Prms
		o.Model, prms = GetAndInitSolidModel(edat.Mat, ndim)
		if o.Model == nil {
			return nil
		}

		// model specialisations
		switch m := o.Model.(type) {
		case msolid.Small:
			o.MdlSmall = m
		case msolid.Large:
			o.MdlLarge = m
		}

		// parameters
		for _, p := range prms {
			switch p.N {
			case "rho":
				o.Rho = p.V
			case "Cdam":
				o.Cdam = p.V
			}
		}

		// local starred variables
		o.ζs = la.MatAlloc(nip, ndim)
		o.χs = la.MatAlloc(nip, ndim)
		o.divχs = make([]float64, nip)

		// scratchpad. computed @ each ip
		nsig := 2 * ndim
		o.grav = make([]float64, ndim)
		o.us = make([]float64, ndim)
		o.fi = make([]float64, o.Nu)
		o.D = la.MatAlloc(nsig, nsig)
		o.K = la.MatAlloc(o.Nu, o.Nu)
		if o.UseB {
			o.B = la.MatAlloc(nsig, o.Nu)
		}

		// strains
		o.ε = make([]float64, nsig)
		o.Δε = make([]float64, nsig)

		// variables for debugging
		if o.Debug {
			o.fex = make([]float64, o.Shp.Nverts)
			o.fey = make([]float64, o.Shp.Nverts)
			if ndim == 3 {
				o.fez = make([]float64, o.Shp.Nverts)
			}
		}

		// surface loads (natural boundary conditions)
		for _, fc := range faceConds {
			o.NatBcs = append(o.NatBcs, &NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o ElemU) Id() int { return o.Cid }

// SetEqs set equations
func (o *ElemU) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
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

// SetEleConds set element conditions
func (o *ElemU) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	if key == "g" { // gravity
		o.Gfcn = f
	}
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemU) InterpStarVars(sol *Solution) (ok bool) {

	// skip steady cases
	if Global.Sim.Data.Steady {
		return true
	}

	// for each integration point
	ndim := Global.Ndim
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		// interpolate starred variables
		o.divχs[idx] = 0
		for i := 0; i < ndim; i++ {
			o.ζs[idx][i] = 0
			o.χs[idx][i] = 0
			for m := 0; m < o.Shp.Nverts; m++ {
				r := o.Umap[i+m*ndim]
				o.ζs[idx][i] += o.Shp.S[m] * sol.Zet[r]
				o.χs[idx][i] += o.Shp.S[m] * sol.Chi[r]
				o.divχs[idx] += o.Shp.G[m][i] * sol.Chi[r]
			}
		}
	}
	return
}

// AddToRhs adds -R to global residual vector fb
func (o *ElemU) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// clear fi vector if using B matrix
	if o.UseB {
		la.VecFill(o.fi, 0)
	}

	// for each integration point
	dc := Global.DynCoefs
	ndim := Global.Ndim
	nverts := o.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// auxiliary
		coef := o.Shp.J * ip.W * o.Thickness
		S := o.Shp.S
		G := o.Shp.G

		// add internal forces to fb
		if o.UseB {
			radius := 1.0
			if Global.Sim.Data.Axisym {
				radius = o.Shp.AxisymGetRadius(o.X)
				coef *= radius
			}
			IpBmatrix(o.B, ndim, nverts, G, Global.Sim.Data.Axisym, radius, S)
			la.MatTrVecMulAdd(o.fi, coef, o.B, o.States[idx].Sig) // fi += coef * tr(B) * σ
		} else {
			for m := 0; m < nverts; m++ {
				for i := 0; i < ndim; i++ {
					r := o.Umap[i+m*ndim]
					for j := 0; j < ndim; j++ {
						fb[r] -= coef * tsr.M2T(o.States[idx].Sig, i, j) * G[m][j] // -fi
					}
				}
			}
		}

		// dynamic term
		if !Global.Sim.Data.Steady {
			for m := 0; m < nverts; m++ {
				for i := 0; i < ndim; i++ {
					r := o.Umap[i+m*ndim]
					fb[r] -= coef * S[m] * (o.Rho*(dc.α1*o.us[i]-o.ζs[idx][i]-o.grav[i]) + o.Cdam*(dc.α4*o.us[i]-o.χs[idx][i])) // -RuBar
				}
			}
		}
	}

	// assemble fb if using B matrix
	if o.UseB {
		for i, I := range o.Umap {
			fb[I] -= o.fi[i]
		}
	}

	// external forces
	return o.add_surfloads_to_rhs(fb, sol)
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o *ElemU) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// zero K matrix
	la.MatFill(o.K, 0)

	// for each integration point
	dc := Global.DynCoefs
	ndim := Global.Ndim
	nverts := o.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// check Jacobian
		if o.Shp.J < 0 {
			LogErrCond(true, "ElemU: eid=%d: Jacobian is negative = %g\n", o.Id(), o.Shp.J)
			return
		}

		// auxiliary
		coef := o.Shp.J * ip.W * o.Thickness
		S := o.Shp.S
		G := o.Shp.G

		// consistent tangent model matrix
		if LogErr(o.MdlSmall.CalcD(o.D, o.States[idx], firstIt), "AddToKb") {
			return
		}

		// add contribution to consistent tangent matrix
		if o.UseB {
			radius := 1.0
			if Global.Sim.Data.Axisym {
				radius = o.Shp.AxisymGetRadius(o.X)
				coef *= radius
			}
			IpBmatrix(o.B, ndim, nverts, G, Global.Sim.Data.Axisym, radius, S)
			la.MatTrMulAdd3(o.K, coef, o.B, o.D, o.B) // K += coef * tr(B) * D * B
		} else {
			IpAddToKt(o.K, nverts, ndim, coef, G, o.D)
		}

		// dynamic term
		if !Global.Sim.Data.Steady {
			for m := 0; m < nverts; m++ {
				for i := 0; i < ndim; i++ {
					r := i + m*ndim
					for n := 0; n < nverts; n++ {
						c := i + n*ndim
						o.K[r][c] += coef * S[m] * S[n] * (o.Rho*dc.α1 + o.Cdam*dc.α4)
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
func (o *ElemU) Update(sol *Solution) (ok bool) {

	// for each integration point
	for idx, _ := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// update internal state
		if !o.ipupdate(idx, o.Shp.S, o.Shp.G, sol) {
			return
		}
	}
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// Ipoints returns the real coordinates of integration points [nip][ndim]
func (o ElemU) Ipoints() (coords [][]float64) {
	coords = la.MatAlloc(len(o.IpsElem), Global.Ndim)
	for idx, ip := range o.IpsElem {
		coords[idx] = o.Shp.IpRealCoords(o.X, ip)
	}
	return
}

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *ElemU) SetIniIvs(sol *Solution, ivs map[string][]float64) (ok bool) {

	// allocate slices of states
	nip := len(o.IpsElem)
	o.States = make([]*msolid.State, nip)
	o.StatesBkp = make([]*msolid.State, nip)

	// for each integration point
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Model.InitIntVars()
		o.StatesBkp[i] = o.States[i].GetCopy()
	}

	// initial stresses
	if _, ok := ivs["sx"]; ok {
		for i := 0; i < nip; i++ {
			Ivs2sigmas(o.States[i].Sig, i, ivs)
			copy(o.StatesBkp[i].Sig, o.States[i].Sig)
		}
	}
	return true
}

// BackupIvs create copy of internal variables
func (o *ElemU) BackupIvs() (ok bool) {
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return true
}

// RestoreIvs restore internal variables from copies
func (o *ElemU) RestoreIvs() (ok bool) {
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return true
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *ElemU) Ureset(sol *Solution) (ok bool) {
	for idx, _ := range o.IpsElem {
		if len(o.States[idx].F) > 0 {
			la.MatFill(o.States[idx].F, 0)
			la.MatFill(o.StatesBkp[idx].F, 0)
		}
	}
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElemU) Encode(enc Encoder) (ok bool) {
	return !LogErr(enc.Encode(o.States), "Encode")
}

// Decode decodes internal variables
func (o ElemU) Decode(dec Decoder) (ok bool) {
	if LogErr(dec.Decode(&o.States), "Decode") {
		return
	}
	return o.BackupIvs()
}

// OutIpsData returns data from all integration points for output
func (o ElemU) OutIpsData() (data []*OutIpData) {
	sigmas := StressKeys()
	for idx, ip := range o.IpsElem {
		s := o.States[idx]
		x := o.Shp.IpRealCoords(o.X, ip)
		v := make(map[string]*float64)
		for i, key := range sigmas {
			v[key] = &s.Sig[i]
		}
		data = append(data, &OutIpData{o.Id(), x, v})
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipudpate updates internal state
func (o *ElemU) ipupdate(idx int, S []float64, G [][]float64, sol *Solution) (ok bool) {

	// compute strains
	ndim := Global.Ndim
	nverts := o.Shp.Nverts
	if o.UseB {
		radius := 1.0
		if Global.Sim.Data.Axisym {
			radius = o.Shp.AxisymGetRadius(o.X)
		}
		IpBmatrix(o.B, ndim, nverts, G, Global.Sim.Data.Axisym, radius, S)
		IpStrainsAndIncB(o.ε, o.Δε, 2*ndim, o.Nu, o.B, sol.Y, sol.ΔY, o.Umap)
	} else {
		IpStrainsAndInc(o.ε, o.Δε, nverts, ndim, sol.Y, sol.ΔY, o.Umap, G)
	}

	// call model update => update stresses
	if LogErr(o.MdlSmall.Update(o.States[idx], o.ε, o.Δε), "ipupdate") {
		return
	}
	return true
}

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemU) ipvars(idx int, sol *Solution) (ok bool) {

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

// surfloads_keys returns the keys that can be used to specify surface loads
func (o *ElemU) surfloads_keys() map[string]bool {
	return map[string]bool{"qn": true, "qn0": true, "aqn": true}
}

// add_surfloads_to_rhs adds surfaces loads to rhs
func (o *ElemU) add_surfloads_to_rhs(fb []float64, sol *Solution) (ok bool) {

	// debugging variables
	ndim := Global.Ndim
	if o.Debug {
		la.VecFill(o.fex, 0)
		la.VecFill(o.fey, 0)
		if ndim == 3 {
			la.VecFill(o.fez, 0)
		}
	}

	// compute surface integral
	for _, load := range o.NatBcs {
		for _, ip := range o.IpsFace {
			if LogErr(o.Shp.CalcAtFaceIp(o.X, ip, load.IdxFace), "add_surfloads_to_rhs") {
				return
			}
			switch load.Key {
			case "qn", "qn0", "aqn":
				coef := ip.W * load.Fcn.F(sol.T, nil) * o.Thickness
				nvec := o.Shp.Fnvec
				Sf := o.Shp.Sf
				if Global.Sim.Data.Axisym && load.Key == "aqn" {
					coef *= o.Shp.AxisymGetRadiusF(o.X, load.IdxFace)
				}
				for j, m := range o.Shp.FaceLocalV[load.IdxFace] {
					for i := 0; i < ndim; i++ {
						r := o.Umap[i+m*ndim]
						fb[r] += coef * Sf[j] * nvec[i] // +fe
					}
					if o.Debug {
						o.fex[m] += coef * Sf[j] * nvec[0]
						o.fey[m] += coef * Sf[j] * nvec[1]
						if ndim == 3 {
							o.fez[m] += coef * Sf[j] * nvec[2]
						}
					}
				}
			}
		}
	}
	return true
}
