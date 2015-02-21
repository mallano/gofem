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
	"github.com/cpmech/gosl/tsr"
)

// ElemU represents a solid element with displacements u as primary variables
type ElemU struct {

	// basic data
	Cell *inp.Cell   // cell
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Ndim int         // space dimension
	Nu   int         // total number of unknowns

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
	Model    msolid.Solid // material model
	MdlSmall msolid.Small // model specialisation for small strains
	MdlLarge msolid.Large // model specialisation for large deformations

	// internal variables
	Sig0      [][]float64     // [nip][nsig] initial stresses
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
	iallocators["u"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info {

		// new info
		var info Info

		// number of nodes in element
		cell := msh.Cells[cid]
		nverts := cell.Shp.Nverts

		// solution variables
		ykeys := []string{"ux", "uy"}
		if msh.Ndim == 3 {
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
	eallocators["u"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {

		// basic data
		var o ElemU
		o.Cell = msh.Cells[cid]
		o.X = BuildCoordsMatrix(o.Cell, msh)
		o.Ndim = msh.Ndim
		o.Nu = o.Ndim * o.Cell.Shp.Nverts

		// flag: use B matrix
		if s_useB, found := io.Keycode(edat.Extra, "useB"); found {
			o.UseB = io.Atob(s_useB)
		}

		if Global.Sim.Data.Axisym {
			o.UseB = true
		}

		// flag: thickess => plane-stress
		o.Thickness = 1.0
		if s_thick, found := io.Keycode(edat.Extra, "thick"); found {
			o.Thickness = io.Atof(s_thick)
			if LogErrCond(!Global.Sim.Data.Pstress, "cannot specify 'thick' in element if Pstress=false in Global Data") {
				return nil
			}
		}

		// flag: debug
		if s_debug, found := io.Keycode(edat.Extra, "debug"); found {
			o.Debug = io.Atob(s_debug)
		}

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

		// material model name
		matname := edat.Mat
		matdata := Global.Mdb.Get(matname)
		if LogErrCond(matdata == nil, "materials database failed on getting %q material\n", matname) {
			return nil
		}
		mdlname := matdata.Model

		// model and its specialisations
		o.Model = msolid.GetModel(Global.Sim.Data.FnameKey, matname, mdlname, false)
		if LogErrCond(o.Model == nil, "cannot find model named %s\n", mdlname) {
			return nil
		}
		err = o.Model.Init(o.Ndim, Global.Sim.Data.Pstress, matdata.Prms)
		if LogErr(err, "Model.Init failed") {
			return nil
		}
		switch m := o.Model.(type) {
		case msolid.Small:
			o.MdlSmall = m
		case msolid.Large:
			o.MdlLarge = m
		}

		// parameters
		for _, p := range matdata.Prms {
			switch p.N {
			case "rho":
				o.Rho = p.V
			case "Cdam":
				o.Cdam = p.V
			}
		}

		// local starred variables
		o.ζs = la.MatAlloc(nip, o.Ndim)
		o.χs = la.MatAlloc(nip, o.Ndim)
		o.divχs = make([]float64, nip)

		// scratchpad. computed @ each ip
		nsig := 2 * o.Ndim
		o.grav = make([]float64, o.Ndim)
		o.us = make([]float64, o.Ndim)
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
			o.fex = make([]float64, o.Cell.Shp.Nverts)
			o.fey = make([]float64, o.Cell.Shp.Nverts)
			if o.Ndim == 3 {
				o.fez = make([]float64, o.Cell.Shp.Nverts)
			}
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *ElemU) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	o.Umap = make([]int, o.Nu)
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
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

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *ElemU) SetNatBcs(key string, idxface int, f fun.Func, extra string) (ok bool) {
	o.NatBcs = append(o.NatBcs, &NaturalBc{key, idxface, f, extra})
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *ElemU) InterpStarVars(sol *Solution) (ok bool) {

	// skip steady cases
	if Global.Sim.Data.Steady {
		return true
	}

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Cell.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		// interpolate starred variables
		o.divχs[idx] = 0
		for i := 0; i < o.Ndim; i++ {
			o.ζs[idx][i] = 0
			o.χs[idx][i] = 0
			for m := 0; m < o.Cell.Shp.Nverts; m++ {
				r := o.Umap[i+m*o.Ndim]
				o.ζs[idx][i] += o.Cell.Shp.S[m] * sol.Zet[r]
				o.χs[idx][i] += o.Cell.Shp.S[m] * sol.Chi[r]
				o.divχs[idx] += o.Cell.Shp.G[m][i] * sol.Chi[r]
			}
		}
	}
	return
}

// adds -R to global residual vector fb
func (o *ElemU) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// clear fi vector if using B matrix
	if o.UseB {
		la.VecFill(o.fi, 0)
	}

	// for each integration point
	dc := Global.DynCoefs
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef := o.Cell.Shp.J * ip.W * o.Thickness
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// add internal forces to fb
		if o.UseB {
			radius := 1.0
			if Global.Sim.Data.Axisym {
				radius = o.Cell.Shp.AxisymGetRadius(o.X)
				coef *= radius
			}
			IpBmatrix(o.B, o.Ndim, nverts, G, Global.Sim.Data.Axisym, radius, S)
			la.MatTrVecMulAdd(o.fi, coef, o.B, o.States[idx].Sig) // fi += coef * tr(B) * σ
		} else {
			for m := 0; m < nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					r := o.Umap[i+m*o.Ndim]
					for j := 0; j < o.Ndim; j++ {
						fb[r] -= coef * tsr.M2T(o.States[idx].Sig, i, j) * G[m][j] // -fi
					}
				}
			}
		}

		// dynamic term
		if !Global.Sim.Data.Steady {
			for m := 0; m < nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					r := o.Umap[i+m*o.Ndim]
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

// adds element K to global Jacobian matrix Kb
func (o *ElemU) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// zero K matrix
	la.MatFill(o.K, 0)

	// for each integration point
	dc := Global.DynCoefs
	nverts := o.Cell.Shp.Nverts
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}

		// check Jacobian
		if o.Cell.Shp.J < 0 {
			LogErrCond(true, "ElemU: eid=%d: Jacobian is negative = %g\n", o.Cell.Id, o.Cell.Shp.J)
			return
		}

		// coefficients
		coef := o.Cell.Shp.J * ip.W * o.Thickness
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// consistent tangent model matrix
		if LogErr(o.MdlSmall.CalcD(o.D, o.States[idx], firstIt), "AddToKb") {
			return
		}

		// add contribution to consistent tangent matrix
		if o.UseB {
			radius := 1.0
			if Global.Sim.Data.Axisym {
				radius = o.Cell.Shp.AxisymGetRadius(o.X)
				coef *= radius
			}
			IpBmatrix(o.B, o.Ndim, nverts, G, Global.Sim.Data.Axisym, radius, S)
			la.MatTrMulAdd3(o.K, coef, o.B, o.D, o.B) // K += coef * tr(B) * D * B
		} else {
			IpAddToKt(o.K, nverts, o.Ndim, coef, G, o.D)
		}

		// dynamic term
		if !Global.Sim.Data.Steady {
			for m := 0; m < nverts; m++ {
				for i := 0; i < o.Ndim; i++ {
					r := i + m*o.Ndim
					for n := 0; n < nverts; n++ {
						c := i + n*o.Ndim
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
	var σ0 []float64
	nverts := o.Cell.Shp.Nverts
	for idx, _ := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		S := o.Cell.Shp.S
		G := o.Cell.Shp.G

		// compute strains
		if o.UseB {
			radius := 1.0
			if Global.Sim.Data.Axisym {
				radius = o.Cell.Shp.AxisymGetRadius(o.X)
			}
			IpBmatrix(o.B, o.Ndim, nverts, G, Global.Sim.Data.Axisym, radius, S)
			IpStrainsAndIncB(o.ε, o.Δε, 2*o.Ndim, o.Nu, o.B, sol.Y, sol.ΔY, o.Umap)
		} else {
			IpStrainsAndInc(o.ε, o.Δε, nverts, o.Ndim, sol.Y, sol.ΔY, o.Umap, G)
		}

		// call model update => update stresses
		if len(o.Sig0) > 0 {
			σ0 = o.Sig0[idx]
		}
		if LogErr(o.MdlSmall.Update(o.States[idx], σ0, o.ε, o.Δε), "Update") {
			return
		}
	}
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *ElemU) InitIvs(sol *Solution) (ok bool) {
	nip := len(o.IpsElem)
	o.States = make([]*msolid.State, nip)
	o.StatesBkp = make([]*msolid.State, nip)
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Model.InitIntVars()
		o.StatesBkp[i] = o.States[i].GetCopy()
	}
	return true
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *ElemU) SetIvs(ivs map[string][]float64) (ok bool) {
	nip := len(o.IpsElem)
	o.Sig0 = Ivs2sigmas(nip, o.Ndim, ivs)
	for i := 0; i < nip; i++ {
		copy(o.States[i].Sig, o.Sig0[i])
		copy(o.StatesBkp[i].Sig, o.Sig0[i])
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

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o ElemU) Encode(enc Encoder) (ok bool) {
	return !LogErr(enc.Encode(o.States), "Encode")
}

// Decode decodes internal variables
func (o ElemU) Decode(dec Decoder) (ok bool) {
	return !LogErr(dec.Decode(&o.States), "Decode")
}

// OutIpsData returns data from all integration points for output
func (o ElemU) OutIpsData() (labels []string, data []*OutIpData) {
	labels = []string{"sx", "sy", "sz", "sxy"}
	if o.Ndim == 3 {
		labels = append(labels, "syz", "szx")
	}
	for idx, ip := range o.IpsElem {
		s := o.States[idx]
		vals := []*float64{&s.Sig[0], &s.Sig[1], &s.Sig[2], &s.Sig[3]}
		if o.Ndim == 3 {
			vals = append(vals, &s.Sig[4], &s.Sig[5])
		}
		x := o.Cell.Shp.IpRealCoords(o.X, ip)
		data = append(data, &OutIpData{vals, x})
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemU) ipvars(idx int, sol *Solution) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.Cell.Shp.CalcAtIp(o.X, o.IpsElem[idx], true), "ipvars") {
		return
	}

	// skip if steady (this must be after CalcAtIp, because callers will need S and G)
	if Global.Sim.Data.Steady {
		return true
	}

	// clear variables
	for i := 0; i < o.Ndim; i++ {
		o.us[i] = 0
	}

	// recover u-variables @ ip
	for m := 0; m < o.Cell.Shp.Nverts; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := o.Umap[i+m*o.Ndim]
			o.us[i] += o.Cell.Shp.S[m] * sol.Y[r]
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
	if o.Debug {
		la.VecFill(o.fex, 0)
		la.VecFill(o.fey, 0)
		if o.Ndim == 3 {
			la.VecFill(o.fez, 0)
		}
	}

	// compute surface integral
	for _, load := range o.NatBcs {
		for _, ip := range o.IpsFace {
			if LogErr(o.Cell.Shp.CalcAtFaceIp(o.X, ip, load.IdxFace), "add_surfloads_to_rhs") {
				return
			}
			switch load.Key {
			case "qn", "qn0", "aqn":
				coef := ip.W * load.Fcn.F(sol.T, nil) * o.Thickness
				nvec := o.Cell.Shp.Fnvec
				Sf := o.Cell.Shp.Sf
				if Global.Sim.Data.Axisym && load.Key == "aqn" {
					coef *= o.Cell.Shp.AxisymGetRadiusF(o.X, load.IdxFace)
				}
				for j, m := range o.Cell.Shp.FaceLocalV[load.IdxFace] {
					for i := 0; i < o.Ndim; i++ {
						r := o.Umap[i+m*o.Ndim]
						fb[r] += coef * Sf[j] * nvec[i] // +fe
					}
					if o.Debug {
						o.fex[m] += coef * Sf[j] * nvec[0]
						o.fey[m] += coef * Sf[j] * nvec[1]
						if o.Ndim == 3 {
							o.fez[m] += coef * Sf[j] * nvec[2]
						}
					}
				}
			}
		}
	}
	return true
}
