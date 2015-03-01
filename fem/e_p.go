// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/fun"
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
	Cid  int         // cell/element id
	X    [][]float64 // matrix of nodal coordinates [ndim][nnode]
	Ndim int         // space dimension
	Shp  *shp.Shape  // shape structure
	Np   int         // total number of unknowns == number of vertices

	// integration points
	IpsElem []*shp.Ipoint // integration points of element
	IpsFace []*shp.Ipoint // integration points corresponding to faces

	// material model
	Mdl *mporous.Model // model
	γl  float64        // unit weight of liquid

	// problem variables
	Pmap []int // assembly map (location array/element equations)

	// internal variables
	States    []*mporous.State
	StatesBkp []*mporous.State

	// gravity
	Gfcn fun.Func // gravity function

	// natural boundary conditions
	NatBcs []*NaturalBc // natural boundary conditions

	// flux boundary conditions (qb == \bar{q})
	ρl_ex     []float64   // [nverts] ρl extrapolted to nodes => if has qb (flux)
	dρldpl_ex [][]float64 // [nverts][nverts] Cpl extrapolted to nodes => if has qb (flux)
	Emat      [][]float64 // [nverts][nips] extrapolator matrix
	DoExtrap  bool        // do extrapolation of ρl and Cpl => for use with flux and seepage conditions

	// seepage face (seepP or seepH)
	Nf         int     // number of fl variables
	HasSeep    bool    // indicates if this element has seepage faces
	Vid2seepId []int   // [nverts] maps local vertex id to index in Fmap
	SeepId2vid []int   // [nf] maps seepage face variable id to local vertex id
	Fmap       []int   // [nf] map of "fl" variables (seepage face)
	Macaulay   bool    // use discrete ramp function instead of smooth ramp
	βrmp       float64 // coefficient for Sramp
	κ          float64 // κ coefficient to normalise equation for seepage face modelling

	// local starred variables
	ψl []float64 // [nip] ψl* = β1.p + β2.dpdt

	// scratchpad. computed @ each ip
	g   []float64   // [ndim] gravity vector
	pl  float64     // pl: liquid pressure
	gpl []float64   // [ndim] ∇pl: gradient of liquid pressure
	ρwl []float64   // [ndim] ρl*wl: weighted liquid relative velocity
	tmp []float64   // [ndim] temporary (auxiliary) vector
	Kpp [][]float64 // [np][np] Kpp := dRpl/dpl consistent tangent matrix
	Kpf [][]float64 // [np][nf] Kpf := dRpl/dfl consistent tangent matrix
	Kfp [][]float64 // [nf][np] Kfp := dRfl/dpl consistent tangent matrix
	Kff [][]float64 // [nf][nf] Kff := dRfl/dfl consistent tangent matrix
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["p"] = func(ndim int, cellType string, faceConds []*FaceCond) *Info {

		// new info
		var info Info

		// number of nodes in element
		nverts := shp.GetNverts(cellType)

		// solution variables
		ykeys := []string{"pl"}
		info.Dofs = make([][]string, nverts)
		for m := 0; m < nverts; m++ {
			info.Dofs[m] = ykeys
		}

		// vertices on seepage faces
		lverts := GetVertsWithCond(faceConds, "seepP", "seepH")
		for _, m := range lverts {
			if m < nverts { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
				info.Dofs[m] = append(info.Dofs[m], "fl")
			}
		}
		if len(lverts) > 0 {
			ykeys = append(ykeys, "fl")
		}

		// maps
		info.Y2F = map[string]string{"pl": "ql"}

		// t1 and t2 variables
		info.T1vars = ykeys
		return &info
	}

	// element allocator
	eallocators["p"] = func(ndim int, cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem {

		// basic data
		var o ElemP
		o.Cid = cid
		o.X = x
		o.Ndim = ndim
		o.Shp = shp.Get(cellType)
		o.Np = o.Shp.Nverts

		// integration points
		o.IpsElem, o.IpsFace = GetIntegrationPoints(edat.Extra, cellType)
		if o.IpsElem == nil || o.IpsFace == nil {
			return nil
		}
		nip := len(o.IpsElem)

		// models
		o.Mdl = GetAndInitPorousModel(edat.Mat)
		if o.Mdl == nil {
			return nil
		}

		// unit weight of liquid
		o.γl = o.Mdl.RhoL0 * o.Mdl.Gref

		// local starred variables
		o.ψl = make([]float64, nip)

		// scratchpad. computed @ each ip
		o.g = make([]float64, o.Ndim)
		o.gpl = make([]float64, o.Ndim)
		o.ρwl = make([]float64, o.Ndim)
		o.tmp = make([]float64, o.Ndim)
		o.Kpp = la.MatAlloc(o.Np, o.Np)

		// vertices on seepage faces
		var seepverts []int
		lverts := GetVertsWithCond(faceConds, "seepP", "seepH")
		for _, m := range lverts {
			if m < o.Np { // avoid adding vertices of superelement (e.g. qua8 vertices in this qua4 cell)
				seepverts = append(seepverts, m)
			}
		}

		o.Nf = len(seepverts)
		o.HasSeep = o.Nf > 0
		if o.HasSeep {

			// vertices on seepage face; numbering
			o.SeepId2vid = seepverts
			o.Vid2seepId = utl.IntVals(o.Np, -1)
			o.Fmap = make([]int, o.Nf)
			for μ, m := range o.SeepId2vid {
				o.Vid2seepId[m] = μ
			}

			// flags
			o.Macaulay, o.βrmp, o.κ = GetSeepFaceFlags(edat.Extra)

			// allocate coupling matrices
			o.Kpf = la.MatAlloc(o.Np, o.Nf)
			o.Kfp = la.MatAlloc(o.Nf, o.Np)
			o.Kff = la.MatAlloc(o.Nf, o.Nf)
		}

		// set natural boundary conditions
		for _, fc := range faceConds {
			o.NatBcs = append(o.NatBcs, &NaturalBc{fc.Cond, fc.FaceId, fc.Func, fc.Extra})
			if fc.Cond == "qb" || fc.Cond == "seepP" || fc.Cond == "seepH" {
				nv := o.Shp.Nverts
				nip := len(o.IpsElem)
				o.ρl_ex = make([]float64, nv)
				o.dρldpl_ex = la.MatAlloc(nv, nv)
				o.Emat = la.MatAlloc(nv, nip)
				o.DoExtrap = true
				if LogErr(o.Shp.Extrapolator(o.Emat, o.IpsElem), "element allocation") {
					return nil
				}
			}
		}

		// return new element
		return &o
	}
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// Id returns the cell Id
func (o ElemP) Id() int { return o.Cid }

// SetEqs sets equations
func (o *ElemP) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	o.Pmap = make([]int, o.Np)
	for m := 0; m < o.Shp.Nverts; m++ {
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

// InterpStarVars interpolates star variables to integration points
func (o *ElemP) InterpStarVars(sol *Solution) (ok bool) {

	// for each integration point
	for idx, ip := range o.IpsElem {

		// interpolation functions and gradients
		if LogErr(o.Shp.CalcAtIp(o.X, ip, true), "InterpStarVars") {
			return
		}

		// interpolate starred variables
		o.ψl[idx] = 0
		for m := 0; m < o.Shp.Nverts; m++ {
			o.ψl[idx] += o.Shp.S[m] * sol.Psi[o.Pmap[m]]
		}
	}
	return true
}

// AddToRhs adds -R to global residual vector fb
func (o ElemP) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// clear variables
	if o.DoExtrap {
		la.VecFill(o.ρl_ex, 0)
	}

	// for each integration point
	β1 := Global.DynCoefs.β1
	nverts := o.Shp.Nverts
	var coef, plt, klr, RhoL, ρl, Cpl float64
	var err error
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef = o.Shp.J * ip.W
		S := o.Shp.S
		G := o.Shp.G

		// tpm variables
		plt = β1*o.pl - o.ψl[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].Sl)
		RhoL = o.States[idx].RhoL
		ρl, Cpl, err = o.States[idx].Lvars(o.Mdl)
		if LogErr(err, "calc of tpm variables failed") {
			return
		}

		// compute ρwl. see Eq. (6) of [1]
		for i := 0; i < o.Ndim; i++ {
			o.ρwl[i] = 0
			for j := 0; j < o.Ndim; j++ {
				o.ρwl[i] += klr * o.Mdl.Klsat[i][j] * (RhoL*o.g[i] - o.gpl[i]) // TODO: fix this
			}
		}

		// add negative of residual term to fb. see Eqs. (12) and (17) of [1]
		for m := 0; m < nverts; m++ {
			r := o.Pmap[m]
			fb[r] -= coef * S[m] * Cpl * plt
			for i := 0; i < o.Ndim; i++ {
				fb[r] += coef * G[m][i] * o.ρwl[i] // += coef * div(ρl*wl)
			}
			if o.DoExtrap { // Eq. (19)
				o.ρl_ex[m] += o.Emat[m][idx] * ρl
			}
		}
		//io.Pf(">>> %3d : pc=%13.10f sl=%13.10f\n", o.Cell.Id, o.States[idx].Pg-o.States[idx].Pl, o.States[idx].Sl)
	}

	// contribution from natural boundary conditions
	if len(o.NatBcs) > 0 {
		return o.add_natbcs_to_rhs(fb, sol)
	}
	return true
}

// AddToKb adds element K to global Jacobian matrix Kb
func (o ElemP) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// clear matrices
	la.MatFill(o.Kpp, 0)
	nverts := o.Shp.Nverts
	if o.DoExtrap {
		for i := 0; i < nverts; i++ {
			o.ρl_ex[i] = 0
			for j := 0; j < nverts; j++ {
				o.dρldpl_ex[i][j] = 0
			}
		}
	}

	// for each integration point
	Cl := o.Mdl.Cl
	β1 := Global.DynCoefs.β1
	var coef, plt, klr, RhoL, ρl, Cpl, dCpldpl, dklrdpl float64
	var err error
	for idx, ip := range o.IpsElem {

		// interpolation functions, gradients and variables @ ip
		if !o.ipvars(idx, sol) {
			return
		}
		coef = o.Shp.J * ip.W
		S := o.Shp.S
		G := o.Shp.G

		// tpm variables
		plt = β1*o.pl - o.ψl[idx]
		klr = o.Mdl.Cnd.Klr(o.States[idx].Sl)
		RhoL = o.States[idx].RhoL
		ρl, Cpl, dCpldpl, dklrdpl, err = o.States[idx].Lderivs(o.Mdl)
		if LogErr(err, "calc of tpm derivatives failed") {
			return
		}

		// Kpp := dRpl/dpl. see Eqs. (18), (A.2) and (A.3) of [1]
		for n := 0; n < nverts; n++ {
			for j := 0; j < o.Ndim; j++ {
				o.tmp[j] = S[n]*dklrdpl*(RhoL*o.g[j]-o.gpl[j]) + klr*(S[n]*Cl*o.g[j]-G[n][j])
			}
			for m := 0; m < nverts; m++ {
				o.Kpp[m][n] += coef * S[m] * S[n] * (dCpldpl*plt + β1*Cpl)
				for i := 0; i < o.Ndim; i++ {
					for j := 0; j < o.Ndim; j++ {
						o.Kpp[m][n] -= coef * G[m][i] * o.Mdl.Klsat[i][j] * o.tmp[j]
					}
				}
				if o.DoExtrap { // inner summation term in Eq. (22)
					o.dρldpl_ex[m][n] += o.Emat[m][idx] * Cpl * S[n]
				}
			}
			if o.DoExtrap { // Eq. (19)
				o.ρl_ex[n] += o.Emat[n][idx] * ρl
			}
		}
	}

	// add to Kb
	if o.HasSeep {

		// contribution from natural boundary conditions
		if !o.add_natbcs_to_jac(sol) {
			return
		}

		// add to sparse matrix Kb
		for i, I := range o.Pmap {
			for j, J := range o.Pmap {
				Kb.Put(I, J, o.Kpp[i][j])
			}
			for j, J := range o.Fmap {
				Kb.Put(I, J, o.Kpf[i][j])
				Kb.Put(J, I, o.Kfp[j][i])
			}
		}
		for i, I := range o.Fmap {
			for j, J := range o.Fmap {
				Kb.Put(I, J, o.Kff[i][j])
			}
		}

	} else {

		// add to sparse matrix Kb
		for i, I := range o.Pmap {
			for j, J := range o.Pmap {
				Kb.Put(I, J, o.Kpp[i][j])
			}
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
		if LogErr(o.Shp.CalcAtIp(o.X, o.IpsElem[idx], false), "Update") { // TODO: no need for false here
			return
		}

		// compute Δpl @ ip by means of interpolating from nodes
		Δpl = 0
		for m := 0; m < o.Shp.Nverts; m++ {
			r := o.Pmap[m]
			Δpl += o.Shp.S[m] * sol.ΔY[r]
		}

		// update state
		if LogErr(o.Mdl.Update(o.States[idx], Δpl, 0, 0), "update failed") {
			return
		}
		//io.Pf("%3d : Δpl=%13.10f pc=%13.10f sl=%13.10f RhoL=%13.10f Wet=%v\n", o.Cell.Id, Δpl, o.States[idx].Pg-o.States[idx].Pl, o.States[idx].Sl, o.States[idx].RhoL, o.States[idx].Wet)
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
		if LogErr(o.Shp.CalcAtIp(o.X, o.IpsElem[idx], false), "InitIvs") {
			return
		}

		// compute pl @ ip by means of interpolating from nodes
		o.pl = 0
		for m := 0; m < o.Shp.Nverts; m++ {
			r := o.Pmap[m]
			o.pl += o.Shp.S[m] * sol.Y[r]
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
func (o *ElemP) SetIvs(ivs map[string][]float64) (ok bool) {
	if pl, ok := ivs["pl"]; ok {
		for idx, _ := range o.IpsElem {
			o.States[idx].Pl = pl[idx]
		}
	}
	if pg, ok := ivs["pg"]; ok {
		for idx, _ := range o.IpsElem {
			o.States[idx].Pg = pg[idx]
		}
	}
	for idx, _ := range o.IpsElem {
		pl := o.States[idx].Pl
		pg := o.States[idx].Pg
		if LogErr(o.Mdl.InitState(o.States[idx], pl, pg, 0), "state initialisation failed") {
			return
		}
	}
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
func (o ElemP) OutIpsData() (data []*OutIpData) {
	for idx, ip := range o.IpsElem {
		s := o.States[idx]
		x := o.Shp.IpRealCoords(o.X, ip)
		v := map[string]*float64{"sl": &s.Sl}
		data = append(data, &OutIpData{o.Id(), x, v})
	}
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *ElemP) ipvars(idx int, sol *Solution) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.Shp.CalcAtIp(o.X, o.IpsElem[idx], true), "ipvars") {
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
	for m := 0; m < o.Shp.Nverts; m++ {
		r := o.Pmap[m]
		o.pl += o.Shp.S[m] * sol.Y[r]
		for i := 0; i < o.Ndim; i++ {
			o.gpl[i] += o.Shp.G[m][i] * sol.Y[r]
		}
	}
	return true
}

// fipvars computes current values @ face integration points
func (o *ElemP) fipvars(fidx int, sol *Solution) (ρl, z, pl, fl float64) {
	Sf := o.Shp.Sf
	iz := o.Ndim - 1 // index of z-coordinate (elevation)
	for i, m := range o.Shp.FaceLocalV[fidx] {
		μ := o.Vid2seepId[m]
		ρl += Sf[i] * o.ρl_ex[m]
		z += Sf[i] * o.X[iz][m]
		pl += Sf[i] * sol.Y[o.Pmap[m]]
		fl += Sf[i] * sol.Y[o.Fmap[μ]]
	}
	return
}

// ramp implements the ramp function
func (o *ElemP) ramp(x float64) float64 {
	if o.Macaulay {
		return fun.Ramp(x)
	}
	return fun.Sramp(x, o.βrmp)
}

// rampderiv returns the ramp function first derivative
func (o *ElemP) rampD1(x float64) float64 {
	if o.Macaulay {
		return fun.Heav(x)
	}
	return fun.SrampD1(x, o.βrmp)
}

// add_natbcs_to_rhs adds natural boundary conditions to rhs
func (o ElemP) add_natbcs_to_rhs(fb []float64, sol *Solution) (ok bool) {

	// compute surface integral
	var tmp float64
	var ρl, z, pl, fl, plmax, g, rmp, rx, rf float64
	for _, nbc := range o.NatBcs {

		// temporary value == function evaluation
		tmp = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			if LogErr(o.Shp.CalcAtFaceIp(o.X, ipf, iface), "add_natbcs_to_rhs") {
				return
			}
			Sf := o.Shp.Sf

			// select natural boundary condition type
			switch nbc.Key {
			case "seepP", "seepH":

				// variables extrapolated to face
				ρl, z, pl, fl = o.fipvars(iface, sol)

				// specified condition
				plmax = tmp
				if nbc.Key == "seepH" {
					plmax = max(o.γl*(tmp-z), 0.0)
				}

				// compute residuals
				g = pl - plmax // Eq. (24)
				rmp = o.ramp(fl + o.κ*g)
				rx = ρl * rmp // Eq. (30)
				rf = fl - rmp // Eq. (26)
				if math.Abs(g) < 1e-13 {
					g = 0
				}
				for i, m := range o.Shp.FaceLocalV[iface] {
					μ := o.Vid2seepId[m]
					fb[o.Pmap[m]] -= ipf.W * Sf[i] * rx
					fb[o.Fmap[μ]] -= ipf.W * Sf[i] * rf
				}
			}
		}
	}
	return true
}

// add_natbcs_to_jac adds contribution from natural boundary conditions to Jacobian
func (o ElemP) add_natbcs_to_jac(sol *Solution) (ok bool) {

	// clear matrices
	if o.HasSeep {
		for i := 0; i < o.Np; i++ {
			for j := 0; j < o.Nf; j++ {
				o.Kpf[i][j] = 0
				o.Kfp[j][i] = 0
			}
		}
		la.MatFill(o.Kff, 0)
	}

	// compute surface integral
	nverts := o.Shp.Nverts
	var tmp float64
	var ρl, z, pl, fl, plmax, g, rmp, rmpD float64
	var drxdpl, drxdfl, drfdpl, drfdfl float64
	for _, nbc := range o.NatBcs {

		// temporary value == function evaluation
		tmp = nbc.Fcn.F(sol.T, nil)

		// loop over ips of face
		for _, ipf := range o.IpsFace {

			// interpolation functions and gradients @ face
			iface := nbc.IdxFace
			if LogErr(o.Shp.CalcAtFaceIp(o.X, ipf, iface), "add_natbcs_to_jac") {
				return
			}
			Sf := o.Shp.Sf

			// select natural boundary condition type
			switch nbc.Key {
			case "seepP", "seepH":

				// variables extrapolated to face
				ρl, z, pl, fl = o.fipvars(iface, sol)

				// specified condition
				plmax = tmp
				if nbc.Key == "seepH" {
					plmax = max(o.γl*(tmp-z), 0.0)
				}

				// compute derivatives
				g = pl - plmax // Eq. (24)
				rmp = o.ramp(fl + o.κ*g)
				rmpD = o.rampD1(fl + o.κ*g)
				drxdpl = ρl * o.κ * rmpD // first term in Eq. (A.4) (without Sn)
				drxdfl = ρl * rmpD       // Eq. (A.5) (without Sn)
				drfdpl = -o.κ * rmpD     // Eq. (A.6) (corrected with κ and without Sn)
				drfdfl = 1.0 - rmpD      // Eq. (A.7) (without Sn)
				for i, m := range o.Shp.FaceLocalV[iface] {
					μ := o.Vid2seepId[m]
					for j, n := range o.Shp.FaceLocalV[iface] {
						ν := o.Vid2seepId[n]
						o.Kpp[m][n] += ipf.W * Sf[i] * Sf[j] * drxdpl
						o.Kpf[m][ν] += ipf.W * Sf[i] * Sf[j] * drxdfl
						o.Kfp[μ][n] += ipf.W * Sf[i] * Sf[j] * drfdpl
						o.Kff[μ][ν] += ipf.W * Sf[i] * Sf[j] * drfdfl
					}
					for n := 0; n < nverts; n++ { // Eqs. (18) and (22)
						for l, r := range o.Shp.FaceLocalV[iface] {
							o.Kpp[m][n] += ipf.W * Sf[i] * Sf[l] * o.dρldpl_ex[r][n] * rmp
						}
					}
				}
			}
		}
	}
	return true
}
