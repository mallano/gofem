// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/msolid"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
	"github.com/cpmech/gosl/utl"
)

// Rjoint implements the rod-joint (interface/link) element for reinforced solids.
//  The following convention is considered:
//   n or N   -- means [N]odes
//   p or P   -- means integratioin [P]oints
//   nn or Nn -- number of nodes
//   np or Np -- number of integration [P]points
//   ndim     -- space dimension
//   nsig     -- number of stress/strain components == 2 * ndim
//   rod      -- means rod element
//   rodH     -- rod shape structure
//   rodNn    -- rod number of nodes
//   rodNp    -- rod number of integration points
//   rodS     -- rod shape functions
//   sld      -- means solid element
//   sldH     -- rod shape structure
//   sldNn    -- solid number of nodes
//   sldNp    -- solid number of integration points
//   sldS     -- solid shape functions
//   rodYn    -- rod's (real) coordinates of node
//   rodYp    -- rod's (real) coordinates of integration point
//   r or R   -- means natural coordinates in the solids' system
//   z or Z   -- means natural coordinates in the rod's system
//   s or S   -- parametric coordinate along rod
//   rodRn    -- natural coordinates or rod's nodes w.r.t solid's system
//   rodRp    -- natural coordinates of rod's integration point w.r.t to solid's system
//   Nmat     -- solid shape functions evaluated at rod nodes
//   Pmat     -- solid shape functions evaluated at rod integration points
//  References:
//   [1] R Durand, MM Farias, DM Pedroso. Modelling the strengthening of solids with
//       incompatible line finite elements, Computers and Structures (2014). Submitted.
//   [2] R Durand, MM Farias, DM Pedroso, Computing intersections between non-compatible
//       curves and finite elements, Computational Mechanics (2014). Submitted.
//   [3] R Durand, MM Farias. A local extrapolation method for finite elements,
//       Advances in Engineering Software 67 (2014) 1-9.
//       http://dx.doi.org/10.1016/j.advengsoft.2013.07.002
type Rjoint struct {

	// basic data
	Edat *inp.ElemData // element data; stored in allocator to be used in Connect
	Cid  int           // cell/element id
	Ny   int           // total number of dofs == rod.Nu + sld.Nu

	// essential
	Rod *Rod            // rod element
	Sld *ElemU          // solid element
	Mdl msolid.RjointM1 // material model

	// parameters
	h  float64 // perimeter of rod element; Eq (34)
	k1 float64 // lateral stiffness; Eq (37)
	k2 float64 // lateral stiffness; Eq (37)

	// optional data
	Ncns bool // use non-consistent model

	// shape functions evaluations and extrapolator matrices
	Nmat [][]float64 // [sldNn][rodNn] shape functions of solids @ [N]odes of rod element
	Pmat [][]float64 // [sldNn][rodNp] shape functions of solids @ integration [P]oints of rod element (for Coulomb model)
	Emat [][]float64 // [sldNn][sldNp] solid's extrapolation matrix (for Coulomb model)

	// variables for Coulomb model
	Coulomb bool            // use Coulomb model
	rodRp   [][]float64     // [rodNp][ndim] natural coordinates of ips of rod w.r.t. solid's system
	σNo     [][]float64     // [nneSld][nsig] σ at nodes of solid
	σIp     []float64       // [nsig] σ at ips of rod
	t1      []float64       // [ndim] traction vectors for σc
	t2      []float64       // [ndim] traction vectors for σc
	T1      [][]float64     // [rodNp][nsig] tensor (e1 dy e1)
	T2      [][]float64     // [rodNp][nsig] tensor (e2 dy e2)
	DσNoDu  [][][][]float64 // [sldNn][nsig][sldNn][ndim] ∂σSldNod/∂uSldNod : derivatives of σ @ nodes of solid w.r.t displacements of solid
	DσDun   [][]float64     // [nsig][ndim] ∂σIp/∂us : derivatives of σ @ ip of solid w.r.t displacements of solid

	// corotational system aligned with rod element
	e0 [][]float64 // [rodNp][ndim] local directions at each integration point of rod
	e1 [][]float64 // [rodNp][ndim] local directions at each integration point of rod
	e2 [][]float64 // [rodNp][ndim] local directions at each integration point of rod

	// auxiliary variables
	ΔuC [][]float64 // [rodNn][ndim] relative displ. increment of solid @ nodes of rod; Eq (30)
	Δw  []float64   // [ndim] relative velocity; Eq (32)
	qb  []float64   // [ndim] resultant traction vector 'holding' the rod @ ip; Eq (34)
	fC  []float64   // [rodNu] internal/contact forces vector; Eq (34)

	// temporary Jacobian matrices. see Eq. (57)
	Krr [][]float64 // [rodNu][rodNu] Eq. (58)
	Krs [][]float64 // [rodNu][sldNu] Eq. (59)
	Ksr [][]float64 // [sldNu][rodNu] Eq. (60)
	Kss [][]float64 // [sldNu][sldNu] Eq. (61)

	// internal values
	States    []*msolid.OnedState // [nip] internal states
	StatesBkp []*msolid.OnedState // [nip] backup internal states
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	infogetters["rjoint"] = func(cellType string, faceConds []*FaceCond) *Info {
		return &Info{}
	}

	// element allocator
	eallocators["rjoint"] = func(cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem {
		var o Rjoint
		o.Edat = edat
		o.Cid = cid
		return &o
	}
}

// Id returns the cell Id
func (o Rjoint) Id() int { return o.Cid }

// Connect connects rod/solid elements in this Rjoint
func (o *Rjoint) Connect(cid2elem []Elem, c *inp.Cell) (nnzK int, ok bool) {

	// get rod and solid elements
	rodId := c.JlinId
	sldId := c.JsldId
	o.Rod = cid2elem[rodId].(*Rod)
	o.Sld = cid2elem[sldId].(*ElemU)
	if LogErrCond(o.Rod == nil, "cannot find joint's rod cell with id == %d", rodId) {
		return
	}
	if LogErrCond(o.Sld == nil, "cannot find joint's solid cell with id == %d", sldId) {
		return
	}

	// total number of dofs
	o.Ny = o.Rod.Nu + o.Sld.Nu

	// material model name
	matname := o.Edat.Mat
	matdata := Global.Sim.Mdb.Get(matname)
	if LogErrCond(matdata == nil, "materials database failed on getting %q material\n", matname) {
		return
	}

	// initialise model
	if LogErr(o.Mdl.Init(matdata.Prms), "cannot initialise model for Rjoint element") {
		return
	}

	// parameters
	for _, p := range matdata.Prms {
		switch p.N {
		case "h":
			o.h = p.V
		case "k1":
			o.k1 = p.V
		case "k2":
			o.k2 = p.V
		case "mu":
			if p.V > 0.0 {
				o.Coulomb = true
			}
		}
	}

	// auxiliary
	ndim := Global.Ndim
	nsig := 2 * ndim

	// rod data
	rodH := o.Rod.Shp
	rodNp := len(o.Rod.IpsElem)
	rodNn := rodH.Nverts
	rodNu := o.Rod.Nu

	// solid data
	sldH := o.Sld.Shp
	sldS := sldH.S
	sldNp := len(o.Sld.IpsElem)
	sldNn := sldH.Nverts
	sldNu := o.Sld.Nu

	// shape functions of solid @ nodes of rod
	o.Nmat = la.MatAlloc(sldNn, rodNn)
	rodYn := make([]float64, ndim)
	rodRn := make([]float64, 3)
	for m := 0; m < rodNn; m++ {
		for i := 0; i < ndim; i++ {
			rodYn[i] = o.Rod.X[i][m]
		}
		if LogErr(sldH.InvMap(rodRn, rodYn, o.Sld.X), "inverse map failed") {
			return
		}
		if LogErr(sldH.CalcAtR(o.Sld.X, rodRn, false), "shape functions calculation failed") {
			return
		}
		for n := 0; n < sldNn; n++ {
			o.Nmat[n][m] = sldH.S[n]
		}
	}

	// coulomb model => σc depends on p values of solid
	if o.Coulomb {

		// allocate variables
		o.Pmat = la.MatAlloc(sldNn, rodNp)
		o.Emat = la.MatAlloc(sldNn, sldNp)
		o.rodRp = la.MatAlloc(rodNp, 3)
		o.σNo = la.MatAlloc(sldNn, nsig)
		o.σIp = make([]float64, nsig)
		o.t1 = make([]float64, ndim)
		o.t2 = make([]float64, ndim)
		o.T1 = la.MatAlloc(rodNp, nsig)
		o.T2 = la.MatAlloc(rodNp, nsig)
		o.DσNoDu = utl.Deep4alloc(sldNn, nsig, sldNn, ndim)
		o.DσDun = la.MatAlloc(nsig, ndim)

		// extrapolator matrix
		if LogErr(sldH.Extrapolator(o.Emat, o.Sld.IpsElem), "Extrapolator of solid failed") {
			return
		}

		// shape function of solid @ ips of rod
		for idx, ip := range o.Rod.IpsElem {
			rodYp := rodH.IpRealCoords(o.Rod.X, ip)
			if LogErr(sldH.InvMap(o.rodRp[idx], rodYp, o.Sld.X), "inverse map failed") {
				return
			}
			if LogErr(sldH.CalcAtR(o.Sld.X, o.rodRp[idx], false), "shape functions calculation failed") {
				return
			}
			for n := 0; n < sldNn; n++ {
				o.Pmat[n][idx] = sldS[n]
			}
		}
	}

	// joint direction @ ip[idx]; corotational system aligned with rod element
	o.e0 = la.MatAlloc(rodNp, ndim)
	o.e1 = la.MatAlloc(rodNp, ndim)
	o.e2 = la.MatAlloc(rodNp, ndim)
	π := make([]float64, ndim) // Eq. (27)
	Q := la.MatAlloc(ndim, ndim)
	α := 666.0
	Jvec := rodH.Jvec3d[:ndim]
	for idx, ip := range o.Rod.IpsElem {

		// auxiliary
		e0, e1, e2 := o.e0[idx], o.e1[idx], o.e2[idx]

		// interpolation functions and gradients
		if LogErr(rodH.CalcAtIp(o.Rod.X, ip, true), "shape functions calculation failed") {
			return
		}

		// compute basis vectors
		J := rodH.J
		π[0] = Jvec[0] + α
		π[1] = Jvec[1]
		e0[0] = Jvec[0] / J
		e0[1] = Jvec[1] / J
		if ndim == 3 {
			π[2] = Jvec[2]
			e0[2] = Jvec[2] / J
		}
		la.MatSetDiag(Q, 1)
		la.VecOuterAdd(Q, -1, e0, e0) // Q := I - e0 dyad e0
		la.MatVecMul(e1, 1, Q, π)     // Eq. (29) * norm(E1)
		la.VecScale(e1, 0, 1.0/la.VecNorm(e1), e1)
		if ndim == 3 {
			e2[0] = e0[1]*e1[2] - e0[2]*e1[1]
			e2[1] = e0[2]*e1[0] - e0[0]*e1[2]
			e2[2] = e0[0]*e1[1] - e0[1]*e1[0]
		}

		// compute auxiliary tensors
		if o.Coulomb {
			e1_dy_e1 := tsr.Alloc2()
			e2_dy_e2 := tsr.Alloc2()
			for i := 0; i < ndim; i++ {
				for j := 0; j < ndim; j++ {
					e1_dy_e1[i][j] = e1[i] * e1[j]
					e2_dy_e2[i][j] = e2[i] * e2[j]
				}
			}
			tsr.Ten2Man(o.T1[idx], e1_dy_e1)
			tsr.Ten2Man(o.T2[idx], e2_dy_e2)
		}
	}

	// auxiliary variables
	o.ΔuC = la.MatAlloc(rodNn, ndim)
	o.Δw = make([]float64, ndim)
	o.qb = make([]float64, ndim)
	o.fC = make([]float64, rodNu)

	// temporary Jacobian matrices. see Eq. (57)
	o.Krr = la.MatAlloc(rodNu, rodNu)
	o.Krs = la.MatAlloc(rodNu, sldNu)
	o.Ksr = la.MatAlloc(sldNu, rodNu)
	o.Kss = la.MatAlloc(sldNu, sldNu)

	// debugging
	//if true {
	if false {
		o.debug_print_init()
	}

	// success
	return o.Ny * o.Ny, true
}

// implementation ///////////////////////////////////////////////////////////////////////////////////

// SetEqs set equations
func (o *Rjoint) SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) {
	return true
}

// SetEleConds set element conditions
func (o *Rjoint) SetEleConds(key string, f fun.Func, extra string) (ok bool) {
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *Rjoint) InterpStarVars(sol *Solution) (ok bool) {
	return true
}

// adds -R to global residual vector fb
func (o *Rjoint) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// auxiliary
	ndim := Global.Ndim
	rodH := o.Rod.Shp
	rodS := rodH.S
	rodNn := rodH.Nverts
	sldH := o.Sld.Shp
	sldNn := sldH.Nverts

	// internal forces vector
	la.VecFill(o.fC, 0)

	// loop over rod's integration points
	var coef, τ, qn1, qn2 float64
	for idx, ip := range o.Rod.IpsElem {

		// auxiliary
		e0, e1, e2 := o.e0[idx], o.e1[idx], o.e2[idx]

		// interpolation functions and gradients
		if LogErr(rodH.CalcAtIp(o.Rod.X, ip, true), "AddToRhs") {
			return
		}
		coef = ip.W * rodH.J

		// state variables
		τ = o.States[idx].Sig
		qn1 = o.States[idx].Phi[0]
		qn2 = o.States[idx].Phi[1]

		// fC vector. Eq. (34)
		for i := 0; i < ndim; i++ {
			o.qb[i] = τ*o.h*e0[i] + qn1*e1[i] + qn2*e2[i]
			for m := 0; m < rodNn; m++ {
				r := i + m*ndim
				o.fC[r] += coef * rodS[m] * o.qb[i]
			}
		}
	}

	// fb = -Resid;  fR = -fC  and  fS = Nmat*fC  =>  fb := {fC, -Nmat*fC}
	for i := 0; i < ndim; i++ {
		for m := 0; m < rodNn; m++ {
			r := i + m*ndim
			I := o.Rod.Umap[r]
			fb[I] += o.fC[r] // fb := - (fR == -fC Eq (35))
			for n := 0; n < sldNn; n++ {
				s := i + n*ndim
				J := o.Sld.Umap[s]
				fb[J] -= o.Nmat[n][m] * o.fC[r] // fb := - (fS Eq (36))
			}
		}
	}
	return true
}

// adds element K to global Jacobian matrix Kb
func (o *Rjoint) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {

	// auxiliary
	ndim := Global.Ndim
	nsig := 2 * ndim
	rodH := o.Rod.Shp
	rodS := rodH.S
	rodNn := rodH.Nverts
	sldH := o.Sld.Shp
	sldNn := sldH.Nverts

	// compute DσNoDu
	if o.Coulomb {

		// clear deep4 structure
		utl.Deep4set(o.DσNoDu, 0)

		// loop over solid's integration points
		for idx, ip := range o.Sld.IpsElem {

			// interpolation functions, gradients and variables @ ip
			if LogErr(sldH.CalcAtIp(o.Sld.X, ip, true), "AddToKb") {
				return
			}

			// consistent tangent model matrix
			if LogErr(o.Sld.MdlSmall.CalcD(o.Sld.D, o.Sld.States[idx], firstIt), "AddToKb") {
				return
			}

			// extrapolate derivatives
			for n := 0; n < sldNn; n++ {
				DerivSig(o.DσDun, n, ndim, sldH.G, o.Sld.D)
				for m := 0; m < sldNn; m++ {
					for i := 0; i < nsig; i++ {
						for j := 0; j < ndim; j++ {
							o.DσNoDu[m][i][n][j] += o.Emat[m][idx] * o.DσDun[i][j]
						}
					}
				}
			}
		}
	}

	// debugging
	//if true {
	if false {
		utl.PrintDeep4("DσNoDu", o.DσNoDu, "%20.10f")
	}

	// zero K matrices
	for i, _ := range o.Rod.Umap {
		for j, _ := range o.Rod.Umap {
			o.Krr[i][j] = 0
		}
		for j, _ := range o.Sld.Umap {
			o.Krs[i][j] = 0
			o.Ksr[j][i] = 0
		}
	}
	la.MatFill(o.Kss, 0)

	// auxiliary
	var coef float64
	var DτDω, DτDσc, DσcDu_nj float64
	var Dp1Du_nj, Dp2Du_nj float64
	var Dwb0Du_nj, Dwb1Du_nj, Dwb2Du_nj float64
	var DτDu_nj, DqbDu_nij float64
	var Dwb0Dur_nj, Dwb1Dur_nj, Dwb2Dur_nj float64
	var DqbDur_nij float64

	// loop over rod's integration points
	var err error
	for idx, ip := range o.Rod.IpsElem {

		// auxiliary
		e0, e1, e2 := o.e0[idx], o.e1[idx], o.e2[idx]

		// interpolation functions and gradients
		if LogErr(rodH.CalcAtIp(o.Rod.X, ip, true), "AddToKb") {
			return
		}
		coef = ip.W * rodH.J

		// model derivatives
		DτDω, DτDσc, err = o.Mdl.CalcD(o.States[idx], firstIt)
		if LogErr(err, "AddToKb") {
			return
		}

		// compute derivatives
		for j := 0; j < ndim; j++ {

			// Krr and Ksr; derivatives with respect to ur_nj
			for n := 0; n < rodNn; n++ {

				// ∂wb/∂ur Eq (A.4)
				Dwb0Dur_nj = -rodS[n] * e0[j]
				Dwb1Dur_nj = -rodS[n] * e1[j]
				Dwb2Dur_nj = -rodS[n] * e2[j]

				// compute ∂■/∂ur derivatives
				c := j + n*ndim
				for i := 0; i < ndim; i++ {

					// ∂qb/∂ur Eq (A.2)
					DqbDur_nij = o.h*e0[i]*(DτDω*Dwb0Dur_nj) + o.k1*e1[i]*Dwb1Dur_nj + o.k2*e2[i]*Dwb2Dur_nj

					// Krr := ∂fr/∂ur Eq (58)
					for m := 0; m < rodNn; m++ {
						r := i + m*ndim
						o.Krr[r][c] -= coef * rodS[m] * DqbDur_nij
					}

					//  Ksr := ∂fs/∂ur Eq (60)
					for m := 0; m < sldNn; m++ {
						r := i + m*ndim
						for p := 0; p < rodNn; p++ {
							o.Ksr[r][c] += coef * o.Nmat[m][p] * rodS[p] * DqbDur_nij
						}
					}
				}
			}

			// Krs and Kss
			for n := 0; n < sldNn; n++ {

				// ∂σc/∂us_nj
				DσcDu_nj = 0
				if o.Coulomb {

					// Eqs (A.10) (A.11) and (A.12)
					Dp1Du_nj, Dp2Du_nj = 0, 0
					for m := 0; m < sldNn; m++ {
						for i := 0; i < nsig; i++ {
							Dp1Du_nj += o.Pmat[m][idx] * o.T1[idx][i] * o.DσNoDu[m][i][n][j]
							Dp2Du_nj += o.Pmat[m][idx] * o.T2[idx][i] * o.DσNoDu[m][i][n][j]
						}
					}
					DσcDu_nj = (Dp1Du_nj + Dp2Du_nj) / 2.0
				}

				// ∂wb/∂us Eq (A.5)
				Dwb0Du_nj, Dwb1Du_nj, Dwb2Du_nj = 0, 0, 0
				for m := 0; m < rodNn; m++ {
					Dwb0Du_nj += rodS[m] * o.Nmat[n][m] * e0[j]
					Dwb1Du_nj += rodS[m] * o.Nmat[n][m] * e1[j]
					Dwb2Du_nj += rodS[m] * o.Nmat[n][m] * e2[j]
				}

				// ∂τ/∂us_nj hightlighted term in Eq (A.3)
				DτDu_nj = DτDω*Dwb0Du_nj + DτDσc*DσcDu_nj
				if o.Ncns {
					DτDu_nj = DτDω * Dwb0Du_nj
				}

				// compute ∂■/∂us derivatives
				c := j + n*ndim
				for i := 0; i < ndim; i++ {

					// ∂qb/∂us Eq (A.3)
					DqbDu_nij = o.h*e0[i]*DτDu_nj + o.k1*e1[i]*Dwb1Du_nj + o.k2*e2[i]*Dwb2Du_nj

					// Krs := ∂fr/∂us Eq (59)
					for m := 0; m < rodNn; m++ {
						r := i + m*ndim
						o.Krs[r][c] -= coef * rodS[m] * DqbDu_nij
					}

					// Kss := ∂fs/∂us Eq (61)
					for m := 0; m < sldNn; m++ {
						r := i + m*ndim
						for p := 0; p < rodNn; p++ {
							o.Kss[r][c] += coef * o.Nmat[m][p] * rodS[p] * DqbDu_nij
						}
					}
				}
			}
		}
	}

	// debug
	//if true {
	if false {
		o.debug_print_K()
	}

	// add K to sparse matrix Kb
	for i, I := range o.Rod.Umap {
		for j, J := range o.Rod.Umap {
			Kb.Put(I, J, o.Krr[i][j])
		}
		for j, J := range o.Sld.Umap {
			Kb.Put(I, J, o.Krs[i][j])
			Kb.Put(J, I, o.Ksr[j][i])
		}
	}
	for i, I := range o.Sld.Umap {
		for j, J := range o.Sld.Umap {
			Kb.Put(I, J, o.Kss[i][j])
		}
	}
	return true
}

// Update perform (tangent) update
func (o *Rjoint) Update(sol *Solution) (ok bool) {

	// auxiliary
	ndim := Global.Ndim
	nsig := 2 * ndim
	rodH := o.Rod.Shp
	rodS := rodH.S
	rodNn := rodH.Nverts
	sldH := o.Sld.Shp
	sldNn := sldH.Nverts

	// extrapolate stresses at integration points of solid element to its nodes
	if o.Coulomb {
		la.MatFill(o.σNo, 0)
		for idx, _ := range o.Sld.IpsElem {
			σ := o.Sld.States[idx].Sig
			for i := 0; i < nsig; i++ {
				for m := 0; m < sldNn; m++ {
					o.σNo[m][i] += o.Emat[m][idx] * σ[i]
				}
			}
		}
	}

	// interpolate Δu of solid to find ΔuC @ rod node; Eq (30)
	var r, I int
	for m := 0; m < rodNn; m++ {
		for i := 0; i < ndim; i++ {
			o.ΔuC[m][i] = 0
			for n := 0; n < sldNn; n++ {
				r = i + n*ndim
				I = o.Sld.Umap[r]
				o.ΔuC[m][i] += o.Nmat[n][m] * sol.ΔY[I] // Eq (30)
			}
		}
	}

	// loop over ips of rod
	var Δwb0, Δwb1, Δwb2, σc float64
	for idx, ip := range o.Rod.IpsElem {

		// auxiliary
		e0, e1, e2 := o.e0[idx], o.e1[idx], o.e2[idx]

		// interpolation functions and gradients
		if LogErr(rodH.CalcAtIp(o.Rod.X, ip, true), "Update") {
			return
		}

		// interpolated relative displacements @ ip of join; Eqs (31) and (32)
		for i := 0; i < ndim; i++ {
			o.Δw[i] = 0
			for m := 0; m < rodNn; m++ {
				r = i + m*ndim
				I = o.Rod.Umap[r]
				o.Δw[i] += rodS[m] * (o.ΔuC[m][i] - sol.ΔY[I]) // Eq (31) and (32)
			}
		}

		// relative displacents in the coratational system
		Δwb0, Δwb1, Δwb2 = 0, 0, 0
		for i := 0; i < ndim; i++ {
			Δwb0 += e0[i] * o.Δw[i]
			Δwb1 += e1[i] * o.Δw[i]
			Δwb2 += e2[i] * o.Δw[i]
		}

		// new confining stress
		σc = 0.0
		if o.Coulomb {

			// calculate σIp
			for j := 0; j < nsig; j++ {
				o.σIp[j] = 0
				for n := 0; n < sldNn; n++ {
					o.σIp[j] += o.Pmat[n][idx] * o.σNo[n][j]
				}
			}

			// calculate t1 and t2
			for i := 0; i < ndim; i++ {
				o.t1[i], o.t2[i] = 0, 0
				for j := 0; j < ndim; j++ {
					o.t1[i] += tsr.M2T(o.σIp, i, j) * e1[j]
					o.t2[i] += tsr.M2T(o.σIp, i, j) * e2[j]
				}
			}

			// calculate p1, p2 and σcNew
			p1, p2 := 0.0, 0.0
			for i := 0; i < ndim; i++ {
				p1 += o.t1[i] * e1[i]
				p2 += o.t2[i] * e2[i]
			}

			// σcNew
			σc = -(p1 + p2) / 2.0
		}

		// update model
		if LogErr(o.Mdl.Update(o.States[idx], σc, Δwb0), "Update") {
			return
		}
		o.States[idx].Phi[0] += o.k1 * Δwb1 // qn1
		o.States[idx].Phi[1] += o.k2 * Δwb2 // qn2

		// debugging
		//if true {
		if false {
			o.debug_update(idx, Δwb0, Δwb1, Δwb2, σc)
		}
	}
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// Ipoints returns the real coordinates of integration points [nip][ndim]
func (o Rjoint) Ipoints() (coords [][]float64) {
	return o.Rod.Ipoints()
}

// SetIniIvs sets initial ivs for given values in sol and ivs map
func (o *Rjoint) SetIniIvs(sol *Solution, ivs map[string][]float64) (ok bool) {
	nip := len(o.Rod.IpsElem)
	o.States = make([]*msolid.OnedState, nip)
	o.StatesBkp = make([]*msolid.OnedState, nip)
	for i := 0; i < nip; i++ {
		o.States[i], _ = o.Mdl.InitIntVars()
		o.StatesBkp[i] = o.States[i].GetCopy()
	}
	return true
}

// BackupIvs create copy of internal variables
func (o *Rjoint) BackupIvs() (ok bool) {
	for i, s := range o.StatesBkp {
		s.Set(o.States[i])
	}
	return true
}

// RestoreIvs restore internal variables from copies
func (o *Rjoint) RestoreIvs() (ok bool) {
	for i, s := range o.States {
		s.Set(o.StatesBkp[i])
	}
	return true
}

// Ureset fixes internal variables after u (displacements) have been zeroed
func (o *Rjoint) Ureset(sol *Solution) (ok bool) {
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o Rjoint) Encode(enc Encoder) (ok bool) {
	return !LogErr(enc.Encode(o.States), "Encode")
}

// Decode decodes internal variables
func (o Rjoint) Decode(dec Decoder) (ok bool) {
	if LogErr(dec.Decode(&o.States), "Decode") {
		return
	}
	return o.BackupIvs()
}

// OutIpsData returns data from all integration points for output
func (o Rjoint) OutIpsData() (data []*OutIpData) {
	return
}

// debugging ////////////////////////////////////////////////////////////////////////////////////////

func (o Rjoint) debug_print_init() {
	sldNn := o.Sld.Shp.Nverts
	rodNn := o.Rod.Shp.Nverts
	rodNp := len(o.Rod.IpsElem)
	io.Pf("Nmat =\n")
	for i := 0; i < sldNn; i++ {
		for j := 0; j < rodNn; j++ {
			io.Pf("%g ", o.Nmat[i][j])
		}
		io.Pf("\n")
	}
	io.Pf("\nPmat =\n")
	for i := 0; i < sldNn; i++ {
		for j := 0; j < rodNp; j++ {
			io.Pf("%g ", o.Pmat[i][j])
		}
		io.Pf("\n")
	}
	io.Pf("\n")
	la.PrintMat("e0", o.e0, "%20.13f", false)
	io.Pf("\n")
	la.PrintMat("e1", o.e1, "%20.13f", false)
	io.Pf("\n")
	la.PrintMat("e2", o.e2, "%20.13f", false)
	io.Pf("\n")
	la.PrintMat("T1", o.T1, "%20.13f", false)
	io.Pf("\n")
	la.PrintMat("T2", o.T2, "%20.13f", false)
}

func (o Rjoint) debug_print_K() {
	ndim := Global.Ndim
	sldNn := o.Sld.Shp.Nverts
	rodNn := o.Rod.Shp.Nverts
	K := la.MatAlloc(o.Ny, o.Ny)
	start := o.Sld.Nu
	for i := 0; i < ndim; i++ {
		for m := 0; m < sldNn; m++ {
			r := i + m*ndim
			for j := 0; j < ndim; j++ {
				for n := 0; n < sldNn; n++ {
					c := j + n*ndim
					K[r][c] = o.Kss[r][c]
				}
				for n := 0; n < rodNn; n++ {
					c := j + n*ndim
					K[r][start+c] = o.Ksr[r][c]
					K[start+c][r] = o.Krs[c][r]
				}
			}
		}
	}
	for i := 0; i < ndim; i++ {
		for m := 0; m < rodNn; m++ {
			r := i + m*ndim
			for j := 0; j < ndim; j++ {
				for n := 0; n < rodNn; n++ {
					c := j + n*ndim
					K[start+r][start+c] = o.Krr[r][c]
				}
			}
		}
	}
	la.PrintMat("K", K, "%20.10f", false)
}

func (o Rjoint) debug_update(idx int, Δwb0, Δwb1, Δwb2, σc float64) {
	τ := o.States[idx].Sig
	qn1 := o.States[idx].Phi[0]
	qn2 := o.States[idx].Phi[1]
	la.PrintVec("Δw", o.Δw, "%13.10f", false)
	io.Pf("Δwb0=%13.10f Δwb1=%13.10f Δwb2=%13.10f\n", Δwb0, Δwb1, Δwb2)
	la.PrintVec("σIp", o.σIp, "%13.10f", false)
	io.Pf("σc=%13.10f t1=%13.10f t2=%13.10f\n", σc, o.t1, o.t2)
	io.Pf("τ=%13.10f qn1=%13.10f qn2=%13.10f\n", τ, qn1, qn2)
}
