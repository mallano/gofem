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
	"github.com/cpmech/gosl/utl"
)

// Rjoint implements the rod-joint (interface/link) element for reinforced solids.
//  The following convention is considered:
//   n or N    -- means [N]odes
//   p or P    -- means integratioin [P]oints
//   nn or Nn  -- number of nodes
//   np or Np  -- number of integration [P]points
//   ndim      -- space dimension
//   nsig      -- number of stress/strain components == 2 * ndim
//   rod       -- means rod element
//   rodH      -- rod shape structure
//   rodNn     -- rod number of nodes
//   rodNp     -- rod number of integration points
//   rodS      -- rod shape functions
//   sld       -- means solid element
//   sldH      -- rod shape structure
//   sldNn     -- solid number of nodes
//   sldNp     -- solid number of integration points
//   sldS      -- solid shape functions
//   rodYn     -- rod's (real) coordinates of node
//   rodYp     -- rod's (real) coordinates of integration point
//   r or R    -- means natural coordinates in the solids' system
//   z or Z    -- means natural coordinates in the rod's system
//   s or S    -- parametric coordinate along rod
//   rodRn     -- natural coordinates or rod's nodes w.r.t solid's system
//   rodRp     -- natural coordinates of rod's integration point w.r.t to solid's system
//   sldS_rodN -- solid shape functions evaluated at rod nodes
//   sldS_rodP -- solid shape functions evaluated at rot integration points
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
	Cell *inp.Cell     // cell
	Ndim int           // space dimension
	Ny   int           // total number of dofs == rod.Nu + sld.Nu

	// underlying elements
	Rod *Rod   // rod element
	Sld *ElemU // solid element

	// integration points
	IpsElem []*shp.Ipoint // [nip] integration points of (this) element

	// material model
	Mdl msolid.RjointM1 // material model

	// parameters
	h  float64 // perimeter of rod element
	k1 float64 // lateral stiffness
	k2 float64 // lateral stiffness

	// shape functions evaluations (see convention above)
	sldS_rodN [][]float64 // [rodNn][sldNn] shape functions of solids @ [N]odes of rod element

	// corotational system aligned with rod element
	e0 [][]float64 // [rodNp][ndim] local directions at each integration point of rod
	e1 [][]float64 // [rodNp][ndim] local directions at each integration point of rod
	e2 [][]float64 // [rodNp][ndim] local directions at each integration point of rod

	// variables for Coulomb model
	Coulomb   bool        // use Coulomb model
	sldEmat   [][]float64 // [sldNn][sldNp] solid's extrapolation matrix
	sldSigN   [][]float64 // [sldNn][nsig] σ at nodes of solid (extrapolated using Emat)
	rodRp     [][]float64 // [rodNp][ndim] natural coordinates of ips of rod w.r.t. solid's system
	sldS_rodP [][]float64 // [rodNp][sldNn] shape functions of solids @ integration [P]oints of rod element

	// scratchpad for Coulomb model. computed @ each ip
	sldSigN_rodP []float64       // [nsig] σ of solid interpolated to ips of rod
	t1           []float64       // [ndim] traction vectors for σc
	t2           []float64       // [ndim] traction vectors for σc
	T1           [][]float64     // [rodNp][nsig] tensor (e1 dy e1)
	T2           [][]float64     // [rodNp][nsig] tensor (e2 dy e2)
	DσNoDu       [][][][]float64 // [sldNn][nsig][sldNn][ndim] ∂σSldNod/∂uSldNod : derivatives of σ @ nodes of solid w.r.t displacements of solid
	DσDun        [][]float64     // [nsig][ndim] ∂σSldIp/∂uSldNod : derivatives of σ @ ip of solid w.r.t displacements of solid

	// problem variables
	Ymap []int // assembly map (location array/element equations)

	// internal variables
	States    []*msolid.OnedState // [nip] internal states
	StatesBkp []*msolid.OnedState // [nip] backup internal states

	// scratchpad. computed @ each ip
	Δw   []float64 // relative velocity
	qb   []float64 // resultant traction vector 'holding' the rod @ ip
	σc   float64   // confining stress
	Δwb0 float64   // increment of relative displacement along 0-direction
	Δwb1 float64   // increment of relative displacement along 1-direction
	Δwb2 float64   // increment of relative displacement along 2-direction
	τ    float64   // shear stress along interface between rod and solid
	qn1  float64   // orthogonal distributed loads acting on the perimeter of the rod along 1-direction
	qn2  float64   // orthogonal distributed loads acting on the perimeter of the rod along 2-direction
}

// initialisation ///////////////////////////////////////////////////////////////////////////////////

// register element
func init() {

	// information allocator
	iallocators["rjoint"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info {
		return &Info{}
	}

	// element allocator
	eallocators["rjoint"] = func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {
		var o Rjoint
		o.Edat = edat
		o.Cell = msh.Cells[cid]
		o.Ndim = msh.Ndim
		return &o
	}
}

// connect rod and solid ////////////////////////////////////////////////////////////////////////////

// Connect connects rod/solid elements in this Rjoint
func (o *Rjoint) Connect(cid2elem []Elem) (nnzK int, ok bool) {

	// get rod and solid elements
	rodId := o.Cell.JlinId
	sldId := o.Cell.JsldId
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

	// integration points; based on Shape of rod element
	var err error
	var nip int
	if s_nip, found := io.Keycode(o.Edat.Extra, "nip"); found {
		nip = io.Atoi(s_nip)
	}
	o.IpsElem, err = shp.GetIps(o.Rod.Cell.Shp.Type, nip)
	if LogErr(err, "GetIps failed element") {
		return
	}

	// material model name
	matname := o.Edat.Mat
	matdata := Global.Mdb.Get(matname)
	if LogErrCond(matdata == nil, "Mdb.Get failed\n") {
		return
	}

	// initialise model
	o.Mdl.Init(matdata.Prms)

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
	ndim := o.Ndim
	nsig := 2 * o.Ndim

	// rod data
	rodH := o.Rod.Cell.Shp
	rodNp := len(o.Rod.IpsElem)
	rodNn := rodH.Nverts

	// solid data
	sldH := o.Sld.Cell.Shp
	sldS := sldH.S
	sldNp := len(o.Sld.IpsElem)
	sldNn := sldH.Nverts

	// shape functions of solid @ nodes of rod
	o.sldS_rodN = la.MatAlloc(rodNn, sldNn)
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
			o.sldS_rodN[m][n] = sldH.S[n]
		}
	}

	// coulomb model => σc depends on p values of solid
	if o.Coulomb {

		// allocate variables
		o.sldEmat = la.MatAlloc(sldNn, sldNp)
		o.sldSigN = la.MatAlloc(sldNn, nsig)
		o.rodRp = la.MatAlloc(rodNp, 3)
		o.sldS_rodP = la.MatAlloc(rodNp, sldNn)

		// extrapolator matrix
		if LogErr(sldH.Extrapolator(o.sldEmat, o.Sld.IpsElem), "Extrapolator of solid failed") {
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
				o.sldS_rodP[idx][n] = sldS[n]
			}
		}

		// scratchpad for Coulomb model. computed @ each ip
		o.sldSigN_rodP = make([]float64, nsig)
		o.t1 = make([]float64, ndim)
		o.t2 = make([]float64, ndim)
		o.T1 = la.MatAlloc(rodNp, nsig)
		o.T2 = la.MatAlloc(rodNp, nsig)
		o.DσNoDu = utl.Deep4alloc(sldNn, nsig, sldNn, ndim)
		o.DσDun = la.MatAlloc(nsig, ndim)
	}

	// joint direction @ ip[idx]; corotational system aligned with rod element
	o.e0 = la.MatAlloc(rodNp, ndim)
	o.e1 = la.MatAlloc(rodNp, ndim)
	o.e2 = la.MatAlloc(rodNp, ndim)
	π := make([]float64, ndim) // Eq. (27)
	Q := la.MatAlloc(ndim, ndim)
	α := 666.0
	J := rodH.J
	Jvec := rodH.Jvec3d[:ndim]
	for idx, ip := range o.IpsElem {
		e0, e1, e2 := o.e0[idx], o.e1[idx], o.e2[idx]
		if LogErr(rodH.CalcAtIp(o.Rod.X, ip, true), "shape functions calculation failed") {
			return
		}
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

	// full Usmap
	o.Ymap = make([]int, o.Ny)
	for m := 0; m < sldNn; m++ {
		for i := 0; i < ndim; i++ {
			r := i + m*ndim
			o.Ymap[r] = o.Sld.Umap[r]
		}
	}
	start := o.Sld.Nu
	for m := 0; m < rodNn; m++ {
		for i := 0; i < ndim; i++ {
			r := i + m*ndim
			o.Ymap[start+r] = o.Rod.Umap[r]
		}
	}

	// scratchpad. computed @ each ip
	o.Δw = make([]float64, o.Ndim)
	o.qb = make([]float64, o.Ndim)

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

// SetSurfLoads set surface loads (natural boundary conditions)
func (o *Rjoint) SetNatBcs(key string, idxface int, f fun.Func, extra string) (ok bool) {
	return true
}

// InterpStarVars interpolates star variables to integration points
func (o *Rjoint) InterpStarVars(sol *Solution) (ok bool) {
	return true
}

// adds -R to global residual vector fb
func (o *Rjoint) AddToRhs(fb []float64, sol *Solution) (ok bool) {

	// auxiliary
	rodH := o.Rod.Cell.Shp
	sldH := o.Sld.Cell.Shp
	rodS := rodH.S
	rodNv := rodH.Nverts
	sldNv := sldH.Nverts

	start := o.Sld.Nu
	var coef float64
	var rsld, rlin int
	for idx, ip := range o.IpsElem {
		rodH.CalcAtIp(o.Rod.X, ip, true)

		// state variables
		o.τ = o.States[idx].Sig
		o.qn1 = o.States[idx].Phi[0]
		o.qn2 = o.States[idx].Phi[1]

		for i := 0; i < o.Ndim; i++ {
			o.qb[i] = o.τ*o.h*o.e0[idx][i] + o.qn1*o.e1[idx][i] + o.qn2*o.e2[idx][i]
		}
		coef = ip.W * rodH.J
		for m := 0; m < rodNv; m++ {
			for i := 0; i < o.Ndim; i++ {
				rlin = o.Ymap[start+i+m*o.Ndim]
				fb[rlin] += coef * rodS[m] * o.qb[i] // -fi  :  fR = -fC
			}
		}
	}

	for m := 0; m < rodNv; m++ {
		for i := 0; i < o.Ndim; i++ {
			rlin = start + i + m*o.Ndim
			for n := 0; n < sldNv; n++ {
				rsld = o.Ymap[i+n*o.Ndim]
				fb[rsld] += o.sldS_rodN[m][n] //* o.Rus[rlin] // -fi  :  f = fC
			}
		}
	}
	return true
}

// adds element K to global Jacobian matrix Kb
func (o *Rjoint) AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) {
	return true
}

// Update perform (tangent) update
func (o *Rjoint) Update(sol *Solution) (ok bool) {
	return true
}

// internal variables ///////////////////////////////////////////////////////////////////////////////

// InitIvs reset (and fix) internal variables after primary variables have been changed
func (o *Rjoint) InitIvs(sol *Solution) (ok bool) {
	return true
}

// SetIvs set secondary variables; e.g. during initialisation via files
func (o *Rjoint) SetIvs(zvars map[string][]float64) (ok bool) {
	return true
}

// BackupIvs create copy of internal variables
func (o *Rjoint) BackupIvs() (ok bool) {
	return true
}

// RestoreIvs restore internal variables from copies
func (o *Rjoint) RestoreIvs() (ok bool) {
	return true
}

// writer ///////////////////////////////////////////////////////////////////////////////////////////

// Encode encodes internal variables
func (o Rjoint) Encode(enc Encoder) (ok bool) {
	return true
}

// Decode decodes internal variables
func (o Rjoint) Decode(dec Decoder) (ok bool) {
	return true
}

// OutIpsData returns data from all integration points for output
func (o Rjoint) OutIpsData() (labels []string, data []*OutIpData) {
	return
}

// auxiliary ////////////////////////////////////////////////////////////////////////////////////////

// ipvars computes current values @ integration points. idx == index of integration point
func (o *Rjoint) ipvars(idx int, sol *Solution, delta bool) (ok bool) {

	// interpolation functions and gradients
	if LogErr(o.Rod.Cell.Shp.CalcAtIp(o.Rod.X, o.IpsElem[idx], true), "ipvars") {
		return
	}

	// auxiliary
	rodH := o.Rod.Cell.Shp
	sldH := o.Sld.Cell.Shp
	rodS := rodH.S
	rodNv := rodH.Nverts
	sldNv := sldH.Nverts
	nsig := 2 * o.Ndim

	// displacements or increments of displacements
	U := sol.Y
	if delta {
		U = sol.ΔY
	}

	// relative displacement
	var rlin, rsld int
	for i := 0; i < o.Ndim; i++ {
		o.Δw[i] = 0
	}
	for m := 0; m < rodNv; m++ {
		for i := 0; i < o.Ndim; i++ {
			for n := 0; n < sldNv; n++ {
				rsld = o.Sld.Umap[i+n*o.Ndim]
				o.Δw[i] += rodS[m] * o.sldS_rodN[m][n] * U[rsld]
			}
			rlin = o.Rod.Umap[i+m*o.Ndim]
			o.Δw[i] -= rodS[m] * U[rlin]
		}
	}
	o.Δwb0, o.Δwb1, o.Δwb2 = 0, 0, 0
	for i := 0; i < o.Ndim; i++ {
		o.Δwb0 += o.e0[idx][i] * o.Δw[i]
		o.Δwb1 += o.e1[idx][i] * o.Δw[i]
		o.Δwb2 += o.e2[idx][i] * o.Δw[i]
	}

	// state variables
	o.τ = o.States[idx].Sig
	o.qn1 = o.States[idx].Phi[0]
	o.qn2 = o.States[idx].Phi[1]

	// new confining stress
	o.σc = 0.0
	if o.Coulomb {

		// calculate σIp
		for j := 0; j < nsig; j++ {
			o.sldSigN_rodP[j] = 0
			for n := 0; n < sldNv; n++ {
				o.sldSigN_rodP[j] += o.sldS_rodP[idx][n] * o.sldSigN[n][j]
			}
		}

		// calculate t1 and t2
		for i := 0; i < o.Ndim; i++ {
			o.t1[i], o.t2[i] = 0, 0
			for j := 0; j < o.Ndim; j++ {
				o.t1[i] += tsr.M2T(o.sldSigN_rodP, i, j) * o.e1[idx][j]
				o.t2[i] += tsr.M2T(o.sldSigN_rodP, i, j) * o.e2[idx][j]
			}
		}

		// calculate p1, p2 and σcNew
		p1, p2 := 0.0, 0.0
		for i := 0; i < o.Ndim; i++ {
			p1 += o.t1[i] * o.e1[idx][i]
			p2 += o.t2[i] * o.e2[idx][i]
		}

		// σcNew
		o.σc = -(p1 + p2) / 2.0
	}
	return true
}
