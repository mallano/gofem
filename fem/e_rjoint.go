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

// Rjoint implements the rod-joint (interface/link) element for reinforced solids
type Rjoint struct {

	// basic data
	Edat *inp.ElemData // element data; stored in allocator to be used in Connect
	Cell *inp.Cell     // cell
	Ndim int           // space dimension
	Ny   int           // total number of dofs == Nu_rod + Nu_solid

	// supporting elements
	Rod *Rod   // rod element
	Sld *ElemU // solid element

	// integration points
	IpsElem []*shp.Ipoint // [nip] integration points of element

	// material model and internal variables
	Model msolid.RjointM1 // material model

	// parameters
	h  float64 // perimeter of rod element
	k1 float64 // lateral stiffness
	k2 float64 // lateral stiffness

	// shape functions evaluations. nv == nverts
	SslNo [][]float64 // [rodNv][sldNv] shape functions of solids @ nodes of line element
	SslIp [][]float64 // [rodNip][sldNv] shape functions of solids @ ips of line element

	// corotational system aligned with rod element
	e0 [][]float64 // [nipRod][ndim] local directions at each ip of rod
	e1 [][]float64 // [nipRod][ndim] local directions at each ip of rod
	e2 [][]float64 // [nipRod][ndim] local directions at each ip of rod

	// variables for Coulomb model
	Coulomb bool            // uses Coulomb model
	rsIp    [][]float64     // [rodNip][ndim] natural coordinates w.r.t the solids' system of ips of rod
	σNo     [][]float64     // [sldNv][nsig] σ at nodes of solid
	σIp     []float64       // [nsig] σ at ips of rod
	DσNoDu  [][][][]float64 // [sldNv][nσ][sldNv][ndim] ∂σSldNod/∂uSldNod : derivatives of σ @ nodes of solid w.r.t displacements of solid
	DσDun   [][]float64     // [nσ][ndim] ∂σSldIp/∂uSldNod : derivatives of σ @ ip of solid w.r.t displacements of solid
	t1      []float64       // traction vectors for σc
	t2      []float64       // traction vectors for σc
	T1      [][]float64     // [nipRod][nσ] tensor (e1 dy e1)
	T2      [][]float64     // [nipRod][nσ] tensor (e2 dy e2)
	sldE    [][]float64     // [sldNv][sldNip]

	// problem variables
	Ymap []int // assembly map (location array/element equations)

	// internal variables
	States    []*msolid.OnedState // [nip] internal states
	StatesBkp []*msolid.OnedState // [nip] backup internal states

	// scratchpad. computed @ each ip
	Δw   []float64 // relative velocity
	qb   []float64 // resultant traction vector 'holding' the rod @ ip
	σc   float64
	Δwb0 float64
	Δwb1 float64
	Δwb2 float64
	τ    float64
	qn1  float64
	qn2  float64
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
	o.Model.Init(matdata.Prms)

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
	rodNip := len(o.Rod.IpsElem)
	sldNip := len(o.Sld.IpsElem)
	rodH := o.Rod.Cell.Shp
	sldH := o.Sld.Cell.Shp
	rodS := rodH.S
	sldS := sldH.S
	rodNv := rodH.Nverts
	sldNv := sldH.Nverts
	nsig := 2 * o.Ndim

	// shape functions of solid @ nodes of rod
	o.SslNo = la.MatAlloc(rodNv, sldNv)
	yno := make([]float64, o.Ndim) // real coordinate of rod node
	rls := make([]float64, 3)      // natural coordinates of node of rod w.r.t solid's natural system
	for m := 0; m < rodNv; m++ {
		for i := 0; i < o.Ndim; i++ {
			yno[i] = o.Rod.X[i][m]
		}
		sldH.InvMap(rls, yno, o.Sld.X)
		sldH.CalcAtR(o.Sld.X, rls, false)
		for n := 0; n < sldNv; n++ {
			o.SslNo[m][n] = sldH.S[n]
		}
	}

	// coulomb model => σc depends on p values of solid
	if o.Coulomb {

		// variables for σc
		o.σNo = la.MatAlloc(sldNv, nsig)
		o.σIp = make([]float64, nsig)
		o.t1 = make([]float64, o.Ndim)
		o.t2 = make([]float64, o.Ndim)

		// extrapolator matrix
		o.sldE = la.MatAlloc(sldNv, sldNip)
		if LogErr(sldH.Extrapolator(o.sldE, o.Sld.IpsElem), "Extrapolator of solid failed") {
			return
		}

		// shape function of solid @ ips of rod
		o.SslIp = la.MatAlloc(rodNip, sldNv)
		o.rsIp = la.MatAlloc(rodNip, 3)
		yip := make([]float64, o.Ndim) // real coordiantes of ip of rod
		for idx, ip := range o.Rod.IpsElem {
			rodH.CalcAtIp(o.Rod.X, ip, false)
			for i := 0; i < o.Ndim; i++ {
				yip[i] = 0
				for m := 0; m < rodNv; m++ {
					yip[i] += rodS[m] * o.Rod.X[i][m]
				}
			}
			sldH.InvMap(o.rsIp[idx], yip, o.Sld.X)
			sldH.CalcAtR(o.Sld.X, o.rsIp[idx], false)
			for n := 0; n < sldNv; n++ {
				o.SslIp[idx][n] = sldS[n]
			}
		}

		// allocate derivatives of p @ node and ips of solid w.r.t diplacements of solid
		o.DσNoDu = utl.Deep4alloc(sldNv, nsig, sldNv, o.Ndim)
		o.DσDun = la.MatAlloc(nsig, o.Ndim)
		o.T1 = la.MatAlloc(rodNip, nsig)
		o.T2 = la.MatAlloc(rodNip, nsig)
	}

	// joint direction @ ip[idx]; corotational system aligned with rod element
	o.e0 = la.MatAlloc(rodNip, o.Ndim)
	o.e1 = la.MatAlloc(rodNip, o.Ndim)
	o.e2 = la.MatAlloc(rodNip, o.Ndim)
	a := make([]float64, o.Ndim)
	Q := la.MatAlloc(o.Ndim, o.Ndim)
	for idx, ip := range o.Rod.IpsElem {
		if LogErr(rodH.CalcAtIp(o.Rod.X, ip, true), "CalcAtIp of rod failed") {
			return
		}
		la.VecCopy(a, 1.0, rodH.Jvec3d[:o.Ndim])
		a[0] += 666.0
		la.VecCopy(o.e0[idx], 1.0/rodH.J, rodH.Jvec3d[:o.Ndim]) // e0 := Jvec / J
		la.MatSetDiag(Q, 1.0)
		la.VecOuterAdd(Q, -1.0, o.e0[idx], o.e0[idx])
		la.MatVecMul(o.e1[idx], 1.0, Q, a)
		la.VecScale(o.e1[idx], 0, 1.0/la.VecNorm(o.e1[idx]), o.e1[idx])
		if o.Ndim == 3 {
			o.e2[idx][0] = o.e0[idx][1]*o.e1[idx][2] - o.e0[idx][2]*o.e1[idx][1]
			o.e2[idx][1] = o.e0[idx][2]*o.e1[idx][0] - o.e0[idx][0]*o.e1[idx][2]
			o.e2[idx][2] = o.e0[idx][0]*o.e1[idx][1] - o.e0[idx][1]*o.e1[idx][0]
		}
		if o.Coulomb {
			e1_dy_e1 := tsr.Alloc2()
			e2_dy_e2 := tsr.Alloc2()
			for i := 0; i < o.Ndim; i++ {
				for j := 0; j < o.Ndim; j++ {
					e1_dy_e1[i][j] = o.e1[idx][i] * o.e1[idx][j]
					e2_dy_e2[i][j] = o.e2[idx][i] * o.e2[idx][j]
				}
			}
			tsr.Ten2Man(o.T1[idx], e1_dy_e1)
			tsr.Ten2Man(o.T2[idx], e2_dy_e2)
		}
	}

	// full Usmap
	o.Ymap = make([]int, o.Ny)
	for m := 0; m < sldNv; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
			o.Ymap[r] = o.Sld.Umap[r]
		}
	}
	start := o.Sld.Nu
	for m := 0; m < rodNv; m++ {
		for i := 0; i < o.Ndim; i++ {
			r := i + m*o.Ndim
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
				fb[rsld] += o.SslNo[m][n] //* o.Rus[rlin] // -fi  :  f = fC
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
				o.Δw[i] += rodS[m] * o.SslNo[m][n] * U[rsld]
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
			o.σIp[j] = 0
			for n := 0; n < sldNv; n++ {
				o.σIp[j] += o.SslIp[idx][n] * o.σNo[n][j]
			}
		}

		// calculate t1 and t2
		for i := 0; i < o.Ndim; i++ {
			o.t1[i], o.t2[i] = 0, 0
			for j := 0; j < o.Ndim; j++ {
				o.t1[i] += tsr.M2T(o.σIp, i, j) * o.e1[idx][j]
				o.t2[i] += tsr.M2T(o.σIp, i, j) * o.e2[idx][j]
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
