// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"
	"os"
	"path/filepath"
	"sort"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

/* FaceCond holds information of one single face boundary condition. Example:
 *
 *                     -12 => "qn", "seepH"
 *  36    -12     35   -11 => "ux", "seepH"
 *   (3)--------(2)
 *    |    2     |     face id => conditions
 *    |          |           0 => <nil>
 *    |3        1| -11       1 => {"ux", "seepH"} => localVerts={1,2} => globalVerts={34,35}
 *    |          |           2 => {"qn", "seepH"} => localVerts={2,3} => globalVerts={35,36}
 *    |    0     |           3 => <nil>
 *   (0)--------(1)
 *  33            34    "seepH" => localVerts={1,2,3}
 *                      "seepH" => globalVerts={34,35,36}
 *
 *  face id => condition
 *        1 => "ux"    => localVerts={1,2} => globalVerts={34,35}
 *        1 => "seepH" => localVerts={1,2} => globalVerts={34,35}
 *        2 => "qn"    => localVerts={2,3} => globalVerts={35,36}
 *        2 => "seepH" => localVerts={2,3} => globalVerts={35,36}
 */

type FaceCond struct {
	FaceId      int      // msh: cell's face local id
	LocalVerts  []int    // msh: cell's face local vertices ids (sorted)
	GlobalVerts []int    // msh: global vertices ids (sorted)
	Cond        string   // sim: condition; e.g. "qn" or "seepH"
	Func        fun.Func // sim: function to compute boundary condition
	Extra       string   // sim: extra information
}

// GetVertsWithCond gets all vertices with any of the given conditions
//  "seepH" => localVerts={1,2,3}
func GetVertsWithCond(fconds []*FaceCond, conds ...string) (verts []int) {
	for _, fc := range fconds {
		// check if fc has any of conds
		if utl.StrIndexSmall(conds, fc.Cond) < 0 {
			continue
		}

		// add local verts
		for _, lv := range fc.LocalVerts {
			// check if vert was added already
			if utl.IntIndexSmall(verts, lv) < 0 {
				verts = append(verts, lv) // add a vert
			}
		}
	}

	sort.Ints(verts)
	return
}

// Solution holds the solution data @ nodes.
//        / u \         / u \
//        |   | => y =  |   |
//  yb =  | p |         \ p / (ny x 1)
//        |   |
//        \ λ / (nyb x 1)
//
type Solution struct {

	// state
	T      float64   // current time
	Y      []float64 // DOFs (solution variables); e.g. y = {u, p}
	Dydt   []float64 // dy/dt
	D2ydt2 []float64 // d²y/dt²

	// auxiliary
	ΔY  []float64 // total increment (for nonlinear solver)
	Psi []float64 // t1 star vars; e.g. ψ* = β1.p + β2.dpdt
	Zet []float64 // t2 star vars; e.g. ζ* = α1.u + α2.v + α3.a
	Chi []float64 // t2 star vars; e.g. χ* = α4.u + α5.v + α6.a
	L   []float64 // Lagrange multipliers
}

// Domain holds all Nodes and Elements active during a stage in addition to the Solution at nodes.
// Only elements in this processor are recorded here; however information from
// all cells might be recorded as well.
type Domain struct {

	// init: region, mesh, linear solver
	Reg    *inp.Region // region data
	Msh    *inp.Mesh   // mesh data
	LinSol la.LinSol   // linear solver

	// stage: auxiliary maps for setting boundary conditions
	FaceConds map[int][]*FaceCond // maps cell id to its face boundary conditions

	// stage: nodes (active) and elements (active AND in this processor)
	Nodes  []*Node // active nodes (for each stage)
	Elems  []Elem  // [procNcells] only active elements in this processor (for each stage)
	MyCids []int   // [procNcells] the ids of cells in this processor

	// stage: auxiliary maps for dofs and equation types
	F2Y      map[string]string // converts f-keys to y-keys; e.g.: "ux" => "fx"
	YandC    map[string]bool   // y and constraints keys; e.g. "ux", "pl", "H", "incsup", "rigid"
	Dof2Tnum map[string]int    // {t1,t2}-types: dof => t_number; e.g. "ux" => 2, "pl" => 1

	// stage: auxiliary maps for nodes and elements
	Vid2node   []*Node // [nverts] VertexId => index in Nodes. Inactive vertices are 'nil'
	Cid2elem   []Elem  // [ncells] CellId => index in Elems. Cells in other processors or inactive are 'nil'
	Cid2active []bool  // [ncells] CellId => whether cell is active or not in ANY processor

	// stage: subsets of elements
	ElemIntvars []ElemIntvars   // elements with internal vars in this processor
	ElemConnect []ElemConnector // connector elements in this processor

	// stage: coefficients and prescribed forces
	EssenBcs EssentialBcs // constraints (Lagrange multipliers)
	PtNatBcs PtNaturalBcs // point loads such as prescribed forces at nodes

	// stage: t1 and t2 variables
	T1eqs []int // first t-derivative variables; e.g.:  dp/dt vars (subset of ykeys)
	T2eqs []int // second t-derivative variables; e.g.: d²u/dt² vars (subset of ykeys)

	// stage: dimensions
	NnzKb int // number of nonzeros in Kb matrix
	Ny    int // total number of dofs, except λ
	Nlam  int // total number of Lagrange multipliers
	NnzA  int // number of nonzeros in A (constraints) matrix
	Nyb   int // total number of equations: ny + nλ

	// stage: solution and linear solver
	Sol      *Solution   // solution state
	Kb       *la.Triplet // Jacobian == dRdy
	Fb       []float64   // residual == -fb
	Wb       []float64   // workspace
	InitLSol bool        // flag telling that linear solver needs to be initialised prior to any further call
}

// NewDomain returns a new domain
func NewDomain(reg *inp.Region, distr bool) *Domain {
	var dom Domain
	dom.Reg = reg
	dom.Msh = reg.Msh
	if distr {
		if LogErrCond(Global.Nproc != len(dom.Msh.Part2cells), "number of processors must be equal to the number of partitions defined in mesh file. %d != %d", Global.Nproc, len(dom.Msh.Part2cells)) {
			return nil
		}
	}
	dom.LinSol = la.GetSolver(Global.Sim.LinSol.Name)
	return &dom
}

// SetStage set nodes, equation numbers and auxiliary data for given stage
func (o *Domain) SetStage(idxstg int, stg *inp.Stage, distr bool) (setstageisok bool) {

	// backup state
	if idxstg > 0 {
		o.create_stage_copy()
		if !o.fix_inact_flags(stg.Activate, false) {
			return
		}
		if !o.fix_inact_flags(stg.Deactivate, true) {
		}
	}

	// auxiliary maps for setting boundary conditions
	o.FaceConds = make(map[int][]*FaceCond) // cid => conditions

	// nodes (active) and elements (active AND in this processor)
	o.Nodes = make([]*Node, 0)
	o.Elems = make([]Elem, 0)
	o.MyCids = make([]int, 0)

	// auxiliary maps for dofs and equation types
	o.F2Y = make(map[string]string)
	o.YandC = GetIsEssenKeyMap()
	o.Dof2Tnum = make(map[string]int)

	// auxiliary maps for nodes and elements
	o.Vid2node = make([]*Node, len(o.Msh.Verts))
	o.Cid2elem = make([]Elem, len(o.Msh.Cells))
	o.Cid2active = make([]bool, len(o.Msh.Cells))

	// subsets of elements
	o.ElemConnect = make([]ElemConnector, 0)
	o.ElemIntvars = make([]ElemIntvars, 0)

	// allocate nodes and cells (active only) -------------------------------------------------------

	// for each cell
	var eq int // current equation number => total number of equations @ end of loop
	o.NnzKb = 0
	for _, c := range o.Msh.Cells {

		// get element data and information structure
		edat := o.Reg.Etag2data(c.Tag)
		if LogErrCond(edat == nil, "cannot get element's data with etag=%d", c.Tag) {
			return
		}
		if edat.Inact {
			continue
		}
		o.Cid2active[c.Id] = true

		// prepare maps of face conditions
		for faceId, faceTag := range c.FTags {
			if faceTag < 0 {
				faceBc := stg.GetFaceBc(faceTag)
				if faceBc != nil {
					lverts := shp.GetFaceLocalVerts(c.Type, faceId)
					gverts := o.faceLocal2globalVerts(lverts, c)
					for j, key := range faceBc.Keys {
						fcn := Global.Sim.Functions.Get(faceBc.Funcs[j])
						fcond := &FaceCond{faceId, lverts, gverts, key, fcn, faceBc.Extra}
						fconds := o.FaceConds[c.Id]
						o.FaceConds[c.Id] = append(fconds, fcond)
					}
				}
			}
		}

		// get element info (such as DOFs, etc.)
		info := GetElemInfo(c.Type, edat.Type, o.FaceConds[c.Id])
		if info == nil {
			return
		}

		// for non-joint elements, add new DOFs
		if !c.IsJoint {
			chk.IntAssert(len(info.Dofs), len(c.Verts))

			// store y and f information
			for ykey, fkey := range info.Y2F {
				o.F2Y[fkey] = ykey
				o.YandC[ykey] = true
			}

			// t1 and t2 equations
			for _, ykey := range info.T1vars {
				o.Dof2Tnum[ykey] = 1
			}
			for _, ykey := range info.T2vars {
				o.Dof2Tnum[ykey] = 2
			}

			// loop over nodes of this element
			var eNdof int // number of DOFs of this elmeent
			for j, v := range c.Verts {

				// new or existent node
				var nod *Node
				if o.Vid2node[v] == nil {
					nod = NewNode(o.Msh.Verts[v])
					o.Vid2node[v] = nod
					o.Nodes = append(o.Nodes, nod)
				} else {
					nod = o.Vid2node[v]
				}

				// set DOFs and equation numbers
				for _, ukey := range info.Dofs[j] {
					eq = nod.AddDofAndEq(ukey, eq)
					eNdof += 1
				}
			}

			// number of non-zeros
			o.NnzKb += eNdof * eNdof
		}

		// allocate element
		mycell := c.Part == Global.Rank // cell belongs to this processor
		if !distr {
			mycell = true // not distributed simulation => this processor has all cells
		}
		if mycell {

			// new element
			ele := NewElem(edat, c.Id, o.Msh, o.FaceConds[c.Id])
			if ele == nil {
				return
			}
			o.Cid2elem[c.Id] = ele
			o.Elems = append(o.Elems, ele)
			o.MyCids = append(o.MyCids, ele.Id())

			// give equation numbers to new element
			eqs := make([][]int, len(c.Verts))
			for j, v := range c.Verts {
				for _, dof := range o.Vid2node[v].Dofs {
					eqs[j] = append(eqs[j], dof.Eq)
				}
			}
			ele.SetEqs(eqs, nil)

			// subsets of elements
			o.add_element_to_subsets(ele)
		}
	}

	// connect elements (e.g. Joints)
	for _, e := range o.ElemConnect {
		nnz, ok := e.Connect(o.Cid2elem, o.Msh.Cells[e.Id()])
		if LogErrCond(!ok, "Connect failed") {
			return
		}
		o.NnzKb += nnz
	}

	// logging
	log.Printf("dom: stage # %d %s\n", idxstg, stg.Desc)
	log.Printf("dom: nnodes=%d nelems=%d\n", len(o.Nodes), len(o.Elems))

	// element conditions, essential and natural boundary conditions --------------------------------

	// (re)set constraints and prescribed forces structures
	o.EssenBcs.Reset()
	o.PtNatBcs.Reset()

	// element conditions
	for _, ec := range stg.EleConds {
		cells, ok := o.Msh.CellTag2cells[ec.Tag]
		if LogErrCond(!ok, "cannot find cells with tag = %d to assign conditions", ec.Tag) {
			return
		}
		for _, c := range cells {
			e := o.Cid2elem[c.Id]
			if e != nil { // set conditions only for this processor's / active element
				for j, key := range ec.Keys {
					fcn := Global.Sim.Functions.Get(ec.Funcs[j])
					if LogErrCond(fcn == nil, "Functions.Get failed\n") {
						return
					}
					e.SetEleConds(key, fcn, ec.Extra)
				}
			}
		}
	}

	// face boundary conditions
	for cidx, fcs := range o.FaceConds {
		c := o.Msh.Cells[cidx]
		for _, fc := range fcs {
			lverts := fc.LocalVerts
			gverts := o.faceLocal2globalVerts(lverts, c)
			var enodes []*Node
			for _, v := range gverts {
				enodes = append(enodes, o.Vid2node[v])
			}
			if o.YandC[fc.Cond] {
				if !o.EssenBcs.Set(fc.Cond, enodes, fc.Func, fc.Extra) {
					return
				}
			}

		}
	}

	// vertex bounday conditions
	for _, nc := range stg.NodeBcs {
		verts, ok := o.Msh.VertTag2verts[nc.Tag]
		if LogErrCond(!ok, "cannot find vertices with tag = %d to assign node boundary conditions", nc.Tag) {
			return
		}
		for _, v := range verts {
			if o.Vid2node[v.Id] != nil { // set BCs only for active nodes
				n := o.Vid2node[v.Id]
				for j, key := range nc.Keys {
					fcn := Global.Sim.Functions.Get(nc.Funcs[j])
					if LogErrCond(fcn == nil, "Functions.Get failed\n") {
						return
					}
					if o.YandC[key] {
						o.EssenBcs.Set(key, []*Node{n}, fcn, nc.Extra)
					} else {
						o.PtNatBcs.Set(o.F2Y[key], n, fcn, nc.Extra)
					}
				}
			}
		}
	}

	// resize slices --------------------------------------------------------------------------------

	// t1 and t2 equations
	o.T1eqs = make([]int, 0)
	o.T2eqs = make([]int, 0)
	for _, nod := range o.Nodes {
		for _, dof := range nod.Dofs {
			switch o.Dof2Tnum[dof.Key] {
			case 1:
				o.T1eqs = append(o.T1eqs, dof.Eq)
			case 2:
				o.T2eqs = append(o.T2eqs, dof.Eq)
			default:
				LogErrCond(true, "t1 and t2 equations are incorrectly set")
				return
			}
		}
	}

	// size of arrays
	o.Ny = eq
	o.Nlam, o.NnzA = o.EssenBcs.Build(o.Ny)
	o.Nyb = o.Ny + o.Nlam

	// solution structure and linear solver
	o.Sol = new(Solution)
	o.Kb = new(la.Triplet)
	o.Fb = make([]float64, o.Nyb)
	o.Wb = make([]float64, o.Nyb)
	o.Kb.Init(o.Nyb, o.Nyb, o.NnzKb+2*o.NnzA)
	o.InitLSol = true // tell solver that lis has to be initialised before use

	// allocate arrays
	o.Sol.Y = make([]float64, o.Ny)
	o.Sol.ΔY = make([]float64, o.Ny)
	o.Sol.L = make([]float64, o.Nlam)
	if !Global.Sim.Data.Steady {
		o.Sol.Dydt = make([]float64, o.Ny)
		o.Sol.D2ydt2 = make([]float64, o.Ny)
		o.Sol.Psi = make([]float64, o.Ny)
		o.Sol.Zet = make([]float64, o.Ny)
		o.Sol.Chi = make([]float64, o.Ny)
	}

	// initialise internal variables
	if stg.HydroSt {
		if !o.SetHydroSt(stg) {
			return
		}
	} else if stg.GeoSt != nil {
		if !o.SetGeoSt(stg) {
			return
		}
	} else if stg.IniStress != nil {
		if !o.SetIniStress(stg) {
			return
		}
	} else {
		for _, e := range o.ElemIntvars {
			e.SetIniIvs(o.Sol, nil)
		}
	}

	// import results from another set of files
	if stg.Import != "" {
		fnp := os.ExpandEnv(stg.Import)
		dir := filepath.Dir(fnp)
		fnk := io.FnKey(filepath.Base(fnp))
		sum := ReadSum(dir, fnk)
		if LogErrCond(sum == nil, "cannot import state from %s", stg.Import) {
			return
		}
		if !o.In(sum, len(sum.OutTimes)-1, false) {
			return
		}
	}

	// logging
	if Global.LogBcs {
		log.Printf("dom: essential boundary conditions:%v", o.EssenBcs.List(stg.Control.Tf))
		log.Printf("dom: ptnatbcs=%v", o.PtNatBcs.List(stg.Control.Tf))
	}
	log.Printf("dom: ny=%d nlam=%d nnzKb=%d nnzA=%d nt1eqs=%d nt2eqs=%d", o.Ny, o.Nlam, o.NnzKb, o.NnzA, len(o.T1eqs), len(o.T2eqs))

	// success
	return true
}

// auxiliary functions //////////////////////////////////////////////////////////////////////////////

// add_element_to_subsets adds an Elem to many subsets as it fits
func (o *Domain) add_element_to_subsets(ele Elem) {
	if e, ok := ele.(ElemIntvars); ok {
		o.ElemIntvars = append(o.ElemIntvars, e)
	}
	if e, ok := ele.(ElemConnector); ok {
		o.ElemConnect = append(o.ElemConnect, e)
	}
}

// star_vars computes starred variables
func (o *Domain) star_vars(Δt float64) (err error) {

	// skip if steady simulation
	if Global.Sim.Data.Steady {
		return
	}

	// recompute coefficients
	dc := Global.DynCoefs
	err = dc.CalcBoth(Δt)
	if err != nil {
		return
	}

	// compute starred vectors
	for _, I := range o.T1eqs {
		o.Sol.Psi[I] = dc.β1*o.Sol.Y[I] + dc.β2*o.Sol.Dydt[I]
	}
	for _, I := range o.T2eqs {
		o.Sol.Zet[I] = dc.α1*o.Sol.Y[I] + dc.α2*o.Sol.Dydt[I] + dc.α3*o.Sol.D2ydt2[I]
		o.Sol.Chi[I] = dc.α4*o.Sol.Y[I] + dc.α5*o.Sol.Dydt[I] + dc.α6*o.Sol.D2ydt2[I]
	}

	// set internal starred variables
	for _, e := range o.Elems {
		e.InterpStarVars(o.Sol)
	}
	return
}

// create_stage_copy creates a copy of current stage => to be used later when activating/deactivating elements
func (o *Domain) create_stage_copy() {
}

// set_act_deact_flags sets inactive flags for new active/inactive elements
func (o *Domain) fix_inact_flags(eids_or_tags []int, deactivate bool) (ok bool) {
	for _, tag := range eids_or_tags {
		if tag >= 0 { // this meahs that tag == cell.Id
			cell := o.Msh.Cells[tag]
			tag = cell.Tag
		}
		edat := o.Reg.Etag2data(tag)
		if LogErrCond(edat == nil, "cannot get element's data with etag=%d", tag) {
			return
		}
		edat.Inact = deactivate
	}
	return true
}

// faceLocal2globalVerts returns the face global vertices ids
func (o Domain) faceLocal2globalVerts(faceLverts []int, cell *inp.Cell) (faceGverts []int) {
	faceGverts = make([]int, len(faceLverts))
	for i, l := range faceLverts {
		faceGverts[i] = cell.Verts[l]
	}
	return
}
