// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// OutIpData is an auxiliary structure to transfer data from integration points (IP) to output routines.
type OutIpData struct {
	Eid int                 // id of element that owns this ip
	X   []float64           // coordinates
	V   map[string]*float64 // maps label (e.g. "sx") to pointer to value
}

// Elem defines what elements must calculate
type Elem interface {

	// information and initialisation
	Id() int                                           // returns the cell Id
	SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) // set equations

	// conditions (natural BCs and element's)
	SetEleConds(key string, f fun.Func, extra string) (ok bool) // set element conditions

	// called for each time step
	InterpStarVars(sol *Solution) (ok bool) // interpolate star variables to integration points

	// called for each iteration
	AddToRhs(fb []float64, sol *Solution) (ok bool)                // adds -R to global residual vector fb
	AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) // adds element K to global Jacobian matrix Kb
	Update(sol *Solution) (ok bool)                                // perform (tangent) update

	// reading and writing of element data
	Encode(enc Encoder) (ok bool) // encodes internal variables
	Decode(dec Decoder) (ok bool) // decodes internal variables

	// output
	OutIpsData() (data []*OutIpData) // returns data from all integration points for output
}

// ElemConnector defines connector elements; elements that depend upon others
type ElemConnector interface {
	Id() int                                                  // returns the cell Id
	Connect(cid2elem []Elem, c *inp.Cell) (nnzK int, ok bool) // connect multiple elements; e.g.: connect rod/solid elements in Rjoints
}

// ElemIntvars defines elements with {z,q} internal variables
type ElemIntvars interface {
	InitIvs(sol *Solution) (ok bool)           // reset (and fix) internal variables after primary variables have been changed
	SetIvs(ivs map[string][]float64) (ok bool) // set secondary variables; e.g. during initialisation via files
	BackupIvs() (ok bool)                      // create copy of internal variables
	RestoreIvs() (ok bool)                     // restore internal variables from copies
}

// Info holds all information required to set a simulation stage
type Info struct {

	// essential
	Dofs [][]string        // solution variables PER NODE. ex for 2 nodes: [["ux", "uy", "rz"], ["ux", "uy", "rz"]]
	Y2F  map[string]string // maps "y" keys to "f" keys. ex: "ux" => "fx", "pl" => "ql"

	// internal Dofs; e.g. for mixed formulations
	NintDofs int // number of internal dofs

	// t1 and t2 variables (time-derivatives of first and second order)
	T1vars []string // "pl"
	T2vars []string // "ux", "uy"
}

// GetElemInfo returns information about elements/formulations
//  cellType -- e.g. "qua8"
//  elemType -- e.g. "u"
func GetElemInfo(ndim int, cellType, elemType string, faceConds []*FaceCond) *Info {
	infogetter, ok := infogetters[elemType]
	if LogErrCond(!ok, "cannot find element type = %s", elemType) {
		return nil
	}
	info := infogetter(ndim, cellType, faceConds)
	if LogErrCond(info == nil, "cannot find info from %q element", elemType) {
		return nil
	}
	return info
}

// NewElem returns a new element from its type; e.g. "p", "u" or "up"
func NewElem(edat *inp.ElemData, cid int, msh *inp.Mesh, faceConds []*FaceCond) Elem {
	elemType := edat.Type
	allocator, ok := eallocators[elemType]
	if LogErrCond(!ok, "cannot find element type = %s", elemType) {
		return nil
	}
	c := msh.Cells[cid]
	x := BuildCoordsMatrix(c, msh)
	ele := allocator(msh.Ndim, c.Type, faceConds, cid, edat, x)
	if LogErrCond(ele == nil, "cannot allocate %q element", elemType) {
		return nil
	}
	return ele
}

// BuildCoordsMatrix returns the coordinate matrix of a particular Cell
func BuildCoordsMatrix(c *inp.Cell, msh *inp.Mesh) (x [][]float64) {
	x = la.MatAlloc(msh.Ndim, len(c.Verts))
	for i := 0; i < msh.Ndim; i++ {
		for j, v := range c.Verts {
			x[i][j] = msh.Verts[v].C[i]
		}
	}
	return
}

// infogetters holds all available formulations/info; elemType => infogetter
var infogetters = make(map[string]func(ndim int, cellType string, faceConds []*FaceCond) *Info)

// eallocators holds all available elements; elemType => eallocator
var eallocators = make(map[string]func(ndim int, cellType string, faceConds []*FaceCond, cid int, edat *inp.ElemData, x [][]float64) Elem)
