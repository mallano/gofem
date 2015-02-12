// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
)

// Elem defines what elements must calculate
type Elem interface {

	// initialisation
	SetEqs(eqs [][]int, mixedform_eqs []int) (ok bool) // set equations

	// conditions (natural BCs and element's)
	SetEleConds(key string, f fun.Func, extra string) (ok bool)            // set element conditions
	SetNatBcs(key string, idxface int, f fun.Func, extra string) (ok bool) // set natural boundary conditions

	// called for each time step
	InterpStarVars(sol *Solution) (ok bool) // interpolate star variables to integration points

	// called for each iteration
	AddToRhs(fb []float64, sol *Solution) (ok bool)                // adds -R to global residual vector fb
	AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) (ok bool) // adds element K to global Jacobian matrix Kb
	Update(sol *Solution) (ok bool)                                // perform (tangent) update
}

// ElemConnector defines connector elements; elements that depend upon others
type ElemConnector interface {
	Connect(elems []Elem, cid2elem []Elem) (ok bool) // connect multiple elements; e.g.: connect rod/solid elements in Rjoints
}

// ElemWriter defines elements that can write output
type ElemWriter interface {
	Encode(enc Encoder) (ok bool) // encodes internal variables
	Decode(dec Decoder) (ok bool) // decodes internal variables
	//Write(enc Encoder, sol *Solution) // writes internal variables to Writer
}

// ElemIntvars defines elements with {z,q} internal variables
type ElemIntvars interface {
	InitIvs(sol *Solution) (ok bool)             // reset (and fix) internal variables after primary variables have been changed
	SetIvs(zvars map[string][]float64) (ok bool) // set secondary variables; e.g. during initialisation via files
	BackupIvs() (ok bool)                        // create copy of internal variables
	RestoreIvs() (ok bool)                       // restore internal variables from copies
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
func GetElemInfo(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info {
	allocator, ok := iallocators[edat.Type]
	if LogErrCond(!ok, "cannot find element type = %s", edat.Type) {
		return nil
	}
	info := allocator(edat, cid, msh)
	if LogErrCond(info == nil, "cannot find info from %q element", edat.Type) {
		return nil
	}
	return info
}

// NewElem returns a new element from its type; e.g. "p", "u" or "up"
func NewElem(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {
	allocator, ok := eallocators[edat.Type]
	if LogErrCond(!ok, "cannot find element type = %s", edat.Type) {
		return nil
	}
	ele := allocator(edat, cid, msh)
	if LogErrCond(ele == nil, "cannot allocate %q element", edat.Type) {
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

// iallocators holds all available formulations/info; elemType => iallocator
var iallocators = make(map[string]func(edat *inp.ElemData, cid int, msh *inp.Mesh) *Info)

// eallocators holds all available elements; elemType => eallocator
var eallocators = make(map[string]func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem)
