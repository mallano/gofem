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
	SetEqs(eqs [][]int, mixedform_eqs []int) // set equations

	// conditions (natural BCs and element's)
	SetEleConds(key string, f fun.Func, extra string)               // set element conditions
	SetSurfLoads(key string, idxface int, f fun.Func, extra string) // set surface loads (natural boundary conditions)

	// called for each time step
	InterpStarVars(sol *Solution) error // interpolate star variables to integration points

	// called for each iteration
	AddToRhs(fb []float64, sol *Solution) error                // adds -R to global residual vector fb
	AddToKb(Kb *la.Triplet, sol *Solution, firstIt bool) error // adds element K to global Jacobian matrix Kb
	Update(sol *Solution) error                                // perform (tangent) update
}

// ElemConnector defines connector elements; elements that depend upon others
type ElemConnector interface {
	Connect(elems []Elem, cid2elem []Elem) // connect multiple elements; e.g.: connect rod/solid elements in Rjoints
}

// ElemWriter defines elements that can write output
type ElemWriter interface {
	Encode(enc Encoder) (err error) // encodes internal variables
	Decode(dec Decoder) (err error) // decodes internal variables
	//Write(enc Encoder, sol *Solution) // writes internal variables to Writer
}

// ElemIntvars defines elements with {z,q} internal variables
type ElemIntvars interface {
	InitIvs(sol *Solution)             // reset (and fix) internal variables after primary variables have been changed
	SetIvs(zvars map[string][]float64) // set secondary variables; e.g. during initialisation via files
	BackupIvs()                        // create copy of internal variables
	RestoreIvs()                       // restore internal variables from copies
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
func GetElemInfo(edat *inp.ElemData, cid int, msh *inp.Mesh) Info {
	allocator, ok := iallocators[edat.Type]
	PanicOrNot(!ok, "cannot find element type = %s", edat.Type)
	return allocator(edat, cid, msh)
}

// NewElem returns a new element from its type; e.g. "p", "u" or "up"
func NewElem(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem {
	allocator, ok := eallocators[edat.Type]
	PanicOrNot(!ok, "cannot find element type = %s", edat.Type)
	return allocator(edat, cid, msh)
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
var iallocators = make(map[string]func(edat *inp.ElemData, cid int, msh *inp.Mesh) Info)

// eallocators holds all available elements; elemType => eallocator
var eallocators = make(map[string]func(edat *inp.ElemData, cid int, msh *inp.Mesh) Elem)
