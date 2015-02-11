// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/utl"
)

// EssentialBc holds information about essential bounday conditions such as constrained nodes.
// Lagrange multipliers are used to implement both single- and multi-point constraints.
//  In general, essential bcs / constraints are defined by means of:
//
//      A * y = c
//
//  The resulting Kb matrix will then have the following form:
//      _       _
//     |  K  At  | / δy \   / -R - At*λ \
//     |         | |    | = |           |
//     |_ A   0 _| \ δλ /   \  c - A*y  /
//         Kb       δyb          fb
//
type EssentialBc struct {
	Key   string    // ux, uy, rigid, incsup
	Eqs   []int     // equations
	ValsA []float64 // values for matrix A
	Fcn   fun.Func  // function that implements the "c" in A * y = c
	Inact bool      // inactive
}

// EssentialBcs implements a structure to record the definition of essential bcs / constraints.
// Each constraint will have a unique Lagrange multiplier index.
type EssentialBcs struct {
	Ndim   int            // space dimension
	Eq2idx map[int][]int  // maps eq number to indices in cstmp
	Bcs    []*EssentialBc // active essential bcs / constraints
	A      la.Triplet     // matrix of coefficients 'A'
	Am     *la.CCMatrix   // compressed form of A matrix

	// temporary
	BcsTmp []*EssentialBc // temporary essential bcs / constraints, including inactive ones
}

// Reset initialises this structure. It also performs a reset of internal structures.
func (o *EssentialBcs) Reset(ndim int) {
	o.Ndim = ndim
	o.BcsTmp = make([]*EssentialBc, 0)
	o.Eq2idx = make(map[int][]int)
	o.Bcs = make([]*EssentialBc, 0)
}

// Build builds this structure and its iternal data
//  nλ -- is the number of essential bcs / constraints == number of Lagrange multipliers
//  nnzA -- is the number of non-zeros in matrix 'A'
func (o *EssentialBcs) Build(ny int) (nλ, nnzA int) {

	// count number of active constraints and non-zeros in matrix A
	for _, c := range o.BcsTmp {
		if !c.Inact {
			o.Bcs = append(o.Bcs, c)
			nλ += 1
			nnzA += len(c.ValsA)
		}
	}

	// skip if there are no constraints
	if nλ == 0 {
		return
	}

	// set matrix A
	o.A.Init(nλ, ny, nnzA)
	for i, c := range o.Bcs {
		for j, eq := range c.Eqs {
			o.A.Put(i, eq, c.ValsA[j])
		}
	}
	o.Am = o.A.ToMatrix(nil)
	return
}

// AddtoRhs adds the essential bcs / constraints terms to the augmented fb vector
func (o EssentialBcs) AddToRhs(fb []float64, sol *Solution) {

	// skip if there are no constraints
	if len(o.Bcs) == 0 {
		return
	}

	// add -At*λ to fb
	la.SpMatTrVecMulAdd(fb, -1, o.Am, sol.L) // fb += -1 * At * λ

	// assemble -rc = c - A*y into fb
	ny := len(sol.Y)
	for i, c := range o.Bcs {
		fb[ny+i] = c.Fcn.F(sol.T, nil)
	}
	la.SpMatVecMulAdd(fb[ny:], -1, o.Am, sol.Y) // fb += -1 * A * y
}

// add adds new essential bcs / constraint and sets map eq2idx
func (o *EssentialBcs) add(key string, eqs []int, valsA []float64, fcn fun.Func) {
	idx := len(o.BcsTmp)
	o.BcsTmp = append(o.BcsTmp, &EssentialBc{key, eqs, valsA, fcn, false})
	for _, eq := range eqs {
		utl.IntIntsMapAppend(&o.Eq2idx, eq, idx)
	}
}

// add_single adds single-point constraint
func (o *EssentialBcs) add_single(key string, eq int, fcn fun.Func) {
	for _, idx := range o.Eq2idx[eq] {
		c := o.BcsTmp[idx]
		if c.Key == "rigid" || c.Key == "incsup" {
			return
		}
		c.Inact = true
	}
	o.add(key, []int{eq}, []float64{1}, fcn)
}

// GetFirstYandCmap returns the initial "yandc" map with additional keys that EssentialBcs can handle
func GetIsEssenKeyMap() map[string]bool {
	return map[string]bool{"rigid": true, "incsup": true, "H": true}
}

// Set sets a constraint if it does NOT exist yet.
//  key   -- can be Dof key such as "ux", "uy" or constraint type such as "mpc" or "rigid"
//  extra -- is a keycode-style data. e.g. "!type:incsup2d !alp:30"
//  Notes: 1) the default for key is single point constraint; e.g. "ux", "uy", ...
//         2) hydraulic head can be set with key == "H"
func (o *EssentialBcs) Set(key string, nodes []*Node, fcn fun.Func, extra string) (setisok bool) {

	// len(nod) must be greater than 0
	utl.IntAssertLessThan(0, len(nodes)) // 0 < len(nod)

	// rigid element
	if key == "rigid" {
		a := nodes[0].dofs
		for i := 1; i < len(nodes); i++ {
			for j, b := range nodes[i].dofs {
				o.add(key, []int{a[j].Eq, b.Eq}, []float64{1, -1}, &fun.Zero)
			}
		}
		return true // success
	}

	// inclined support
	if key == "incsup" {

		// check
		if LogErrCond(o.Ndim != 2, "inclined support works only in 2D for now") {
			return false // problem
		}

		// get data
		var α float64
		if val, found := utl.Keycode(extra, "alp"); found {
			α = utl.Atof(val) * math.Pi / 180.0
		}
		co, si := math.Cos(α), math.Sin(α)

		// set for all nodes
		for _, nod := range nodes {

			// find existent constraints and deactivate them
			eqx := nod.dofs[0].Eq
			eqy := nod.dofs[1].Eq
			for _, eq := range []int{eqx, eqy} {
				for _, idx := range o.Eq2idx[eq] {
					c := o.BcsTmp[idx]
					if c.Key != "rigid" {
						c.Inact = true
					}
				}
			}

			// set constraint
			o.add(key, []int{eqx, eqy}, []float64{co, si}, &fun.Zero)
		}
		return true // success
	}

	// hydraulic head
	if key == "H" {

		// get γl
		var γl float64
		if val, found := utl.Keycode(extra, "gamL"); found {
			γl = utl.Atof(val)
		} else {
			LogErrCond(true, "gamL (unit weight of liquid) must be provided when using H (hydraulic head) as essential boundary condition")
			return false // problem
		}

		// set for all nodes
		for _, nod := range nodes {

			// create function: pl(t) = γl * H(t) - γl * z
			d := nod.GetDof("pl")
			if d == nil {
				continue // node doesn't have key. ex: pl in qua8/qua4 elements
			}
			z := nod.vert.C[1] // 2D
			if o.Ndim == 3 {
				z = nod.vert.C[2] // 3D
			}
			pl := fun.Add{A: γl, Fa: fcn, B: -γl, Fb: &fun.Cte{C: z}}

			// set constraint
			o.add_single("pl", d.Eq, &pl)
		}
		return true // success
	}

	// single-point constraint
	for _, nod := range nodes {
		d := nod.GetDof(key)
		if d == nil {
			return true // success
		}
		o.add_single(key, d.Eq, fcn)
	}

	// success
	return true
}

// List returns a simple list logging bcs at time t
func (o *EssentialBcs) List(t float64) (l string) {
	for i, bc := range o.Bcs {
		if i > 0 {
			l += " "
		}
		l += utl.Sf("[%s eqs=%v f(%g)=%g]", bc.Key, bc.Eqs, t, bc.Fcn.F(t, nil))
	}
	return
}
