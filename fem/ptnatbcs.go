// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// PtNaturalBc holds information on point natural boundary conditions such as
// prescribed forces or fluxes) at nodes
type PtNaturalBc struct {
	Key   string    // key such as fux, fpl, etc...
	Eq    int       // equation
	X     []float64 // location
	Fcn   fun.Func  // function
	Extra string    // extra information
}

// PointLoads is a set of prescribed forces
type PtNaturalBcs struct {
	Eq2idx map[int]int    // maps eq number to indices in Bcs
	Bcs    []*PtNaturalBc //active boundary conditions such as prescribed forces
}

// Reset initialises internal structures
func (o *PtNaturalBcs) Reset() {
	o.Eq2idx = make(map[int]int)
	o.Bcs = make([]*PtNaturalBc, 0)
}

// AddToRhs adds the boundary conditions terms to the augmented fb vector
func (o PtNaturalBcs) AddToRhs(fb []float64, t float64) {
	for _, p := range o.Bcs {
		fb[p.Eq] += p.Fcn.F(t, p.X)
	}
}

// Set sets new point natural boundary condition data
func (o *PtNaturalBcs) Set(key string, nod *Node, fcn fun.Func, extra string) (setisok bool) {
	d := nod.GetDof(key)
	if LogErrCond(d == nil, "cannot find dof named %q", key) {
		return
	}
	if idx, ok := o.Eq2idx[d.Eq]; ok {
		o.Bcs[idx].Key = "f" + key
		o.Bcs[idx].Eq = d.Eq
		o.Bcs[idx].X = nod.Vert.C
		o.Bcs[idx].Fcn = fcn
		o.Bcs[idx].Extra = extra
	} else {
		o.Eq2idx[d.Eq] = len(o.Bcs)
		o.Bcs = append(o.Bcs, &PtNaturalBc{"f" + key, d.Eq, nod.Vert.C, fcn, extra})
	}
	return true
}

// List returns a simple list logging bcs at time t
func (o *PtNaturalBcs) List(t float64) (l string) {
	for i, bc := range o.Bcs {
		if i > 0 {
			l += " "
		}
		l += io.Sf("[%s eq=%d f(%g)=%g x=%v]", bc.Key, bc.Eq, t, bc.Fcn.F(t, bc.X), bc.X)
	}
	return
}
