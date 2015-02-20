// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/fun"

// OnedLinElast implements a linear elastic model for 1D elements
type OnedLinElast struct {
	E float64 // Young modulus
}

// add model to factory
func init() {
	onedallocators["oned-elast"] = func() OnedSolid { return new(OnedLinElast) }
}

// Init initialises model
func (o *OnedLinElast) Init(ndim int, prms fun.Prms) (err error) {
	for _, p := range prms {
		switch p.N {
		case "E":
			o.E = p.V
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o OnedLinElast) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "E", V: 2.0e8},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o OnedLinElast) InitIntVars() (s *OnedState, err error) {
	s = NewOnedState(0, 0)
	return
}

// Update updates stresses for given strains
func (o OnedLinElast) Update(s *OnedState, ε, Δε float64) (err error) {
	s.Sig += o.E * Δε
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o OnedLinElast) CalcD(s *OnedState, firstIt bool) (float64, error) {
	return o.E, nil
}
