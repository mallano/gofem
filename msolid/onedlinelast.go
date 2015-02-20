// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/fun"

// OnedElast implements a linear elastic model for 1D elements
type OnedElast struct {
	E float64 // Young modulus
}

// add model to factory
func init() {
	onedallocators["oned-elast"] = func() OnedModel { return new(OnedElast) }
}

// Init initialises model
func (o *OnedElast) Init(ndim int, prms fun.Prms) (err error) {
	for _, p := range prms {
		switch p.N {
		case "E":
			o.E = p.V
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o OnedElast) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "E", V: 2.0e8},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o OnedElast) InitIntVars() (s *OnedState, err error) {
	s = NewOnedState(0, 0)
	return
}

// Update updates stresses for given strains
func (o OnedElast) Update(s *OnedState, ε, Δε float64) (err error) {
	s.Sig += o.E * Δε
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o OnedElast) CalcD(s *OnedState, firstIt bool) float64 {
	return o.E
}
