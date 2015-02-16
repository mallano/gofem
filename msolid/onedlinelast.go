// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/fun"

// OneDElast implements a linear elastic model
type OneDElast struct {
	E    float64 // Young modulus
	Nsig int
}

// add model to factory
func init() {
	rodallocators["oned-elast"] = func() RodModel { return new(OneDElast) }
}

// Init initialises model
func (o *OneDElast) Init(ndim int, prms fun.Prms) (err error) {
	o.Nsig = 1

	for _, p := range prms {
		switch p.N {
		case "E":
			o.E = p.V
		}
	}

	return
}

// GetPrms gets (an example) of parameters
func (o OneDElast) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "E", V: 2.0e8},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o OneDElast) InitIntVars() (s *State, err error) {
	s = NewState(o.Nsig, 0, 0, false)
	return
}

// Update updates stresses for given strains
func (o OneDElast) Update(s *State, ε, Δε float64) (err error) {
	s.Sig[0] += o.E * Δε
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o OneDElast) CalcD(s *State, firstIt bool) float64 {
	return o.E
}
