// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/fun"

// RjointM1 implements a 1D plasticity model for rod-joints (links/interface)
type RjointM1 struct {
}

// add model to factory
func init() {
	onedallocators["rjoint-m1"] = func() OnedSolid { return new(RjointM1) }
}

// Init initialises model
func (o *RjointM1) Init(ndim int, prms fun.Prms) (err error) {
	return
}

// GetPrms gets (an example) of parameters
func (o RjointM1) GetPrms() fun.Prms {
	return []*fun.Prm{}
}

// InitIntVars initialises internal (secondary) variables
func (o RjointM1) InitIntVars() (s *OnedState, err error) {
	s = NewOnedState(1, 0)
	return
}

// Update updates stresses for given strains
func (o *RjointM1) Update(s *OnedState, ε, Δε float64) (err error) {
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *RjointM1) CalcD(s *OnedState, firstIt bool) (D float64, err error) {
	return
}
