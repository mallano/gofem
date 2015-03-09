// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
)

// State holds all continuum mechanics data, including for updating the state
type State struct {

	// essential
	Sig0 []float64 // σ0: initial Cauchy stress tensor (effective) [nsig]
	Sig  []float64 // σ: current Cauchy stress tensor (effective) [nsig]

	// for plasticity
	Alp        []float64 // α: internal variables of rate type [nalp]
	Dgam       float64   // Δγ: increment of Lagrange multiplier (for plasticity only)
	Loading    bool      // unloading flag (for plasticity only)
	ApexReturn bool      // return-to-apex (for plasticity only)

	// additional internal variables
	Phi []float64 // additional internal variables; e.g. for holding Δσ in the general stress updater

	// for large deformation
	F [][]float64 // deformation gradient [3][3]
}

// NewState allocates state structure for small or large deformation analyses
//  large  -- large deformation analyses; otherwise small strains
func NewState(nsig, nalp, nphi int, large bool) *State {
	var state State
	state.Sig0 = make([]float64, nsig)
	state.Sig = make([]float64, nsig)
	if nalp > 0 {
		state.Alp = make([]float64, nalp)
	}
	if nphi > 0 {
		state.Phi = make([]float64, nphi)
	}
	if large {
		state.F = la.MatAlloc(3, 3)
	}
	return &state
}

// Set copies states
//  Note: 1) this and other states must have been pre-allocated with the same sizes
//        2) this method does not check for errors
func (o *State) Set(other *State) {
	o.Dgam = other.Dgam
	o.Loading = other.Loading
	o.ApexReturn = other.ApexReturn
	chk.IntAssert(len(o.Sig0), len(other.Sig0))
	chk.IntAssert(len(o.Sig), len(other.Sig))
	chk.IntAssert(len(o.Alp), len(other.Alp))
	chk.IntAssert(len(o.Phi), len(other.Phi))
	copy(o.Sig0, other.Sig0)
	copy(o.Sig, other.Sig)
	copy(o.Alp, other.Alp)
	copy(o.Phi, other.Phi)
	if len(o.F) > 0 {
		la.MatCopy(o.F, 1, other.F)
	}
}

// GetCopy returns a copy of this state
func (o *State) GetCopy() *State {
	large := len(o.F) > 0
	other := NewState(len(o.Sig), len(o.Alp), len(o.Phi), large)
	other.Set(o)
	return other
}
