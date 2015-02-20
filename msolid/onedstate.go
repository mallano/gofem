// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import "github.com/cpmech/gosl/chk"

// OnedState holds data for 1D models
type OnedState struct {

	// essential
	Sig float64 // σ: Cauchy stress component

	// for plasticity
	Alp     []float64 // α: internal variables of rate type [nalp]
	Dgam    float64   // Δγ: increment of Lagrange multiplier (for plasticity only)
	Loading bool      // unloading flag (for plasticity only)

	// additional internal variables
	Phi []float64 // additional internal variables; e.g. for holding Δσ in the general stress updater

	// for large deformation
	F float64 // deformation gradient
}

// NewOnedState allocates 1D state structure for small or large deformation analyses
//  large  -- large deformation analyses; otherwise small strains
func NewOnedState(nalp, nphi int) *OnedState {
	var state OnedState
	if nalp > 0 {
		state.Alp = make([]float64, nalp)
	}
	if nphi > 0 {
		state.Phi = make([]float64, nphi)
	}
	return &state
}

// Set copies states
//  Note: 1) this and other states must have been pre-allocated with the same sizes
//        2) this method does not check for errors
func (o *OnedState) Set(other *OnedState) {
	o.Sig = other.Sig
	o.Dgam = other.Dgam
	o.Loading = other.Loading
	chk.IntAssert(len(o.Alp), len(other.Alp))
	chk.IntAssert(len(o.Phi), len(other.Phi))
	copy(o.Alp, other.Alp)
	copy(o.Phi, other.Phi)
	o.F = other.F
}

// GetCopy returns a copy of this state
func (o *OnedState) GetCopy() *OnedState {
	other := NewOnedState(len(o.Alp), len(o.Phi))
	other.Set(o)
	return other
}
