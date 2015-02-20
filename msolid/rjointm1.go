// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/fun"
)

// RjointM1 implements a 1D plasticity model for rod-joints (links/interface)
type RjointM1 struct {
	ks  float64 // elasticity constant
	τy0 float64 // initial yield stress
	kh  float64 // hardening modulus
	μ   float64 // frictioin coefficient
}

// Init initialises model
func (o *RjointM1) Init(prms fun.Prms) (err error) {
	for _, p := range prms {
		switch p.N {
		case "ks":
			o.ks = p.V
		case "tauy0":
			o.τy0 = p.V
		case "kh":
			o.kh = p.V
		case "mu":
			o.μ = p.V
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o RjointM1) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "ks", V: 1e4},
		&fun.Prm{N: "tauy0", V: 20},
		&fun.Prm{N: "kh", V: 0},
		&fun.Prm{N: "mu", V: 0.5},
	}
}

// InitIntVars initialises internal (secondary) variables
func (o RjointM1) InitIntVars() (s *OnedState, err error) {
	s = NewOnedState(1, 2) // 1:{ωpb}  2:{qn1,qn2}
	return
}

// Update updates stresses for given strains
func (o *RjointM1) Update(s *OnedState, σcNew, Δω float64) (err error) {

	// limit σcNew
	if σcNew < 0 {
		σcNew = 0
	}

	// internal values
	τ := &s.Sig
	ωpb := &s.Alp[0]

	// trial stress
	τ_tr := (*τ) + o.ks*Δω
	f_tr := math.Abs(τ_tr) - (o.τy0 + o.kh*(*ωpb) + o.μ*σcNew)

	// elastic update
	if f_tr <= 0.0 {
		*τ = τ_tr
		s.Loading = false
		return
	}

	// plastic update
	Δγ := f_tr / (o.ks + o.kh)
	*τ = τ_tr - o.ks*Δγ*fun.Sign(τ_tr)
	*ωpb += Δγ
	s.Loading = true
	return
}

// CalcD computes D = dσ_new/dε_new consistent with StressUpdate
func (o *RjointM1) CalcD(s *OnedState, firstIt bool) (DτDω, DτDσc float64, err error) {

	// elastic
	if !s.Loading {
		return o.ks, 0, nil
	}

	// plastic
	τ := s.Sig
	DτDω = o.ks * o.kh / (o.ks + o.kh)
	DτDσc = o.ks * o.μ * fun.Sign(τ) / (o.ks + o.kh)
	return
}
