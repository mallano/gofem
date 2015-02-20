// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package msolid implements models for solids based on continuum mechanics
/*
 *            |   Nonrate     |    Rate
 *  ==========================|=================================
 *            |               |
 *            |  σ=f(ε)       | dσdt = f(σ,dεdt)
 *    Small   |  CalcStressE  | σ_(n+1) = σ_(n) + Δt * f_(n+1)
 *            |  D=dσ/dε|_any | StressUpdate
 *            |  ContinuousD  | D = dσ/dε_(n+1)
 *            |               | ConsistentD
 *            |               |
 *  --------------------------|---------------------------------
 *            |               |
 *    Large   | σ=f(F)        | dσdt = f(σ,F,dFdt)
 *            | CalcStressF   | D = dσdF_(n+1)
 *            | ContinuousD   |
 *            |               |
 */
package msolid

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// Solid defines the interface for solid models
type Solid interface {
	Init(ndim int, pstress bool, prms fun.Prms) error // initialises model
	GetPrms() fun.Prms                                // gets (an example) of parameters
	InitIntVars() (*State, error)                     // initialises AND allocates internal (secondary) variables
}

// Small defines rate type solid models for small strain analyses
type Small interface {
	Update(s *State, ε, Δε []float64) error            // updates stresses for given strains
	CalcD(D [][]float64, s *State, firstIt bool) error // computes D = dσ_new/dε_new consistent with StressUpdate
	ContD(D [][]float64, s *State) error               // computes D = dσ_new/dε_new continuous
}

// Large defines rate type solid models for large deformation analyses
type Large interface {
	Update(s *State, F, FΔ [][]float64) error              // updates stresses for new deformation F and FΔ
	CalcA(A [][][][]float64, s *State, firstIt bool) error // computes tangent modulus A = (2/J) * ∂τ/∂b . b - σ palm I
}

// SmallStrainUpdater define small-strain models that can update strains for given stresses
type SmallStrainUpdater interface {
	StrainUpdate(s *State, Δσ []float64) error // updates strains for given stresses (small strains formulation)
}

// GetModel returns (existent or new) solid model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetModel(simfnk, matname, modelname string, getnew bool) Solid {

	// get new model, regardless wheter it exists in database or not
	if getnew {
		allocator, ok := allocators[modelname]
		if !ok {
			return nil
		}
		return allocator()
	}

	// search database
	key := io.Sf("%s_%s_%s", simfnk, matname, modelname)
	if model, ok := _models[key]; ok {
		return model
	}

	// if not found, get new
	allocator, ok := allocators[modelname]
	if !ok {
		return nil
	}
	model := allocator()
	_models[key] = model
	return model
}

// allocators holds all available solid models; modelname => allocator
var allocators = map[string]func() Solid{}

// _models holds pre-allocated solid models (internal); key => Solid
var _models = map[string]Solid{}
