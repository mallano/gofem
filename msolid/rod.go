// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package msolid implements models for rod elements
package msolid

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// RodModel defines the interface for rod models
type RodModel interface {
	Init(ndim int, prms fun.Prms) error   // initialises model
	GetPrms() fun.Prms                    // gets (an example) of parameters
	CalcD(s *State, firstIt bool) float64 // computes D = dσ_new/dε_new consistent with StressUpdate
	Update(s *State, ε, Δε float64) (err error)
	InitIntVars() (*State, error) // initialises AND allocates internal (secondary) variables
}

// GetModel returns (existent or new) rod model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetRodModel(simfnk, matname, modelname string, getnew bool) RodModel {

	// get new model, regardless wheter it exists in database or not
	if getnew {
		rodallocator, ok := rodallocators[modelname]
		if !ok {
			return nil
		}
		return rodallocator()
	}

	// search database
	key := utl.Sf("%s_%s_%s", simfnk, matname, modelname)
	if model, ok := _rodmodels[key]; ok {
		return model
	}

	// if not found, get new
	rodallocator, ok := rodallocators[modelname]
	if !ok {
		return nil
	}
	model := rodallocator()
	_rodmodels[key] = model
	return model
}

// rodallocators holds all available rod models; modelname => allocator
var rodallocators = map[string]func() RodModel{}

// _models holds pre-allocated rod models (internal); key => RodModel
var _rodmodels = map[string]RodModel{}
