// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package msolid implements models for oned elements
package msolid

import (
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// OnedSolid defines the interface for 1D models
type OnedSolid interface {
	Init(ndim int, prms fun.Prms) error                // initialises model
	GetPrms() fun.Prms                                 // gets (an example) of parameters
	InitIntVars() (*OnedState, error)                  // initialises AND allocates internal (secondary) variables
	Update(s *OnedState, ε, Δε float64) error          // update state
	CalcD(s *OnedState, firstIt bool) (float64, error) // computes D = dσ_new/dε_new consistent with StressUpdate
}

// GetOnedSolid returns (existent or new) 1D model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetOnedSolid(simfnk, matname, modelname string, getnew bool) OnedSolid {

	// get new model, regardless wheter it exists in database or not
	if getnew {
		onedallocator, ok := onedallocators[modelname]
		if !ok {
			return nil
		}
		return onedallocator()
	}

	// search database
	key := io.Sf("%s_%s_%s", simfnk, matname, modelname)
	if model, ok := _onedmodels[key]; ok {
		return model
	}

	// if not found, get new
	onedallocator, ok := onedallocators[modelname]
	if !ok {
		return nil
	}
	model := onedallocator()
	_onedmodels[key] = model
	return model
}

// onedLogModels prints to log information on existent and allocated Models
func onedLogModels() {
	l := "msolid: 1D: available:"
	for name, _ := range onedallocators {
		l += " " + name
	}
	log.Println(l)
	l = "msolid: 1D: allocated:"
	for key, _ := range _onedmodels {
		l += " " + key
	}
	log.Println(l)
}

// onedallocators holds all available oned models; modelname => allocator
var onedallocators = map[string]func() OnedSolid{}

// _models holds pre-allocated oned models (internal); key => OnedSolid
var _onedmodels = map[string]OnedSolid{}
