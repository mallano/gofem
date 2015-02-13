// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mconduct implements models for liquid and gas conductivity in porous media
package mconduct

import (
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// Model defines liquid-gas conductivity models
type Model interface {
	Init(prms fun.Prms) error      // Init initialises this structure
	GetPrms(example bool) fun.Prms // gets (an example) of parameters
	Klr(sl float64) float64        // Klr returns klr
	Kgr(sg float64) float64        // Kgr returns kgr
	DklrDsl(sl float64) float64    // DklrDsl returns ∂klr/∂sl
	DkgrDsg(sg float64) float64    // DkgrDsl returns ∂kgr/∂sl
}

// GetModel returns (existent or new) conductivity model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetModel(simfnk, matname, modelname string, getnew bool) Model {

	// get new model, regardless whether it exists in database or not
	if getnew {
		allocator, ok := allocators[modelname]
		if !ok {
			return nil
		}
		return allocator()
	}

	// search database
	key := utl.Sf("%s_%s_%s", simfnk, matname, modelname)
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

// LogModels prints to log information on existent and allocated Models
func LogModels() {
	log.Printf("mconduct: available models:")
	for name, _ := range allocators {
		log.Printf(" " + name)
	}
	log.Printf("\nmconduct: allocated models:")
	for key, _ := range _models {
		log.Printf(" " + key)
	}
}

// allocators holds all available models
var allocators = map[string]func() Model{}

// _models holds pre-allocated models
var _models = map[string]Model{}
