// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mporous implements models for porous media based on the Theory of Porous Media
package mporous

import (
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// CndModel defines liquid-gas conductivity models
type CndModel interface {
	Init(prms fun.Prms) error   // Init initialises this structure
	GetPrms() fun.Prms          // gets (an example) of parameters
	Klr(sl float64) float64     // Klr returns klr
	Kgr(sg float64) float64     // Kgr returns kgr
	DklrDsl(sl float64) float64 // DklrDsl returns ∂klr/∂sl
	DkgrDsg(sg float64) float64 // DkgrDsl returns ∂kgr/∂sl
}

// GetCndModel returns (existent or new) conductivity model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetCndModel(simfnk, matname, modelname string, getnew bool) CndModel {

	// get new model, regardless whether it exists in database or not
	if getnew {
		allocator, ok := cndallocators[modelname]
		if !ok {
			return nil
		}
		return allocator()
	}

	// search database
	key := utl.Sf("%s_%s_%s", simfnk, matname, modelname)
	if model, ok := _cndmodels[key]; ok {
		return model
	}

	// if not found, get new
	allocator, ok := cndallocators[modelname]
	if !ok {
		return nil
	}
	model := allocator()
	_cndmodels[key] = model
	return model
}

// LogCndModels prints to log information on existent and allocated CndModels
func LogCndModels() {
	log.Printf("mporous: cndmodels: available:")
	for name, _ := range cndallocators {
		log.Printf(" " + name)
	}
	log.Printf("\nmporous: cndmodels: allocated:")
	for key, _ := range _cndmodels {
		log.Printf(" " + key)
	}
}

// cndallocators holds all available conductivity modles
var cndallocators = map[string]func() CndModel{}

// _cndmodels holds pre-allocated conductivity models
var _cndmodels = map[string]CndModel{}
