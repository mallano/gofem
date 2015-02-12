// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

type Derivs struct {
	DCcDpc     float64 // ∂Cc/∂pc
	DCcDsl     float64 // ∂Cc/∂sl
	D2CcDpc2   float64 // ∂²Cc/∂pc²
	D2CcDsl2   float64 // ∂²Cc/∂sl²
	D2CcDslDpc float64 // ∂²Cc/(∂pc ∂sl)
}

func (o *Derivs) SetZero() {
	o.DCcDpc, o.DCcDsl, o.D2CcDpc2, o.D2CcDsl2, o.D2CcDslDpc = 0, 0, 0, 0, 0
}

type Model interface {
	Init(prms fun.Prms) error // initialises retention model
	GetPrms() fun.Prms        // gets (an example) of parameters
}

type Direct interface {
	Sl(pc float64) float64              // compute sl directly from pc
	Cc(pc float64) float64              // compute Cc(pc) := dsl/dpc
	Derivs(d *Derivs, pc float64) error // compute ∂Cc/∂pc and ∂²Cc/∂pc²
}

type RateType interface {
	Cc(pc, sl float64, wetting bool) (float64, error)     // compute Cc(pc,sl) := dsl/dpc
	Derivs(d *Derivs, pc, sl float64, wetting bool) error // derivatives
}

// GetModel returns (existent or new) liquid retention model
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetModel(simfnk, matname, modelname string, getnew bool) Model {

	// get new model, regardless wheter it exists in database or not
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

// allocators holds all available solid models; modelname => allocator
var allocators = map[string]func() Model{}

// _models holds pre-allocated solid models (internal); key => Solid
var _models = map[string]Model{}
