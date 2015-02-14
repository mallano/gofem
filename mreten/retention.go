// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
	"github.com/cpmech/gosl/utl"
)

// D holds all derivatives related to liquid retention models
var D struct {
	DCcDpc     float64 // ∂Cc/∂pc
	DCcDsl     float64 // ∂Cc/∂sl
	D2CcDpc2   float64 // ∂²Cc/∂pc²
	D2CcDsl2   float64 // ∂²Cc/∂sl²
	D2CcDpcDsl float64 // ∂²Cc/(∂pc ∂sl)
}

// Model implements a liquid retention model (LRM)
type Model interface {
	Init(prms fun.Prms) error                         // initialises retention model
	GetPrms(example bool) fun.Prms                    // gets (an example) of parameters
	SlMin() float64                                   // returns sl_min
	Cc(pc, sl float64, wet bool) (float64, error)     // compute Cc(pc,sl) := dsl/dpc
	DCcDsl(pc, sl float64, wet bool) (float64, error) // computes DCcDsl only
	Derivs(pc, sl float64, wet bool) error            // derivatives
}

// Nonrate is a subset of LRM that directly computes saturation from capillary pressure
type Nonrate interface {
	Sl(pc float64) float64 // compute sl directly from pc
}

// Update updates pc and sl for given Δpc. An implicit ODE solver is used.
func Update(mdl Model, pc0, sl0, Δpc float64) (slNew float64, err error) {

	// wetting flag
	wet := Δpc < 0

	// callback functions
	//   x      = [0.0, 1.0]
	//   pc     = pc0 + x * Δpc
	//   y[0]   = sl
	//   f(x,y) = dy/dx = dsl/dpc * dpc/dx = Cc * Δpc
	//   J(x,y) = df/dy = DCcDsl * Δpc
	fcn := func(f []float64, x float64, y []float64, args ...interface{}) (e error) {
		f[0], e = mdl.Cc(pc0+x*Δpc, y[0], wet)
		f[0] *= Δpc
		return nil
	}
	jac := func(dfdy *la.Triplet, x float64, y []float64, args ...interface{}) (e error) {
		if dfdy.Max() == 0 {
			dfdy.Init(1, 1, 1)
		}
		e = mdl.Derivs(pc0+x*Δpc, y[0], wet)
		dfdy.Start()
		dfdy.Put(0, 0, D.DCcDsl)
		return
	}

	// ode solver
	var odesol ode.ODE
	odesol.Init("Radau5", 1, fcn, jac, nil, nil, true)
	odesol.SetTol(1e-10, 1e-7)

	// solve
	y := []float64{sl0}
	err = odesol.Solve(y, 0, 1, 1, false)
	slNew = y[0]
	return
}

// GetModel returns (existent or new) liquid retention model
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
	log.Printf("mporous: available models:")
	for name, _ := range allocators {
		log.Printf(" " + name)
	}
	log.Printf("\nmporous: allocated models:")
	for key, _ := range _models {
		log.Printf(" " + key)
	}
}

// allocators holds all available models
var allocators = map[string]func() Model{}

// _models holds pre-allocated models
var _models = map[string]Model{}
