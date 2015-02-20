// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mreten implements models for liquid retention curves
//  References:
//   [1] Pedroso DM, Sheng D and Zhao, J (2009) The concept of reference curves for constitutive
//       modelling in soil mechanics, Computers and Geotechnics, 36(1-2), 149-165,
//       http://dx.doi.org/10.1016/j.compgeo.2008.01.009
//   [2] Pedroso DM and Williams DJ (2010) A novel approach for modelling soil-water
//       characteristic curves with hysteresis, Computers and Geotechnics, 37(3), 374-380,
//       http://dx.doi.org/10.1016/j.compgeo.2009.12.004
//   [3] Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic
//       curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340,
//       http://dx.doi.org/10.1016/j.compgeo.2010.12.004
package mreten

import (
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/ode"
)

// Model implements a liquid retention model (LRM)
//  Derivs computes (see [1] page 618):
//    L  = ∂Cc/∂pc
//    Lx = ∂²Cc/∂pc²
//    J  = ∂Cc/∂sl
//    Jx == ∂²Cc/(∂pc ∂sl)
//    Jy == ∂²Cc/∂sl²
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
type Model interface {
	Init(prms fun.Prms) error                                              // initialises retention model
	GetPrms(example bool) fun.Prms                                         // gets (an example) of parameters
	SlMin() float64                                                        // returns sl_min
	Cc(pc, sl float64, wet bool) (float64, error)                          // computes Cc = f = ∂sl/∂pc
	L(pc, sl float64, wet bool) (float64, error)                           // computes L = ∂Cc/∂pc
	J(pc, sl float64, wet bool) (float64, error)                           // computes J = ∂Cc/∂sl
	Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) // computes all derivatives
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
		J, e := mdl.J(pc0+x*Δpc, y[0], wet)
		dfdy.Start()
		dfdy.Put(0, 0, J)
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
