// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mporous implements models for porous media based on the Theory of Porous Media
package mporous

import (
	"log"

	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// Model holds material parameters for porous media
type Model struct {

	// constants
	NmaxIt  int     // max number iterations in Update
	Itol    float64 // iterations tolerance in Update
	MEtrial bool    // perform Modified-Euler trial to start update process

	// parameters
	Nf0   float64 // nf0: initial volume fraction of all fluids ~ porosity
	RhoL0 float64 // ρL0: initial liquid real density
	RhoG0 float64 // ρG0: initial gas real density
	RhoS  float64 // real (intrinsic) density of solids
	BulkL float64 // liquid bulk moduli at temperature θini
	RTg   float64 // R*Θ*g: initial gas constant
	Gref  float64 // reference gravity, at time of measuring ksat, kgas

	// derived
	Cl    float64     // liquid compresssibility
	Cg    float64     // gas compressibility
	Klsat [][]float64 // klsat ÷ Gref
	Kgsat [][]float64 // kgsat ÷ Gref

	// conductivity and retention models
	Cnd mconduct.Model // liquid-gas conductivity models
	Lrm mreten.Model   // retention model
}

// Init initialises this structure
func (o *Model) Init(prms fun.Prms, cnd mconduct.Model, lrm mreten.Model) (err error) {

	// constants
	o.NmaxIt = 20
	o.Itol = 1e-8
	o.MEtrial = true

	// saturated conductivities
	var klx, kly, klz float64
	var kgx, kgy, kgz float64

	// read paramaters in
	for _, p := range prms {
		switch p.N {
		case "NmaxIt":
			o.NmaxIt = int(p.V)
		case "Itol":
			o.Itol = p.V
		case "MEtrial":
			o.MEtrial = p.V > 0
		case "nf":
			o.Nf0 = p.V
		case "RhoL":
			o.RhoL0 = p.V
		case "RhoG":
			o.RhoG0 = p.V
		case "RhoS":
			o.RhoS = p.V
		case "BulkL":
			o.BulkL = p.V
		case "RTg":
			o.RTg = p.V
		case "gref":
			o.Gref = p.V
		case "kl":
			klx, kly, klz = p.V, p.V, p.V
		case "kg":
			kgx, kgy, kgz = p.V, p.V, p.V
		default:
			return utl.Err("mporous.Model: parameter named %q is incorrect\n", p.N)
		}
	}

	// derived
	o.Cl = o.RhoL0 / o.BulkL
	o.Cg = 1.0 / o.RTg
	o.Klsat = [][]float64{
		{klx / o.Gref, 0, 0},
		{0, kly / o.Gref, 0},
		{0, 0, klz / o.Gref},
	}
	o.Kgsat = [][]float64{
		{kgx / o.Gref, 0, 0},
		{0, kgy / o.Gref, 0},
		{0, 0, kgz / o.Gref},
	}

	// conductivity and retention models
	o.Cnd = cnd
	o.Lrm = lrm
	return
}

// GetPrms gets (an example) of parameters
func (o Model) GetPrms() fun.Prms {
	return fun.Prms{
		&fun.Prm{N: "nf", V: 0.3},
		&fun.Prm{N: "RhoL", V: 1},
		&fun.Prm{N: "RhoG", V: 0.01},
		&fun.Prm{N: "RhoS", V: 2.7},
		&fun.Prm{N: "BulkL", V: 2.2e6},
		&fun.Prm{N: "RTg", V: 0.02},
		&fun.Prm{N: "gref", V: 10},
		&fun.Prm{N: "kl", V: 1e-3},
		&fun.Prm{N: "kg", V: 1e-2},
	}
}

// Update updates state
func (o Model) Update(s *StateLG, Δpl, Δpg float64) error {

	// initial state and increment
	pc0 := s.Pg - s.Pl
	sl0 := s.Sl
	Δpc := Δpl - Δpg

	// set State with new pressure, densities and flags
	s.Pl += Δpl
	s.Pg += Δpg
	s.RhoL += o.Cl * Δpl
	s.RhoG += o.Cg * Δpg
	s.Dpc = Δpc
	s.Wet = Δpc < 0

	// skip if liquid-retention model is nil
	if o.Lrm == nil {
		return nil
	}

	// trial liquid saturation update
	fA, err := o.Lrm.Cc(pc0, sl0, s.Wet)
	if err != nil {
		return err
	}
	pc := pc0 + Δpc
	if o.MEtrial {
		slFE := sl0 + Δpc*fA
		fB, err := o.Lrm.Cc(pc, slFE, s.Wet)
		if err != nil {
			return err
		}
		s.Sl += 0.5 * Δpc * (fA + fB)
	} else {
		s.Sl += Δpc * fA
	}
	return nil
}

// Cc returns dsl/dpc consistent with the update method
func (o Model) Cc(pc, sl float64, wet bool, Δpc float64) float64 {
	return 0
}

// DCcDpc returns dCc/dpc consistent with the update method
func (o Model) DCcDpc(pc, sl float64, wet bool, Δpc float64) float64 {
	return 0
}

// GetModel returns (existent or new) model for porous media
//  simfnk    -- unique simulation filename key
//  matname   -- name of material
//  modelname -- model name
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetModel(simfnk, matname, modelname string, getnew bool) *Model {

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
var allocators = map[string]func() *Model{}

// _models holds pre-allocated models
var _models = map[string]*Model{}
