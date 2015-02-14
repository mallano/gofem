// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mporous implements models for porous media based on the Theory of Porous Media
package mporous

import (
	"log"
	"math"

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
	ShowR   bool    // show residual values in Update
	AllBE   bool    // use BE for all models, including those that directly implements sl=f(pc)

	// parameters
	Nf0   float64 // nf0: initial volume fraction of all fluids ~ porosity
	RhoL0 float64 // ρL0: initial liquid real density
	RhoG0 float64 // ρG0: initial gas real density
	RhoS0 float64 // real (intrinsic) density of solids
	BulkL float64 // liquid bulk moduli at temperature θini
	RTg   float64 // R*Θ*g: initial gas constant
	Gref  float64 // reference gravity, at time of measuring ksat, kgas
	Pkl   float64 // isotrpic liquid saturated conductivity
	Pkg   float64 // isotrpic gas saturated conductivity

	// derived
	Cl    float64     // liquid compresssibility
	Cg    float64     // gas compressibility
	Klsat [][]float64 // klsat ÷ Gref
	Kgsat [][]float64 // kgsat ÷ Gref

	// conductivity and retention models
	Cnd mconduct.Model // liquid-gas conductivity models
	Lrm mreten.Model   // retention model

	// auxiliary
	nonrateLrm mreten.Nonrate // LRM is of non-rate type
}

// Init initialises this structure
func (o *Model) Init(prms fun.Prms, cnd mconduct.Model, lrm mreten.Model) (err error) {

	// conductivity and retention models
	if cnd == nil || lrm == nil {
		return utl.Err("mporous.Init: conductivity and liquid retention models must be non nil. cnd=%v, lrm=%v\n", cnd, lrm)
	}
	o.Cnd = cnd
	o.Lrm = lrm
	if m, ok := o.Lrm.(mreten.Nonrate); ok {
		o.nonrateLrm = m
	}

	// constants
	o.NmaxIt = 20
	o.Itol = 1e-8
	o.MEtrial = true

	// saturated conductivities
	var klx, kly, klz float64
	var kgx, kgy, kgz float64

	// read paramaters in
	o.RTg = 1.0
	for _, p := range prms {
		switch p.N {
		case "NmaxIt":
			o.NmaxIt = int(p.V)
		case "Itol":
			o.Itol = p.V
		case "MEtrial":
			o.MEtrial = p.V > 0
		case "ShowR":
			o.ShowR = p.V > 0
		case "AllBE":
			o.AllBE = p.V > 0
		case "nf0":
			o.Nf0 = p.V
		case "RhoL0":
			o.RhoL0 = p.V
		case "RhoG0":
			o.RhoG0 = p.V
		case "RhoS0":
			o.RhoS0 = p.V
		case "BulkL":
			o.BulkL = p.V
		case "RTg":
			o.RTg = p.V
		case "gref":
			o.Gref = p.V
		case "kl":
			o.Pkl = p.V
			klx, kly, klz = p.V, p.V, p.V
		case "kg":
			o.Pkg = p.V
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
	return
}

// GetPrms gets (an example) of parameters
func (o Model) GetPrms(example bool) fun.Prms {
	if example {
		return fun.Prms{
			&fun.Prm{N: "nf0", V: 0.3},
			&fun.Prm{N: "RhoL0", V: 1},
			&fun.Prm{N: "RhoG0", V: 0.01},
			&fun.Prm{N: "RhoS0", V: 2.7},
			&fun.Prm{N: "BulkL", V: 2.2e6},
			&fun.Prm{N: "RTg", V: 0.02},
			&fun.Prm{N: "gref", V: 10},
			&fun.Prm{N: "kl", V: 1e-3},
			&fun.Prm{N: "kg", V: 1e-2},
		}
	}
	return fun.Prms{
		&fun.Prm{N: "nf0", V: o.Nf0},
		&fun.Prm{N: "RhoL0", V: o.RhoL0},
		&fun.Prm{N: "RhoG0", V: o.RhoG0},
		&fun.Prm{N: "RhoS0", V: o.RhoS0},
		&fun.Prm{N: "BulkL", V: o.BulkL},
		&fun.Prm{N: "RTg", V: o.RTg},
		&fun.Prm{N: "gref", V: o.Gref},
		&fun.Prm{N: "kl", V: o.Pkl},
		&fun.Prm{N: "kg", V: o.Pkg},
	}
}

// InitState initialises state
func (o Model) InitState(s *StateLG, pl, pg float64) (err error) {
	s.Pl = pl
	s.Pg = pg
	s.Sl = 1
	s.Ns0 = 1 - o.Nf0
	s.RhoL = o.RhoL0
	s.RhoG = o.RhoG0
	s.Dpc = 0
	s.Wet = false
	pc := pg - pl
	if pc > 0 {
		s.Sl, err = mreten.Update(o.Lrm, 0, 1, pc)
	}
	return
}

// Update updates state
func (o Model) Update(s *StateLG, Δpl, Δpg float64) (err error) {

	// auxiliary variables
	pc0 := s.Pg - s.Pl
	sl0 := s.Sl
	Δpc := Δpg - Δpl
	pc := pc0 + Δpc
	slmin := o.Lrm.SlMin()

	// check results upon exist
	defer func() {
		if pc < 0 && s.Sl < 1 {
			err = utl.Err("inconsistent results: saturation must be equal to one when the capillary pressure is ineffective. pc = %g < 0 and sl = %g < 1 is incorrect", pc, s.Sl)
		}
		if s.Sl < slmin {
			err = utl.Err("inconsistent results: saturation must be greater than minimum saturation. sl = %g < %g is incorrect", s.Sl, slmin)
		}
	}()

	// set State with new pressure, densities and flags
	s.Pl += Δpl
	s.Pg += Δpg
	s.RhoL += o.Cl * Δpl
	s.RhoG += o.Cg * Δpg
	s.Dpc = Δpc
	s.Wet = Δpc < 0

	// exit with full liquid saturation if capillary pressure is ineffective
	if pc <= 0.0 {
		s.Sl = 1
		return
	}

	// handle simple retention models
	if o.nonrateLrm != nil && !o.AllBE {
		s.Sl = o.nonrateLrm.Sl(pc)
		return
	}

	// trial liquid saturation update
	fA, err := o.Lrm.Cc(pc0, sl0, s.Wet)
	if err != nil {
		return err
	}
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

	// fix trial sl out-of-range values
	if s.Sl < slmin {
		s.Sl = slmin
	}
	if s.Sl > 1 {
		s.Sl = 1
	}

	// message
	if o.ShowR {
		utl.PfYel("%6s%18s%18s%18s%18s%8s\n", "it", "Cc", "sl", "δsl", "r", "ex(r)")
	}

	// backward-Euler update
	var f, r, J, δsl float64
	var it int
	for it = 0; it < o.NmaxIt; it++ {
		f, err = o.Lrm.Cc(pc, s.Sl, s.Wet)
		if err != nil {
			return
		}
		r = s.Sl - sl0 - Δpc*f
		if o.ShowR {
			utl.Pfyel("%6d%18.14f%18.14f%18.14f%18.10e%8d\n", it, f, s.Sl, δsl, r, utl.Expon(r))
		}
		if math.Abs(r) < o.Itol {
			break
		}
		J, err = o.Lrm.J(pc, s.Sl, s.Wet)
		if err != nil {
			return
		}
		δsl = -r / (1.0 - Δpc*J)
		s.Sl += δsl
		if math.IsNaN(s.Sl) {
			return utl.Err("NaN found: Δpc=%v f=%v r=%v J=%v sl=%v\n", Δpc, f, r, J, s.Sl)
		}
	}

	// message
	if o.ShowR {
		utl.Pfgrey("  pc0=%.6f  sl0=%.6f  Δpl=%.6f  Δpg=%.6f  Δpc=%.6f\n", pc0, sl0, Δpl, Δpg, Δpc)
		utl.Pfgrey("  converged with %d iterations\n", it)
	}

	// check convergence
	if it == o.NmaxIt {
		return utl.Err("saturation updated failed after %d iterations\n", it)
	}
	return
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
//  getnew    -- force a new allocation; i.e. do not use any model found in database
//  Note: returns nil on errors
func GetModel(simfnk, matname string, getnew bool) *Model {

	// get new model, regardless whether it exists in database or not
	if getnew {
		return new(Model)
	}

	// search database
	key := utl.Sf("%s_%s", simfnk, matname)
	if model, ok := _models[key]; ok {
		return model
	}

	// if not found, get new
	model := new(Model)
	_models[key] = model
	return model
}

// LogModels prints to log information on existent and allocated Models
func LogModels() {
	log.Printf("\nmporous: allocated models:")
	for key, _ := range _models {
		log.Printf(" " + key)
	}
}

// _models holds pre-allocated models
var _models = map[string]*Model{}
