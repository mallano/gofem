// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mporous implements models for porous media based on the Theory of Porous Media
package mporous

import (
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// Model holds material parameters for porous media
type Model struct {

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

	// auxiliary models
	Cnd LGcond       // liquid-gas conductivity models
	Lrm mreten.Model // retention model
}

// GetModelPrms gets (an example) of parameters
func GetModelPrms() fun.Prms {
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

// Init initialises this structure
func (o *Model) Init(prms fun.Prms) (err error) {

	// read paramaters in
	var klx, kly, klz float64
	var kgx, kgy, kgz float64
	for _, p := range prms {
		switch p.N {
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

	// auxiliary models
	o.Cnd.Init(prms)
	o.Lrm.Init(prms)
	return
}

// Klr returns klr
func (o *Model) Klr(sl float64) float64 {
	return 0
}

// Kgr returns kgr
func (o *Model) Kgr(sl float64) float64 {
	return 0
}

// DklrDsl returns ∂klr/∂sl
func (o *Model) DklrDsl(sl float64) float64 {
	return 0
}

// DkgrDsg returns ∂kgr/∂sg
func (o *Model) DkgrDsg(sg float64) float64 {
	return 0
}

// Cc returns ∂sl/∂pc consistent with the update method
func (o *Model) Cc(pc, sl float64, wet bool, Δpc float64) float64 {
	return 0
}

// DCcDpc returns ∂Cc/∂pc consistent with the update method
func (o *Model) DCcDpc(pc, sl float64, wet bool, Δpc float64) float64 {
	return 0
}
