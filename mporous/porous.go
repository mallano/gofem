// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mporous implements models for porous media based on the Theory of Porous Media
package mporous

import "github.com/cpmech/gosl/fun"

// Model holds material parameters for porous media
type Model struct {

	// parameters
	Nf0   float64 // nf0: initial volume fraction of all fluids ~ porosity
	RhoL0 float64 // ρL0: initial liquid real density
	RhoG0 float64 // ρG0: initial gas real density
	RhoS  float64 // real (intrinsic) density of solids
	Gamlr float64 // γlr: liquid reference unit weight at temperature θini
	Gamgr float64 // γgr: gas reference unit weight at temperature θini
	Kl    float64 // liquid bulk moduli at temperature θini
	Tini  float64 // T==Θ initial temperature
	RTg   float64 // R*Θ*g: initial gas constant
	R     float64 // gas constant
	Gref  float64 // reference gravity, at time of measuring ksat, kgas
	Cl    float64 // liquid compresssibility
	Cg    float64 // gas compressibility

	// parameters for relative conductivity
	Pλ0l, Pλ1l, Pαl, Pβl float64 // liquid conductivity parameters
	Pλ0g, Pλ1g, Pαg, Pβg float64 // gas conductivity parameters

	// conductivity tensors
	Klsat [][]float64 // klsat ÷ Gref
	Kgsat [][]float64 // kgsat ÷ Gref
}

// Init initialises this structure
func (o *Model) Init(prms fun.Prms) {
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
