// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

// State holds state variables for porous media with liquid and gas
type State struct {
	Pl    float64 // pl: liquid pressure
	Pg    float64 // pg: gas pressure
	Divus float64 // div(us): divergence of solids displacement
	Sl    float64 // sl: liquid saturation
	Ns0   float64 // ns0: initial partial fraction of solids
	RhoL  float64 // ρL: real (intrinsic) density of liquid
	RhoG  float64 // ρG: real (intrinsic) density of gas
	Dpc   float64 // Δpc: step increment of capillary pressure
	Wet   bool    // wetting flag
}

// GetCopy returns a copy of State
func (o State) GetCopy() *State {
	return &State{
		o.Pl,
		o.Pg,
		o.Divus,
		o.Sl,
		o.Ns0,
		o.RhoL,
		o.RhoG,
		o.Dpc,
		o.Wet,
	}
}

// Set sets this State with another State
func (o *State) Set(another *State) {
	o.Pl = another.Pl
	o.Pg = another.Pg
	o.Divus = another.Divus
	o.Sl = another.Sl
	o.Ns0 = another.Ns0
	o.RhoL = another.RhoL
	o.RhoG = another.RhoG
	o.Dpc = another.Dpc
	o.Wet = another.Wet
}
