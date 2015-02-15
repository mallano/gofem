// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

// State holds state variables for porous media with liquid and gas
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis. Int Journal for Numerical Methods in Engineering, 101(8) 606-634 http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media. Computer Methods in Applied Mechanics and Engineering, 285 791-816 http://dx.doi.org/10.1016/j.cma.2014.12.009
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

// Lvars calculates variables for liquid-only simulations
// see Eq. (7) of [2]
func (o State) Lvars(m *Model) (Cpl float64, err error) {

	// n variables
	ns := (1.0 - o.Divus) * o.Ns0
	nf := 1.0 - ns

	// moduli
	Ccb, err := m.Ccb(&o)
	if err != nil {
		return
	}
	Cpl = nf * (o.Sl*m.Cl - o.RhoL*Ccb)
	return
}

// Lderivs calculates derivatives for liquid-only simulations
// see Eq. (A.1) of [2]
func (o State) Lderivs(m *Model) (Cpl, dCpldpl, dklrdpl float64, err error) {

	// n variables
	ns := (1.0 - o.Divus) * o.Ns0
	nf := 1.0 - ns

	// moduli
	Ccb, err := m.Ccb(&o)
	if err != nil {
		return
	}
	Cpl = nf * (o.Sl*m.Cl - o.RhoL*Ccb)

	// derivatives
	Ccd, err := m.Ccd(&o)
	if err != nil {
		return
	}
	dCpldpl = nf * (o.RhoL*Ccd - 2.0*Ccb*m.Cl)

	// conductivity model derivatives
	dklrdpl = -m.Cnd.DklrDsl(o.Sl) * Ccb
	return
}
