// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

// State holds state variables for porous media with liquid and gas
//  References:
//   [1] Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis.
//       Int Journal for Numerical Methods in Engineering, 101(8) 606-634
//       http://dx.doi.org/10.1002/nme.4808
//   [2] Pedroso DM (2015) A solution to transient seepage in unsaturated porous media.
//       Computer Methods in Applied Mechanics and Engineering, 285 791-816
//       http://dx.doi.org/10.1016/j.cma.2014.12.009
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
func (o State) Lvars(m *Model) (ρl, Cpl float64, err error) {

	// n variables
	ns := (1.0 - o.Divus) * o.Ns0
	nf := 1.0 - ns
	nl := nf * o.Sl

	// ρ variables
	ρl = nl * o.RhoL

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
func (o State) Lderivs(m *Model) (ρl, Cpl, dCpldpl, dklrdpl float64, err error) {

	// n variables
	ns := (1.0 - o.Divus) * o.Ns0
	nf := 1.0 - ns
	nl := nf * o.Sl

	// ρ variables
	ρl = nl * o.RhoL

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

// LSvars calculates variables for liquid-solid simulations
func (o State) LSvars(m *Model) (ρl, ρ, p, Cpl, Cvs float64, err error) {

	// n variables; Eqs (13) and (28) of [1]
	ns := (1.0 - o.Divus) * o.Ns0
	nf := 1.0 - ns
	nl := nf * o.Sl

	// ρ variables; Eq (13) of [1]
	ρl = nl * o.RhoL
	ρs := ns * m.RhoS0
	ρ = ρl + ρs

	// pore-fluid pressure
	p = o.Pl * o.Sl // Eq. (16) of [1]

	// moduli
	Ccb, err := m.Ccb(&o)
	if err != nil {
		return
	}
	Cpl = nf * (o.Sl*m.Cl - o.RhoL*Ccb) // Eq. (32a) of [1]
	Cvs = o.Sl * o.RhoL                 // Eq. (32b) of [1]
	return
}

// LSderivs calculates derivatives for liquid-solid simulations
func (o State) LSderivs(m *Model) (ρl, ρ, Cpl, Cvs, dρdpl, dpdpl, dCpldpl, dCvsdpl, dklrdpl, dCpldusM, dρdusM float64, err error) {

	// n variables; Eqs (13) and (28) of [1]
	ns := (1.0 - o.Divus) * o.Ns0
	nf := 1.0 - ns
	nl := nf * o.Sl

	// ρ variables; Eq (13) of [1]
	ρl = nl * o.RhoL
	ρs := ns * m.RhoS0
	ρ = ρl + ρs

	// capillary pressure
	pc := o.Pg - o.Pl

	// moduli
	Ccb, err := m.Ccb(&o)
	if err != nil {
		return
	}
	Cpl = nf * (o.Sl*m.Cl - o.RhoL*Ccb) // Eq (32a) of [1]
	Cvs = o.Sl * o.RhoL                 // Eq (32b) of [1]

	// derivatives
	Ccd, err := m.Ccd(&o)
	if err != nil {
		return
	}

	// derivatives w.r.t pl
	dρdpl = nf * (o.Sl*m.Cl - o.RhoL*Ccb)      // Eq (A.9) of [1]
	dpdpl = o.Sl + pc*Ccb                      // Eq (A.11) of [1]
	dCpldpl = nf * (o.RhoL*Ccd - 2.0*Ccb*m.Cl) // Eq (A.2) of[1]
	dCvsdpl = o.Sl*m.Cl - Ccb*o.RhoL           // Eq (A.4) of [1]
	dklrdpl = -m.Cnd.DklrDsl(o.Sl) * Ccb       // Eq (A.7) of [1]

	// derivatives w.r.t us (multipliers only)
	dρdusM = (o.Sl*o.RhoL - m.RhoS0) * o.Ns0    // Eq (A.10) of [1]
	dCpldusM = (o.Sl*m.Cl - o.RhoL*Ccb) * o.Ns0 // Eq (A.3) of [1]
	return
}
