// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"math"
	"strings"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// VanGen implements van Genuchten's model
type VanGen struct {

	// parameters
	α, m, n float64 // parameters
	slmin   float64 // minimum sl
	pcmin   float64 // pc limit to consider zero slope
}

// add model to factory
func init() {
	allocators["vg"] = func() Model { return new(VanGen) }
}

// Init initialises model
func (o *VanGen) Init(prms fun.Prms) (err error) {
	o.pcmin = 1e-3
	for _, p := range prms {
		switch strings.ToLower(p.N) {
		case "alp":
			o.α = p.V
		case "m":
			o.m = p.V
		case "n":
			o.n = p.V
		case "slmin":
			o.slmin = p.V
		case "pcmin":
			o.pcmin = p.V
		default:
			return utl.Err("vg: parameter named %q is incorrect\n", p.N)
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o VanGen) GetPrms(example bool) fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "alp", V: 0.08},
		&fun.Prm{N: "m", V: 4},
		&fun.Prm{N: "n", V: 4},
		&fun.Prm{N: "slmin", V: 0.01},
		&fun.Prm{N: "pcmin", V: 1e-3},
	}
}

// SlMin returns sl_min
func (o VanGen) SlMin() float64 {
	return o.slmin
}

// Sl computes sl directly from pc
func (o VanGen) Sl(pc float64) float64 {
	if pc <= o.pcmin {
		return 1
	}
	c := math.Pow(o.α*pc, o.n)
	fac := 1.0 - o.slmin
	return fac * math.Pow(1+c, -o.m)
}

// Cc computes Cc(pc) := dsl/dpc
func (o VanGen) Cc(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcmin {
		return 0, nil
	}
	c := math.Pow(o.α*pc, o.n)
	fac := 1.0 - o.slmin
	return -fac * c * math.Pow(c+1.0, -o.m-1.0) * o.m * o.n / pc, nil
}

// J computes J = ∂Cc/∂sl
func (o VanGen) J(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// Derivs compute ∂Cc/∂pc and ∂²Cc/∂pc²
func (o VanGen) Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) {
	if pc <= o.pcmin {
		return
	}
	c := math.Pow(o.α*pc, o.n)
	d := math.Pow(o.α*pc, o.n*2.0)
	mm := o.m * o.m
	nn := o.n * o.n
	mn := o.m * o.n
	ppp := pc * pc * pc
	fac := 1.0 - o.slmin
	L = fac * c * math.Pow(c+1.0, -o.m-2.0) * mn * (c*mn - o.n + c + 1.0) / (pc * pc)
	Lx = -fac * c * math.Pow(c+1.0, -o.m-3.0) * mn * (d*mm*nn - 3.0*c*o.m*nn - c*nn + nn + 3.0*d*mn + 3.0*c*mn - 3.0*c*o.n - 3.0*o.n + 2.0*d + 4.0*c + 2.0) / ppp
	return
}
