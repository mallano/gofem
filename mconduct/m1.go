// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mconduct

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// M1 implements the liquid-gas conductivity model # 1
type M1 struct {

	// parameters for liquid
	λ0L float64
	λ1L float64
	αL  float64
	βL  float64

	// parameters for gas
	λ0G float64
	λ1G float64
	αG  float64
	βG  float64

	// auxiliary functions
	klr fun.RefIncRL1
	kgr fun.RefIncRL1
}

// add model to factory
func init() {
	allocators["m1"] = func() Model { return new(M1) }
}

// GetPrms gets (an example) of parameters
func (o M1) GetPrms() fun.Prms {
	return fun.Prms{
		&fun.Prm{N: "lam0L", V: 0.001},
		&fun.Prm{N: "lam1L", V: 1.2},
		&fun.Prm{N: "alpL", V: 0.01},
		&fun.Prm{N: "betL", V: 10},
		&fun.Prm{N: "lam0G", V: 2.0},
		&fun.Prm{N: "lam1G", V: 0.001},
		&fun.Prm{N: "alpG", V: 0.01},
		&fun.Prm{N: "betG", V: 10},
	}
}

// Init initialises this structure
func (o *M1) Init(prms fun.Prms) (err error) {
	for _, p := range prms {
		switch p.N {
		case "lam0L":
			o.λ0L = p.V
		case "lam1L":
			o.λ1L = p.V
		case "alpL":
			o.αL = p.V
		case "betL":
			o.βL = p.V
		case "lam0G":
			o.λ0G = p.V
		case "lam1G":
			o.λ1G = p.V
		case "alpG":
			o.αG = p.V
		case "betG":
			o.βG = p.V
		default:
			return utl.Err("mconduct.M1: parameter named %q is incorrect\n", p.N)
		}
	}
	o.klr.Init(fun.Prms{
		&fun.Prm{N: "lam0", V: o.λ0L},
		&fun.Prm{N: "lam1", V: o.λ1L},
		&fun.Prm{N: "alp", V: o.αL},
		&fun.Prm{N: "bet", V: o.βL},
	})
	o.kgr.Init(fun.Prms{
		&fun.Prm{N: "lam0", V: o.λ0G},
		&fun.Prm{N: "lam1", V: o.λ1G},
		&fun.Prm{N: "alp", V: o.αG},
		&fun.Prm{N: "bet", V: o.βG},
	})
	return
}

// Klr returns klr
func (o M1) Klr(sl float64) float64 {
	return o.klr.F(sl, nil)
}

// Kgr returns kgr
func (o M1) Kgr(sg float64) float64 {
	return o.kgr.F(sg, nil)
}

// DklrDsl returns ∂klr/∂sl
func (o M1) DklrDsl(sl float64) float64 {
	return o.klr.G(sl, nil)
}

// DkgrDsl returns ∂kgr/∂sg
func (o M1) DkgrDsg(sg float64) float64 {
	return o.kgr.G(sg, nil)
}
