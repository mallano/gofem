// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mconduct

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// M1 implements the liquid-gas conductivity model # 1
type M1 struct {

	// parameters for liquid
	λ0l float64
	λ1l float64
	αl  float64
	βl  float64

	// parameters for gas
	λ0g float64
	λ1g float64
	αg  float64
	βg  float64

	// auxiliary functions
	klr fun.RefIncRL1
	kgr fun.RefIncRL1
}

// add model to factory
func init() {
	allocators["m1"] = func() Model { return new(M1) }
}

// GetPrms gets (an example) of parameters
func (o M1) GetPrms(example bool) fun.Prms {
	return fun.Prms{
		&fun.Prm{N: "lam0l", V: 0.001},
		&fun.Prm{N: "lam1l", V: 1.2},
		&fun.Prm{N: "alpl", V: 0.01},
		&fun.Prm{N: "betl", V: 10},
		&fun.Prm{N: "lam0g", V: 2.0},
		&fun.Prm{N: "lam1g", V: 0.001},
		&fun.Prm{N: "alpg", V: 0.01},
		&fun.Prm{N: "betg", V: 10},
	}
}

// Init initialises this structure
func (o *M1) Init(prms fun.Prms) (err error) {
	for _, p := range prms {
		switch p.N {
		case "lam0l":
			o.λ0l = p.V
		case "lam1l":
			o.λ1l = p.V
		case "alpl":
			o.αl = p.V
		case "betl":
			o.βl = p.V
		case "lam0g":
			o.λ0g = p.V
		case "lam1g":
			o.λ1g = p.V
		case "alpg":
			o.αg = p.V
		case "betg":
			o.βg = p.V
		default:
			return chk.Err("mconduct.M1: parameter named %q is incorrect\n", p.N)
		}
	}
	err = o.klr.Init(fun.Prms{
		&fun.Prm{N: "lam0", V: o.λ0l},
		&fun.Prm{N: "lam1", V: o.λ1l},
		&fun.Prm{N: "alp", V: o.αl},
		&fun.Prm{N: "bet", V: o.βl},
	})
	if err != nil {
		return
	}
	err = o.kgr.Init(fun.Prms{
		&fun.Prm{N: "lam0", V: o.λ0g},
		&fun.Prm{N: "lam1", V: o.λ1g},
		&fun.Prm{N: "alp", V: o.αg},
		&fun.Prm{N: "bet", V: o.βg},
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
