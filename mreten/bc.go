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

// BrooksCorey implements Books and Corey' model
type BrooksCorey struct {

	// parameters
	λ     float64 // slope coefficient
	pcae  float64 // air-entry pressure
	slmin float64 // residual (minium) saturation
}

// add model to factory
func init() {
	allocators["bc"] = func() Model { return new(BrooksCorey) }
}

// Init initialises model
func (o *BrooksCorey) Init(prms fun.Prms) (err error) {
	for _, p := range prms {
		switch strings.ToLower(p.N) {
		case "lam":
			o.λ = p.V
		case "pcae":
			o.pcae = p.V
		case "slmin":
			o.slmin = p.V
		default:
			return utl.Err("bc: parameter named %q is incorrect\n", p.N)
		}
	}
	return
}

// GetPrms gets (an example) of parameters
func (o BrooksCorey) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "lam", V: 0.5},
		&fun.Prm{N: "pcae", V: 0.2},
		&fun.Prm{N: "slmin", V: 0.1},
	}
}

// Sl compute sl directly from pc
func (o BrooksCorey) Sl(pc float64) float64 {
	if pc <= o.pcae {
		return 1
	}
	return o.slmin + (1-o.slmin)*math.Pow(o.pcae/pc, o.λ)
}

// Cc compute Cc(pc) := dsl/dpc
func (o BrooksCorey) Cc(pc float64) float64 {
	if pc <= o.pcae {
		return 0
	}
	return -(1 - o.slmin) * o.λ * math.Pow(o.pcae/pc, o.λ) / pc
}

// Derivs compute ∂Cc/∂pc and ∂²Cc/∂pc²
func (o BrooksCorey) Derivs(d *Derivs, pc float64) error {
	d.SetZero()
	if pc <= o.pcae {
		d.DCcDpc = (1.0 - o.slmin) * o.λ * (o.λ + 1.0) * math.Pow(o.pcae/pc, o.λ) / (pc * pc)
	}
	return nil
}
