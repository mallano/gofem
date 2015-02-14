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
func (o BrooksCorey) GetPrms(example bool) fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "lam", V: 0.5},
		&fun.Prm{N: "pcae", V: 0.2},
		&fun.Prm{N: "slmin", V: 0.1},
	}
}

// SlMin returns sl_min
func (o BrooksCorey) SlMin() float64 {
	return o.slmin
}

// Sl compute sl directly from pc
func (o BrooksCorey) Sl(pc float64) float64 {
	if pc <= o.pcae {
		return 1
	}
	return o.slmin + (1-o.slmin)*math.Pow(o.pcae/pc, o.λ)
}

// Cc compute Cc(pc) := dsl/dpc
func (o BrooksCorey) Cc(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcae {
		return 0, nil
	}
	return -(1 - o.slmin) * o.λ * math.Pow(o.pcae/pc, o.λ) / pc, nil
}

// DCcDsl computes DCcDsl only
func (o BrooksCorey) DCcDsl(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// Derivs compute ∂Cc/∂pc and ∂²Cc/∂pc²
func (o BrooksCorey) Derivs(pc, sl float64, wet bool) error {
	D.DCcDpc, D.D2CcDpc2 = 0, 0
	if pc <= o.pcae {
		return nil
	}
	cf := (1.0 - o.slmin) * o.λ
	pc2 := pc * pc
	pp := math.Pow(o.pcae/pc, o.λ)
	dppdpc := -o.λ * math.Pow(o.pcae/pc, o.λ) / pc
	D.DCcDpc = cf * (o.λ + 1.0) * pp / pc2
	D.D2CcDpc2 = cf * (o.λ + 1.0) * (dppdpc - 2.0*pp/pc) / pc2
	return nil
}
