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

// Lin implements a linear retetion model: sl(pc) := 1 - λ*pc
type Lin struct {

	// parameters
	λ     float64 // slope coefficient
	pcae  float64 // air-entry pressure
	slmin float64 // residual (minium) saturation

	// derived
	pcres float64 // residual pc corresponding to slmin
}

// add model to factory
func init() {
	allocators["lin"] = func() Model { return new(Lin) }
}

// Init initialises model
func (o *Lin) Init(prms fun.Prms) (err error) {
	for _, p := range prms {
		switch strings.ToLower(p.N) {
		case "lam":
			o.λ = p.V
		case "pcae":
			o.pcae = p.V
		case "slmin":
			o.slmin = p.V
		default:
			return utl.Err("lin: parameter named %q is incorrect\n", p.N)
		}
	}
	if o.λ < 1e-15 {
		o.λ = 0
		o.pcres = math.MaxFloat64
	} else {
		o.pcres = o.pcae + (1-o.slmin)/o.λ
	}
	return
}

// GetPrms gets (an example) of parameters
func (o Lin) GetPrms(example bool) fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "lam", V: 0.5},
		&fun.Prm{N: "pcae", V: 0.2},
		&fun.Prm{N: "slmin", V: 0.1},
	}
}

// SlMin returns sl_min
func (o Lin) SlMin() float64 {
	return o.slmin
}

// Sl compute sl directly from pc
func (o Lin) Sl(pc float64) float64 {
	if pc <= o.pcae {
		return 1
	}
	if pc >= o.pcres {
		return o.slmin
	}
	return 1 - o.λ*(pc-o.pcae)
}

// Cc compute Cc(pc) := dsl/dpc
func (o Lin) Cc(pc, sl float64, wet bool) (float64, error) {
	if pc <= o.pcae || pc >= o.pcres {
		return 0, nil
	}
	return -o.λ, nil
}

// DCcDsl computes DCcDsl only
func (o Lin) DCcDsl(pc, sl float64, wet bool) (float64, error) {
	return 0, nil
}

// Derivs compute ∂Cc/∂pc and ∂²Cc/∂pc²
func (o Lin) Derivs(pc, sl float64, wet bool) error {
	D.DCcDpc, D.D2CcDpc2 = 0, 0
	return nil
}
