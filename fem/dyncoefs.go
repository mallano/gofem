// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/utl"

	"github.com/cpmech/gofem/inp"
)

// DynCoefs calculates θ-method, Newmark's or HHT coefficients.
//  Notes:
//   θ1  -- Newmark parameter (gamma)  [0 <= θ1 <= 1]
//   θ2  -- Newmark parameter (2*beta) [0 <= θ2 <= 1]
//   HHT -- use Hilber-Hughes-Taylor method ?
//   α   -- Hilber-Hughes-Taylor parameter [-1/3 <= α <= 0]
//   if HHT==True, θ1 and θ2 are automatically calculated for unconditional stability
type DynCoefs struct {

	// input
	θ, θ1, θ2, α float64
	Ray, HHT     bool
	RayM, RayK   float64

	// derived
	β1, β2     float64
	α1, α2, α3 float64
	α4, α5, α6 float64
	α7, α8     float64
	bm, bk     float64
	hmin       float64
}

// Init initialises this structure
func (o *DynCoefs) Init(dat *inp.SolverData) {

	// hmin
	//o.hmin = dat.Hmin TODO: check this

	// HHT
	o.HHT = dat.HHT

	// θ-method
	o.θ = dat.Theta
	PanicOrNot(o.θ < 1e-5 || o.θ > 1.0, _dyncoefs_err1, o.θ)

	// HHT method
	if dat.HHT {
		o.α = dat.HHTalp
		PanicOrNot(o.α < -1.0/3.0 || o.α > 0.0, _dyncoefs_err2, o.α)
		o.θ1 = (1.0 - 2.0*o.α) / 2.0
		o.θ2 = (1.0 - o.α) * (1.0 - o.α) / 2.0

		// Newmark's method
	} else {
		o.θ1, o.θ2 = dat.Theta1, dat.Theta2
		PanicOrNot(o.θ1 < 0.0001 || o.θ1 > 1.0, _dyncoefs_err3, o.θ1)
		PanicOrNot(o.θ2 < 0.0001 || o.θ2 > 1.0, _dyncoefs_err4, o.θ2)
	}

	// Rayleigh damping
	o.RayM, o.RayK = dat.RayM, dat.RayK
	if dat.RayK > 0.0 || dat.RayM > 0.0 {
		o.Ray = true
	}
}

// CalcBoth computes betas and alphas
func (o *DynCoefs) CalcBoth(Δt float64) (err error) {
	err = o.CalcBetas(Δt)
	if err != nil {
		return
	}
	err = o.CalcAlphas(Δt)
	return
}

// CalcBetas computes only betas
func (o *DynCoefs) CalcBetas(Δt float64) (err error) {

	// timestep
	h := Δt
	if h < o.hmin {
		err = utl.Err(_dyncoefs_err5, o.hmin, h)
		return
	}

	// β coefficients
	o.β1 = 1.0 / (o.θ * h)
	o.β2 = (1.0 - o.θ) / o.θ
	return
}

// CalcAlphas computes only alphas
func (o *DynCoefs) CalcAlphas(Δt float64) (err error) {

	// timestep
	h := Δt
	if h < o.hmin {
		err = utl.Err(_dyncoefs_err6, o.hmin, h)
		return
	}

	// α coefficients
	H := h * h / 2.0
	o.α1, o.α2, o.α3 = 1.0/(o.θ2*H), h/(o.θ2*H), 1.0/o.θ2-1.0
	o.α4, o.α5, o.α6 = o.θ1*h/(o.θ2*H), 2.0*o.θ1/o.θ2-1.0, (o.θ1/o.θ2-1.0)*h

	// HHT method
	o.α7 = o.α4
	o.α8 = 1.0
	if o.HHT {
		o.α7, o.α8 = (1.0+o.α)*o.α4, 1.0+o.α
	}

	// Rayleigh damping
	o.bm, o.bk = o.α1, o.α8
	if o.Ray {
		o.bm += o.α7 * o.RayM
		o.bk += o.α7 * o.RayK
	}
	return
}

// Print prints coefficients
func (o *DynCoefs) Print() {
	utl.Pfgrey("θ=%v, θ1=%v, θ2=%v, α=%v\n", o.θ, o.θ1, o.θ2, o.α)
	utl.Pfgrey("Ray=%v, HHT=%v\n", o.Ray, o.HHT)
	utl.Pfgrey("RayM=%v, RayK=%v\n", o.RayM, o.RayK)
	utl.Pfgrey("β1=%v, β2=%v\n", o.β1, o.β2)
	utl.Pfgrey("α1=%v, α2=%v, α3=%v, α4=%v, α5=%v, α6=%v\n", o.α1, o.α2, o.α3, o.α4, o.α5, o.α6)
	utl.Pfgrey("α7=%v, α8=%v\n", o.α7, o.α8)
	utl.Pfgrey("bm=%v, bk=%v\n", o.bm, o.bk)
}

// error messages
var (
	_dyncoefs_err1 = "θ-method requires 1e-5 <= θ <= 1.0 (θ = %v is incorrect)"
	_dyncoefs_err2 = "HHT method requires: -1/3 <= α <= 0 (α = %v is incorrect)"
	_dyncoefs_err3 = "θ1 must be between 0.0001 and 1.0 (θ1 = %v is incorrect)"
	_dyncoefs_err4 = "θ2 must be between 0.0001 and 1.0 (θ2 = %v is incorrect)"
	_dyncoefs_err5 = "θ-method requires h >= %v (h = %v is incorrect)"
	_dyncoefs_err6 = "Newmark/HHT method requires h >= %v (h = %v is incorrect)"
)
