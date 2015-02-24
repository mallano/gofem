// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
)

// KGcalculator defines calculators of elasticity coefficients K and G
type KGcalculator interface {
	Init(prms fun.Prms) (err error)
	Calc(s *State) (K, G float64)
}

// kgcfactory holds KG calculators
var kgcfactory = map[string]func() KGcalculator{}

// GetKgc returns a KG calculator
// It returns nil on errors
func GetKgc(name string, prms fun.Prms) KGcalculator {
	allocator, ok := kgcfactory[name]
	if !ok {
		return nil
	}
	o := allocator()
	err := o.Init(prms)
	if err != nil {
		return nil
	}
	return o
}

// SmallElasticity implements linear/non-linear elasticity for small strain analyses
type SmallElasticity struct {
	Nsig  int          // number of stress components
	E, Nu float64      // Young modulus and Poisson coefficient
	L, G  float64      // Lame's coefficients. L == λ, G == μ
	K     float64      // Bulk modulus
	Pse   bool         // is plane-stress?
	Kgc   KGcalculator // K and G calculator for non-linear models
}

// Init initialises this structure
func (o *SmallElasticity) Init(ndim int, pstress bool, prms fun.Prms) (err error) {
	o.Nsig = 2 * ndim
	o.Pse = pstress
	var has_E, has_ν, has_l, has_G, has_K bool
	for _, p := range prms {
		switch p.N {
		case "E":
			o.E, has_E = p.V, true
		case "nu":
			o.Nu, has_ν = p.V, true
		case "l":
			o.L, has_l = p.V, true
		case "G":
			o.G, has_G = p.V, true
		case "K":
			o.K, has_K = p.V, true
		}
		if skgc, found := io.Keycode(p.Extra, "kgc"); found {
			o.Kgc = GetKgc(skgc, prms)
			if o.Kgc == nil {
				return chk.Err("cannot find kgc model named %s", skgc)
			}
			err = o.Kgc.Init(prms)
			if err != nil {
				return
			}
		}
	}
	switch {
	case has_E && has_ν:
		o.L = Calc_l_from_Enu(o.E, o.Nu)
		o.G = Calc_G_from_Enu(o.E, o.Nu)
		o.K = Calc_K_from_Enu(o.E, o.Nu)
	case has_l && has_G:
		o.E = Calc_E_from_lG(o.L, o.G)
		o.Nu = Calc_nu_from_lG(o.L, o.G)
		o.K = Calc_K_from_lG(o.L, o.G)
	case has_K && has_G:
		o.E = Calc_E_from_KG(o.K, o.G)
		o.Nu = Calc_nu_from_KG(o.K, o.G)
		o.L = Calc_l_from_KG(o.K, o.G)
	case has_K && has_ν:
		o.E = Calc_E_from_Knu(o.K, o.Nu)
		o.G = Calc_G_from_Knu(o.K, o.Nu)
		o.L = Calc_l_from_Knu(o.K, o.Nu)
	default:
		return chk.Err("combination of Elastic constants is incorrect. options are {E,nu}, {l,G}, {K,G} and {K,nu}\n")
	}
	return
}

// GetPrms gets (an example) of parameters
func (o SmallElasticity) GetPrms() fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "E", V: o.E},
		&fun.Prm{N: "nu", V: o.Nu},
	}
}

// Update computes new stresses for new strains ε
func (o SmallElasticity) Update(s *State, σ0, ε []float64) (err error) {
	σ := s.Sig
	defer func() {
		for i := 0; i < len(σ0); i++ {
			σ[i] += σ0[i]
		}
	}()
	if o.Pse {
		c := o.E / (1.0 - o.Nu*o.Nu)
		σ[0] = c * (ε[0] + o.Nu*ε[1])
		σ[1] = c * (o.Nu*ε[0] + ε[1])
		σ[2] = 0
		σ[3] = c * (1.0 - o.Nu) * ε[3]
		return
	}
	trε := ε[0] + ε[1] + ε[2]
	for i := 0; i < o.Nsig; i++ {
		σ[i] = o.L*trε*tsr.Im[i] + 2.0*o.G*ε[i]
	}
	return
}

// CalcD computes D = dσ_new/dε_new (consistent)
func (o SmallElasticity) CalcD(D [][]float64, s *State) (err error) {
	if o.Pse {
		if o.Nsig != 4 {
			return chk.Err("for plane-stress analyses, D must be 4x4. nsig = %d is incorrect.\n", o.Nsig)
		}
		if o.Kgc != nil {
			return chk.Err("plane-stress analysis does not work with nonlinear K and G\n")
		}
		c := o.E / (1.0 - o.Nu*o.Nu)
		la.MatFill(D, 0)
		D[0][0] = c
		D[0][1] = c * o.Nu
		D[1][0] = c * o.Nu
		D[1][1] = c
		D[3][3] = c * (1.0 - o.Nu)
		return
	}
	if o.Kgc != nil {
		o.K, o.G = o.Kgc.Calc(s)
	}
	for i := 0; i < o.Nsig; i++ {
		for j := 0; j < o.Nsig; j++ {
			D[i][j] = o.K*tsr.Im[i]*tsr.Im[j] + 2*o.G*tsr.Psd[i][j]
		}
	}
	return
}

// converters ///////////////////////////////////////////////////////////////////////////////////////

// -- E, ν -----------------------------------------------------

// Calc_l_from_Enu returns l given E and ν
func Calc_l_from_Enu(E, ν float64) float64 {
	return E * ν / ((1.0 + ν) * (1.0 - 2.0*ν))
}

// Calc_G_from_Enu returns G given E and ν. NOTE: G == μ
func Calc_G_from_Enu(E, ν float64) float64 {
	return E / (2.0 * (1.0 + ν))
}

// Calc_K_from_Enu returns K given E and ν
func Calc_K_from_Enu(E, ν float64) float64 {
	return E / (3.0 * (1.0 - 2.0*ν))
}

// -- l, G -----------------------------------------------------

// Calc_E_from_lG returns E given l and G
func Calc_E_from_lG(l, G float64) float64 {
	return G * (3.0*l + 2.0*G) / (l + G)
}

// Calc_nu_from_lG returns ν given l and G
func Calc_nu_from_lG(l, G float64) float64 {
	return 0.5 * l / (l + G)
}

// Calc_K_from_lG returns K given l and G
func Calc_K_from_lG(l, G float64) float64 {
	return l + 2.0*G/3.0
}

// -- K, G -----------------------------------------------------

// Calc_E_from_KG returns E given K and G
func Calc_E_from_KG(K, G float64) float64 {
	return 9.0 * K * G / (3.0*K + G)
}

// Calc_nu_from_KG returns ν given K and G
func Calc_nu_from_KG(K, G float64) float64 {
	return (3.0*K - 2.0*G) / (6.0*K + 2.0*G)
}

// Calc_l_from_KG returns l given K and G
func Calc_l_from_KG(K, G float64) float64 {
	return K - 2.0*G/3.0
}

// -- K, ν -----------------------------------------------------

// Calc_E_from_Knu returns E given K and ν
func Calc_E_from_Knu(K, ν float64) float64 {
	return 3.0 * K * (1.0 - 2.0*ν)
}

// Calc_G_from_Kν returns G given K and ν
func Calc_G_from_Knu(K, ν float64) float64 {
	return 3.0 * (1.0 - 2.0*ν) * K / (2.0 * (1.0 + ν))
}

// Calc_l_from_Kν returns l given K and ν
func Calc_l_from_Knu(K, ν float64) float64 {
	return 3.0 * K * ν / (1.0 + ν)
}
