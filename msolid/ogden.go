// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"math"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/tsr"
	"github.com/cpmech/gosl/utl"
)

// Ogden implements a linear elastic model
type Ogden struct {

	// basic data
	Nsig int // number of stress components

	// parameters
	Alp []float64 // α parameters
	Mu  []float64 // μ parameters
	K   float64   // Bulk modulus

	// auxiliary
	Fi   [][]float64 // inverse of F [3][3]
	J    float64     // det(F)
	b    [][]float64 // left Cauchy-Green deformation [3][3]
	bm   []float64   // Mandel version of b
	λ    []float64   // eigenvalues of b [3]
	P    [][]float64 // eigenprojectors of b [3][nsig]
	τ    []float64   // eigenvalues Kirchhoff stress [3]
	dτdb [][]float64
}

// add model to factory
func init() {
	allocators["ogden"] = func() Solid { return new(Ogden) }
}

// Init initialises model
func (o *Ogden) Init(ndim int, pstress bool, prms fun.Prms) (err error) {

	// basic data
	o.Nsig = 2 * ndim

	// parameters
	for _, p := range prms {
		if p.N == "K" {
			o.K = p.V
		}
		if p.N[:3] == "alp" {
			o.Alp = append(o.Alp, p.V)
		}
		if p.N[:3] == "mu" {
			o.Mu = append(o.Mu, p.V)
		}
	}
	if len(o.Alp) != len(o.Mu) {
		return utl.Err("number of alp must be equal to number of mu. %d != %d\n", len(o.Alp), len(o.Mu))
	}

	// auxiliary
	o.Fi = tsr.Alloc2()
	o.b = tsr.Alloc2()
	o.bm = make([]float64, o.Nsig)
	o.λ = make([]float64, 3)
	o.P = tsr.M_AllocEigenprojs(o.Nsig)
	o.τ = make([]float64, 3)
	return
}

// GetPrms gets (an example) of parameters
func (o Ogden) GetPrms() fun.Prms {
	return []*fun.Prm{}
}

// InitIntVars initialises internal (secondary) variables
func (o Ogden) InitIntVars() (s *State, err error) {
	s = NewState(o.Nsig, 0, 0, true)
	return
}

// Update updates stresses for given strains
func (o *Ogden) Update(s *State, F [][]float64) (err error) {

	// TODO
	return utl.Err("Ogden model is not implemented yet")

	// spectral decomposition
	err = o.b_and_spectral_decomp(F)
	if err != nil {
		return
	}

	// updated principal Kirchhoff stress
	lnJ := math.Log(o.J)
	for i := 0; i < 3; i++ {
		o.τ[i] = 0
		for p, α := range o.Alp {
			f := (math.Pow(o.λ[0], α) + math.Pow(o.λ[1], α) + math.Pow(o.λ[2], α)) / 3.0
			o.τ[i] += o.Mu[p]*math.Pow(o.J, -α/3.0)*(math.Pow(o.λ[i], α)-f) + o.K*lnJ
		}
	}

	// assemble Cauchy stress
	for i := 0; i < len(s.Sig); i++ {
		s.Sig[i] = (o.τ[0]*o.P[0][i] + o.τ[1]*o.P[1][i] + o.τ[2]*o.P[2][i]) / o.J
	}
	return
}

// CalcA computes tangent modulus A = (2/J) * ∂τ/∂b . b - σ palm I
func (o *Ogden) CalcA(A [][][][]float64, s *State, firstIt bool) (err error) {

	// TODO
	return utl.Err("Ogden model is not implemented yet")

	// spectral decomposition
	err = o.b_and_spectral_decomp(s.F)
	if err != nil {
		return
	}

	// recover principal Kirchhoff
	σ := s.Sig
	for i := 0; i < 3; i++ {
		o.τ[i] = 0
		for j := 0; j < o.Nsig; j++ {
			o.τ[i] += o.J * σ[j] * o.P[i][j]
		}
	}

	// derivatives
	var cf float64
	for _, α := range o.Alp {
		f := (math.Pow(o.λ[0], α) + math.Pow(o.λ[1], α) + math.Pow(o.λ[2], α)) / 3.0
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				// TODO
				o.dτdb[i][j] = (cf)*(f-math.Pow(o.λ[i], α)-math.Pow(o.λ[j], α)) + o.K/(2*o.λ[j]*o.λ[j])
			}
		}
	}

	// assemble

	// compute spatial tangent modulus
	return
}

// spectral_decomp computes the spectral decomposition of b := F*tr(F) tensor
func (o *Ogden) b_and_spectral_decomp(F [][]float64) (err error) {

	// determinant of F
	o.J, err = tsr.Inv(o.Fi, F)
	if err != nil {
		return
	}

	// left Cauchy-Green tensor
	tsr.LeftCauchyGreenDef(o.b, F)

	// eigenvalues and eigenprojectors
	tsr.Ten2Man(o.bm, o.b)
	err = tsr.M_EigenValsProjsNum(o.P, o.λ, o.bm)
	if err != nil {
		return
	}
	o.λ[0] = math.Sqrt(o.λ[0])
	o.λ[1] = math.Sqrt(o.λ[1])
	o.λ[2] = math.Sqrt(o.λ[2])
	return
}
