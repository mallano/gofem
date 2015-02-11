// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/utl"
)

// Driver run simulations with constitutive models for solids
type Driver struct {

	// input
	nsig  int   // number of stress components
	model Solid // solid model

	// settings
	Silent  bool    // do not show error messages
	CheckD  bool    // do check consistent matrix
	UseDfwd bool    // use DerivFwd (forward differences) instead of DerivCen (central differences) when checking D
	TolD    float64 // tolerance to check consistent matrix
	VerD    bool    // verbose check of D
	WithPC  bool    // with predictor-corrector data

	// results
	Res []*State    // stress/ivs results
	Eps [][]float64 // strains

	// for checking consistent matrix
	D [][]float64 // consistent matrix

	// for predictor-corrector plots
	precor []*State // predictor-corrector states
}

// Init initialises driver
func (o *Driver) Init(simfnk, modelname string, ndim int, pstress bool, prms fun.Prms) (err error) {
	o.nsig = 2 * ndim
	getnew := false
	o.model = GetModel(simfnk, "solidmat", modelname, getnew)
	err = o.model.Init(ndim, pstress, prms)
	if err != nil {
		return
	}
	o.D = la.MatAlloc(o.nsig, o.nsig)
	o.TolD = 1e-8
	o.VerD = true
	return
}

// Run runs simulation
func (o *Driver) Run(pth *Path) (err error) {

	// specialised models
	var sml Small
	var eup SmallStrainUpdater
	switch m := o.model.(type) {
	case Small:
		sml = m
	case SmallStrainUpdater:
		eup = m
	default:
		return utl.Err("cannot handle large-deformation models yet\n")
	}

	// allocate results arrays
	nr := 1 + (pth.Size()-1)*pth.Nincs
	if nr < 2 {
		return utl.Err(_driver_err04, pth.Size(), pth.Nincs)
	}
	o.Res = make([]*State, nr)
	o.Eps = la.MatAlloc(nr, o.nsig)
	for i := 0; i < nr; i++ {
		o.Res[i], err = o.model.InitIntVars()
		if err != nil {
			return
		}
	}

	// initial stresses and internal variables
	o.Res[0].Sig[0] = pth.MultS * pth.Sx[0]
	o.Res[0].Sig[1] = pth.MultS * pth.Sy[0]
	o.Res[0].Sig[2] = pth.MultS * pth.Sz[0]

	// auxiliary variables
	Δσ := make([]float64, o.nsig)
	Δε := make([]float64, o.nsig)

	// variables for checking D
	var tmp float64
	var εold, εnew, Δεtmp []float64
	var stmp *State
	derivfcn := num.DerivCen
	if o.CheckD {
		εold = make([]float64, o.nsig)
		εnew = make([]float64, o.nsig)
		Δεtmp = make([]float64, o.nsig)
		stmp, err = o.model.InitIntVars()
		if err != nil {
			return
		}
		if o.UseDfwd {
			derivfcn = num.DerivFwd
		}
	}

	// update states
	k := 1
	for i := 1; i < pth.Size(); i++ {

		// stress path
		if pth.UseS[i] > 0 {
			Δσ[0] = pth.MultS * (pth.Sx[i] - pth.Sx[i-1]) / float64(pth.Nincs)
			Δσ[1] = pth.MultS * (pth.Sy[i] - pth.Sy[i-1]) / float64(pth.Nincs)
			Δσ[2] = pth.MultS * (pth.Sz[i] - pth.Sz[i-1]) / float64(pth.Nincs)
			for inc := 0; inc < pth.Nincs; inc++ {

				// update
				o.Res[k].Set(o.Res[k-1])
				copy(o.Eps[k], o.Eps[k-1])
				if eup != nil {
					err = eup.StrainUpdate(o.Res[k], Δσ)
				}
				if err != nil {
					if !o.Silent {
						utl.Pfred(_driver_err01, err)
					}
					return
				}
				k += 1
			}
		}

		// strain path
		if pth.UseE[i] > 0 {
			Δε[0] = pth.MultE * (pth.Ex[i] - pth.Ex[i-1]) / float64(pth.Nincs)
			Δε[1] = pth.MultE * (pth.Ey[i] - pth.Ey[i-1]) / float64(pth.Nincs)
			Δε[2] = pth.MultE * (pth.Ez[i] - pth.Ez[i-1]) / float64(pth.Nincs)
			for inc := 0; inc < pth.Nincs; inc++ {

				// update strains
				la.VecAdd2(o.Eps[k], 1, o.Eps[k-1], 1, Δε) // εnew = εold + Δε

				// update stresses
				o.Res[k].Set(o.Res[k-1])
				err = sml.Update(o.Res[k], o.Eps[k], Δε)
				if err != nil {
					if !o.Silent {
						utl.Pfred(_driver_err02, err)
					}
					return
				}
				if o.WithPC {
					/* TODO
					tmp := o.Res[k-1].GetCopy()
					o.info.Ec.Update_σ(tmp.σ, o.Δε)
					o.precor = append(o.precor, tmp.GetCopy(), o.Res[k].GetCopy())
					*/
				}

				// check consistent matrix
				if o.CheckD {
					firstIt := false
					err = sml.CalcD(o.D, o.Res[k], firstIt)
					if err != nil {
						return utl.Err(_driver_err03, err)
					}
					copy(εold, o.Eps[k-1])
					copy(εnew, o.Eps[k])
					has_error := false
					if o.VerD {
						utl.Pf("\n")
					}
					for i := 0; i < o.nsig; i++ {
						for j := 0; j < o.nsig; j++ {
							dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
								tmp, εnew[j] = εnew[j], x
								for l := 0; l < o.nsig; l++ {
									Δεtmp[l] = εnew[l] - εold[l]
								}
								stmp.Set(o.Res[k-1])
								err = sml.Update(stmp, εnew, Δεtmp)
								res, εnew[j] = stmp.Sig[i], tmp
								return
							}, εnew[j])
							err := utl.AnaNum(utl.Sf("D[%d][%d]", i, j), o.TolD, o.D[i][j], dnum, o.VerD)
							if err != nil {
								has_error = true
							}
						}
					}
					if has_error {
						return utl.Err(_driver_err03, "ana-num comparison failed\n")
					}
				}
				k += 1
			}
		}
	}
	return
}

// error messages
var (
	_driver_err01 = "strain update failed\n%v\n"
	_driver_err02 = "stress update failed\n%v\n"
	_driver_err03 = "check of consistent matrix failed:\n %v\n\n"
	_driver_err04 = "size of path is incorrect. Size=%d, Nincs=%d\n"
)
