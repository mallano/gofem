// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"encoding/json"
	"math"
	"testing"

	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/msolid"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/num"
)

// T_iteration testing: iteration results
type T_iteration struct {
	It     int     // iteration number
	ResRel float64 // relative residual
	Resid  float64 // absolute residual
}

// T_results testing: results
type T_results struct {
	Status     string        // status message
	LoadFactor float64       // load factor
	Iterations []T_iteration // iterations data
	Kmats      [][][]float64 // [nele][nu][nu] all stiffness matrices
	Disp       [][]float64   // [nnod][ndim] displacements at nodes
	DispMult   float64       // displacements multiplier
	Note       string        // note about number of integration points
	Sigmas     [][][]float64 // [nele][nip][nsig] all stresses @ all ips 2D:{sx, sy, sxy, sz}
}

// T_results_set is a set of comparison results
type T_results_set []*T_results

// testing_compare_results_u compares results with u-formulation
func TestingCompareResultsU(tst *testing.T, simfname, cmpfname string, tolK, tolu, tols float64, skipK, verbose bool) {

	// allocate domain
	d := NewDomain(Global.Sim.Regions[0])
	LogErrCond(!d.SetStage(0, Global.Sim.Stages[0]), "TestingCompareResultsU: SetStage failed")
	if Stop() {
		tst.Errorf("SetStage failed\n")
		return
	}

	// read file
	buf, err := io.ReadFile(cmpfname)
	LogErr(err, "TestingCompareResultsU: ReadFile failed")
	if Stop() {
		tst.Errorf("ReadFile failed\n")
		return
	}

	// unmarshal json
	var cmp_set T_results_set
	err = json.Unmarshal(buf, &cmp_set)
	LogErr(err, "TestingCompareResultsU: Unmarshal failed")
	if Stop() {
		tst.Errorf("Unmarshal failed\n")
		return
	}

	// run comparisons
	for idx, cmp := range cmp_set {

		// displacements multiplier
		dmult := 1.0
		if math.Abs(cmp.DispMult) > 1e-10 {
			dmult = cmp.DispMult
		}

		// time index
		tidx := idx + 1
		if verbose {
			io.PfYel("\n\ntidx = %d . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n", tidx)
		}

		// load gofem results
		LogErrCond(!d.In(tidx), "TestingCompareResultsU: reading of results failed")
		if Stop() {
			tst.Errorf("reading of results failed\n")
			return
		}
		if verbose {
			io.Pfyel("time = %v\n", d.Sol.T)
		}

		// check K matrices
		if !skipK {
			if verbose {
				io.Pfgreen(". . . checking K matrices . . .\n")
			}
			for eid, Ksg := range cmp.Kmats {
				if e, ok := d.Elems[eid].(*ElemU); ok {
					if LogErrCond(!e.AddToKb(d.Kb, d.Sol, true), "TestingCompareResultsU: AddToKb failed") {
						tst.Errorf("AddToKb failed\n")
						break
					}
					chk.Matrix(tst, io.Sf("K%d", eid), tolK, e.K, Ksg)
				}
			}
			if Stop() {
				return
			}
		}

		// check displacements
		if verbose {
			io.Pfgreen(". . . checking displacements . . .\n")
		}
		for nid, usg := range cmp.Disp {
			ix := d.Vid2node[nid].Dofs[0].Eq
			iy := d.Vid2node[nid].Dofs[1].Eq
			chk.AnaNum(tst, "ux", tolu, d.Sol.Y[ix], usg[0]*dmult, verbose)
			chk.AnaNum(tst, "uy", tolu, d.Sol.Y[iy], usg[1]*dmult, verbose)
			if len(usg) == 3 {
				iz := d.Vid2node[nid].Dofs[2].Eq
				chk.AnaNum(tst, "uz", tolu, d.Sol.Y[iz], usg[2]*dmult, verbose)
			}
		}

		// check stresses
		if true {
			if verbose {
				io.Pfgreen(". . . checking stresses . . .\n")
			}
			for eid, sig := range cmp.Sigmas {
				if verbose {
					io.Pforan("eid = %d\n", eid)
				}
				if e, ok := d.Cid2elem[eid].(*ElemU); ok {
					for ip, val := range sig {
						if verbose {
							io.Pfgrey2("ip = %d\n", ip)
						}
						σ := e.States[ip].Sig
						if len(val) == 6 {
							chk.AnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							chk.AnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
						} else {
							chk.AnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							chk.AnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
							chk.AnaNum(tst, "sxy", tols, σ[3]/SQ2, val[2], verbose)
							if len(val) > 3 { // sx, sy, sxy, sz
								chk.AnaNum(tst, "sz ", tols, σ[2], val[3], verbose)
							}
						}
					}
				}
			}
		}
	}
}

func TestingDefineDebugKb(tst *testing.T, eid int, tol float64, verb bool) {
	//derivfcn := num.DerivFwd
	//derivfcn := num.DerivBwd
	derivfcn := num.DerivCen
	Global.DebugKb = func(d *Domain, firstIt bool) {

		// get element structures
		var Ymap []int
		var K [][]float64
		var restorefcn func() bool
		var p_StatesBkp []*mporous.State
		var u_StatesBkp []*msolid.State
		ele := d.Elems[eid]
		if e, ok := ele.(*ElemP); ok {
			Ymap = e.Pmap
			K = e.Kpp
			nip := len(e.IpsElem)
			p_StatesBkp = make([]*mporous.State, nip)
			for i := 0; i < nip; i++ {
				p_StatesBkp[i] = e.StatesBkp[i].GetCopy()
			}
			restorefcn = func() bool {
				for i, s := range e.States {
					s.Set(p_StatesBkp[i])
				}
				return true
			}
		}
		if e, ok := ele.(*ElemU); ok {
			Ymap = e.Umap
			K = e.K
			nip := len(e.IpsElem)
			u_StatesBkp = make([]*msolid.State, nip)
			for i := 0; i < nip; i++ {
				u_StatesBkp[i] = e.StatesBkp[i].GetCopy()
			}
			restorefcn = func() bool {
				for i, s := range e.States {
					s.Set(u_StatesBkp[i])
				}
				return true
			}
		}
		if Ymap == nil || restorefcn == nil {
			tst.Errorf("TestingDefineDebugKb: cannot detect element type")
			return
		}

		// auxliary vectors
		Fbtmp := make([]float64, d.Ny)
		ΔYbkp := make([]float64, d.Ny)
		Yold := make([]float64, d.Ny)
		for i := 0; i < d.Ny; i++ {
			Yold[i] = d.Sol.Y[i] - d.Sol.ΔY[i]
		}
		//io.Pforan("Yold = %v\n", Yold)
		//io.Pforan("Ynew = %v\n", d.Sol.Y)
		//io.Pforan("ΔY   = %v\n", d.Sol.ΔY)
		copy(ΔYbkp, d.Sol.ΔY)
		defer func() {
			restorefcn()
			copy(d.Sol.ΔY, ΔYbkp)
		}()

		// check
		var tmp float64
		for i, I := range Ymap {
			for j, J := range Ymap {
				dnum := derivfcn(func(x float64, args ...interface{}) float64 {
					tmp, d.Sol.Y[J] = d.Sol.Y[J], x
					for l := 0; l < d.Ny; l++ {
						Fbtmp[l] = 0
						d.Sol.ΔY[l] = d.Sol.Y[l] - Yold[l]
					}
					restorefcn()
					ele.Update(d.Sol)
					ele.AddToRhs(Fbtmp, d.Sol)
					d.Sol.Y[J] = tmp
					return -Fbtmp[I]
				}, d.Sol.Y[J])
				if d.Sol.T > 0.7 {
					//chk.PrintAnaNum(io.Sf("K%3d%3d", i, j), tol, K[i][j], dnum, verb)
					chk.AnaNum(tst, io.Sf("K%3d%3d", i, j), tol, K[i][j], dnum, verb)
				}
			}
			//if i > 0 && d.Sol.T > 0.951 {
			//chk.Panic("stop: firstIt=%v", firstIt)
			//}
		}
	}
}
