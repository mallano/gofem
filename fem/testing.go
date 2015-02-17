// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"encoding/json"
	"math/rand"
	"testing"

	"github.com/cpmech/gosl/num"
	"github.com/cpmech/gosl/utl"
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
	buf, err := utl.ReadFile(cmpfname)
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

		// time index
		tidx := idx + 1
		if verbose {
			utl.PfYel("\n\ntidx = %d . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n", tidx)
		}

		// load gofem results
		LogErrCond(!d.In(tidx), "TestingCompareResultsU: reading of results failed")
		if Stop() {
			tst.Errorf("reading of results failed\n")
			return
		}
		if verbose {
			utl.Pfyel("time = %v\n", d.Sol.T)
		}

		// check K matrices
		if !skipK {
			if verbose {
				utl.Pfgreen(". . . checking K matrices . . .\n")
			}
			for eid, Ksg := range cmp.Kmats {
				if e, ok := d.Elems[eid].(*ElemU); ok {
					if LogErrCond(!e.AddToKb(d.Kb, d.Sol, true), "TestingCompareResultsU: AddToKb failed") {
						tst.Errorf("AddToKb failed\n")
						break
					}
					utl.CheckMatrix(tst, utl.Sf("K%d", eid), tolK, e.K, Ksg)
				}
			}
			if Stop() {
				return
			}
		}

		// check displacements
		if verbose {
			utl.Pfgreen(". . . checking displacements . . .\n")
		}
		for nid, usg := range cmp.Disp {
			ix := d.Vid2node[nid].Dofs[0].Eq
			iy := d.Vid2node[nid].Dofs[1].Eq
			utl.CheckAnaNum(tst, "ux", tolu, d.Sol.Y[ix], usg[0], verbose)
			utl.CheckAnaNum(tst, "uy", tolu, d.Sol.Y[iy], usg[1], verbose)
			if len(usg) == 3 {
				iz := d.Vid2node[nid].Dofs[2].Eq
				utl.CheckAnaNum(tst, "uz", tolu, d.Sol.Y[iz], usg[2], verbose)
			}
		}

		// check stresses
		if true {
			if verbose {
				utl.Pfgreen(". . . checking stresses . . .\n")
			}
			for eid, sig := range cmp.Sigmas {
				if verbose {
					utl.Pforan("eid = %d\n", eid)
				}
				if e, ok := d.Cid2elem[eid].(*ElemU); ok {
					for ip, val := range sig {
						if verbose {
							utl.Pfgrey2("ip = %d\n", ip)
						}
						σ := e.States[ip].Sig
						if len(val) == 6 {
							utl.CheckAnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							utl.CheckAnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
						} else {
							utl.CheckAnaNum(tst, "sx ", tols, σ[0], val[0], verbose)
							utl.CheckAnaNum(tst, "sy ", tols, σ[1], val[1], verbose)
							utl.CheckAnaNum(tst, "sxy", tols, σ[3]/SQ2, val[2], verbose)
							if len(val) > 3 { // sx, sy, sxy, sz
								utl.CheckAnaNum(tst, "sz ", tols, σ[2], val[3], verbose)
							}
						}
					}
				}
			}
		}
	}
}

func TestConsistentTangentK(tst *testing.T, eid int, tol float64, verb bool) {

	// read summary
	sum := ReadSum()
	tidx := sum.NumTidx - 1

	// allocate domain
	d := NewDomain(Global.Sim.Regions[0])
	LogErrCond(!d.SetStage(0, Global.Sim.Stages[0]), "TestConsistentTangentK: SetStage failed")
	if Stop() {
		tst.Errorf("SetStage failed\n")
		return
	}

	// load gofem results
	LogErrCond(!d.In(tidx), "TestConsistentTangentK: reading of results failed")
	if Stop() {
		tst.Errorf("reading of results failed\n")
		return
	}

	// create copy of solution array
	ΔY := make([]float64, d.Ny)
	Y := make([]float64, d.Ny)
	copy(Y, d.Sol.Y)

	// update solution with random numbers
	for i := 0; i < d.Ny; i++ {
		d.Sol.ΔY[i] = rand.Float64()
		d.Sol.Y[i] += d.Sol.ΔY[i]
	}

	// create backup copy of all secondary variables
	for _, e := range d.ElemIntvars {
		e.BackupIvs()
	}

	// update secondary variables
	for _, e := range d.Elems {
		if !e.Update(d.Sol) {
			break
		}
	}
	if Stop() {
		tst.Errorf("Update of secondary variables failed\n")
		return
	}

	// get element structures
	var Ymap []int
	var K [][]float64
	var restorefcn func() bool
	ele := d.Elems[eid]
	if e, ok := ele.(*ElemP); ok {
		Ymap = e.Pmap
		K = e.K
		restorefcn = e.RestoreIvs
	}
	if e, ok := ele.(*ElemU); ok {
		Ymap = e.Umap
		K = e.K
		restorefcn = e.RestoreIvs
	}

	// analytical K
	if LogErrCond(!ele.AddToKb(d.Kb, d.Sol, false), "TestConsistentTangentK: AddToKb failed") {
		tst.Errorf("AddToKb failed\n")
		return
	}

	// check
	derivfcn := num.DerivCen
	var tmp float64
	for i, I := range Ymap {
		for j, J := range Ymap {
			dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, d.Sol.Y[J] = d.Sol.Y[J], x
				for l := 0; l < d.Ny; l++ {
					ΔY[l] = d.Sol.Y[J] - Y[J]
				}
				if restorefcn != nil {
					restorefcn()
				}
				ele.Update(d.Sol)
				ele.AddToRhs(d.Fb, d.Sol)
				return d.Fb[I]
			}, d.Sol.Y[J])
			utl.AnaNum(utl.Sf("K%d%d", i, j), tol, K[i][j], dnum, verb)
		}
	}
}
