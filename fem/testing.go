// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"encoding/json"
	"testing"

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
	d := NewDomain(global.Sim.Regions[0])
	if !d.SetStage(0, global.Sim.Stages[0]) {
		tst.Errorf("SetStage failed\n")
	}

	// read file
	buf, err := utl.ReadFile(cmpfname)
	if LogErr(err, "TestingCompareResultsU") {
		tst.Errorf("ReadFile failed\n")
	}

	// unmarshal json
	var cmp_set T_results_set
	err = json.Unmarshal(buf, &cmp_set)
	if LogErr(err, "TestingCompareResultsU") {
		tst.Errorf("Unmarshal failed\n")
	}

	// run comparisons
	for idx, cmp := range cmp_set {

		// time index
		tidx := idx + 1
		if verbose {
			utl.PfYel("\n\ntidx = %d . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n", tidx)
		}

		// load gofem results
		if !d.In(tidx) {
			tst.Errorf("reading results failed\n")
		}
		if verbose {
			utl.Pfyel("solution.t = %v\n", d.Sol.T)
		}

		// check K matrices
		if !skipK {
			if verbose {
				utl.Pfgreen(". . . checking K matrices . . .\n")
			}
			for eid, Ksg := range cmp.Kmats {
				if e, ok := d.Elems[eid].(*ElemU); ok {
					if !e.AddToKb(d.Kb, d.Sol, true) {
						tst.Errorf("AddToKb failed\n")
					}
					utl.CheckMatrix(tst, utl.Sf("K%d", eid), tolK, e.K, Ksg)
				}
			}
		}

		// check displacements
		if verbose {
			utl.Pfgreen(". . . checking displacements . . .\n")
		}
		for nid, usg := range cmp.Disp {
			ix := d.Vid2node[nid].dofs[0].Eq
			iy := d.Vid2node[nid].dofs[1].Eq
			utl.CheckAnaNum(tst, "ux", tolu, d.Sol.Y[ix], usg[0], verbose)
			utl.CheckAnaNum(tst, "uy", tolu, d.Sol.Y[iy], usg[1], verbose)
			if len(usg) == 3 {
				iz := d.Vid2node[nid].dofs[2].Eq
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
