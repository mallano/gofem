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

	// read summary
	if !Global.Sum.Read() {
		tst.Error("cannot read summary file for simulation=%q\n", simfname)
		return
	}

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
	dmult := 1.0
	for idx, cmp := range cmp_set {

		// displacements multiplier
		if idx == 0 && math.Abs(cmp.DispMult) > 1e-10 {
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

func TestingDefineDebugKbP(tst *testing.T, eid int, tmin, tmax, tol float64, verb bool) {
	derivfcn := num.DerivCen
	Global.DebugKb = func(d *Domain, it int) {
		elem := d.Elems[eid]
		if ele, ok := elem.(*ElemP); ok {

			// skip other times
			if tmin >= 0 && tmax >= 0 {
				if d.Sol.T < tmin || d.Sol.T > tmax {
					return
				}
			}

			// message
			if verb {
				io.PfYel("\nit=%2d t=%v\n", it, d.Sol.T)
			}

			// copy states and solution
			nip := len(ele.IpsElem)
			states := make([]*mporous.State, nip)
			statesBkp := make([]*mporous.State, nip)
			for i := 0; i < nip; i++ {
				states[i] = ele.States[i].GetCopy()
				statesBkp[i] = ele.StatesBkp[i].GetCopy()
			}
			Fbtmp, ΔYbkp, Yold := testing_get_aux_vectors(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					ele.States[i].Set(states[i])
					ele.StatesBkp[i].Set(statesBkp[i])
				}
				copy(d.Sol.ΔY, ΔYbkp)
			}()

			// check
			var tmp float64
			for i, I := range ele.Pmap[:1] {
				for j, J := range ele.Pmap[:1] {
					dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
						tmp, d.Sol.Y[J] = d.Sol.Y[J], x
						for k := 0; k < d.Ny; k++ {
							Fbtmp[k] = 0
							d.Sol.ΔY[k] = d.Sol.Y[k] - Yold[k]
						}
						if it == 0 {
							for k := 0; k < nip; k++ {
								ele.States[k].Set(states[k])
							}
						} else {
							for k := 0; k < nip; k++ {
								ele.States[k].Set(statesBkp[k])
							}
						}
						ele.Update(d.Sol)
						ele.AddToRhs(Fbtmp, d.Sol)
						res = -Fbtmp[I]
						d.Sol.Y[J] = tmp
						return res
					}, d.Sol.Y[J])
					chk.AnaNum(tst, io.Sf("K%3d%3d", i, j), tol, ele.Kpp[i][j], dnum, verb)
				}
			}
		}
	}
}

func TestingDefineDebugKbU(tst *testing.T, eid int, tmin, tmax, tol float64, verb bool) {
	//derivfcn := num.DerivFwd
	//derivfcn := num.DerivBwd
	derivfcn := num.DerivCen
	Global.DebugKb = func(d *Domain, it int) {
		elem := d.Elems[eid]
		if ele, ok := elem.(*ElemU); ok {

			// skip other times
			if tmin >= 0 && tmax >= 0 {
				if d.Sol.T < tmin || d.Sol.T > tmax {
					return
				}
			}

			// message
			//if it > 1 {
			//return
			//}
			if verb {
				io.PfYel("\nit=%2d t=%v\n", it, d.Sol.T)
			}

			// copy states and solution
			nip := len(ele.IpsElem)
			states := make([]*msolid.State, nip)
			statesBkp := make([]*msolid.State, nip)
			for i := 0; i < nip; i++ {
				states[i] = ele.States[i].GetCopy()
				statesBkp[i] = ele.StatesBkp[i].GetCopy()
			}
			Fbtmp, ΔYbkp, Yold := testing_get_aux_vectors(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					ele.States[i].Set(states[i])
					ele.StatesBkp[i].Set(statesBkp[i])
				}
				copy(d.Sol.ΔY, ΔYbkp)
			}()

			// check
			var tmp float64
			for i, I := range ele.Umap[:1] {
				for j, J := range ele.Umap[:1] {
					dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
						tmp, d.Sol.Y[J] = d.Sol.Y[J], x
						for k := 0; k < d.Ny; k++ {
							Fbtmp[k] = 0
							d.Sol.ΔY[k] = d.Sol.Y[k] - Yold[k]
						}
						if it == 0 {
							for k := 0; k < nip; k++ {
								ele.States[k].Set(states[k])
							}
						} else {
							for k := 0; k < nip; k++ {
								ele.States[k].Set(statesBkp[k])
							}
						}
						ele.Update(d.Sol)
						ele.AddToRhs(Fbtmp, d.Sol)
						res = -Fbtmp[I]
						d.Sol.Y[J] = tmp
						return res
					}, d.Sol.Y[J])
					chk.AnaNum(tst, io.Sf("K%3d%3d", i, j), tol, ele.K[i][j], dnum, verb)
				}
			}
		}
	}
}

func TestingDefineDebugKbUP(tst *testing.T, eid int, tmin, tmax, tol float64, verb bool) {
	Global.DebugKb = func(d *Domain, it int) {
		elem := d.Elems[eid]
		if e, ok := elem.(*ElemUP); ok {

			// skip other times
			if tmin >= 0 && tmax >= 0 {
				if d.Sol.T < tmin || d.Sol.T > tmax {
					return
				}
			}

			// message
			//if it > 1 {
			//return
			//}
			if verb {
				io.PfYel("\nit=%2d t=%v\n", it, d.Sol.T)
			}

			// copy states and solution
			nip := len(e.U.IpsElem)
			u_states := make([]*msolid.State, nip)
			u_statesBkp := make([]*msolid.State, nip)
			p_states := make([]*mporous.State, nip)
			p_statesBkp := make([]*mporous.State, nip)
			for i := 0; i < nip; i++ {
				u_states[i] = e.U.States[i].GetCopy()
				u_statesBkp[i] = e.U.StatesBkp[i].GetCopy()
				p_states[i] = e.P.States[i].GetCopy()
				p_statesBkp[i] = e.P.StatesBkp[i].GetCopy()
			}
			Fbtmp, ΔYbkp, Yold := testing_get_aux_vectors(d)

			// make sure to restore states and solution
			defer func() {
				for i := 0; i < nip; i++ {
					e.U.States[i].Set(u_states[i])
					e.U.StatesBkp[i].Set(u_statesBkp[i])
					e.P.States[i].Set(p_states[i])
					e.P.StatesBkp[i].Set(p_statesBkp[i])
				}
				copy(d.Sol.ΔY, ΔYbkp)
			}()

			// define restore function
			restore := func() {
				if it == 0 {
					for k := 0; k < nip; k++ {
						e.U.States[k].Set(u_states[k])
						e.P.States[k].Set(p_states[k])
					}
					return
				}
				for k := 0; k < nip; k++ {
					e.U.States[k].Set(u_statesBkp[k])
					e.P.States[k].Set(p_statesBkp[k])
				}
			}

			// check
			ni, nj := 1, 1
			testing_ana_num_K(tst, "Kuu", d, e, e.U.Umap, e.U.Umap, ni, nj, e.U.K, Fbtmp, Yold, tol, verb, restore)
			testing_ana_num_K(tst, "Kup", d, e, e.U.Umap, e.P.Pmap, ni, nj, e.Kup, Fbtmp, Yold, tol, verb, restore)
			testing_ana_num_K(tst, "Kpu", d, e, e.P.Pmap, e.U.Umap, ni, nj, e.Kpu, Fbtmp, Yold, tol, verb, restore)
			testing_ana_num_K(tst, "Kpp", d, e, e.P.Pmap, e.P.Pmap, ni, nj, e.P.Kpp, Fbtmp, Yold, tol, verb, restore)
			testing_ana_num_K(tst, "Kpf", d, e, e.P.Pmap, e.P.Fmap, ni, nj, e.P.Kpf, Fbtmp, Yold, tol, verb, restore)
			testing_ana_num_K(tst, "Kfp", d, e, e.P.Fmap, e.P.Pmap, ni, nj, e.P.Kfp, Fbtmp, Yold, tol, verb, restore)
			testing_ana_num_K(tst, "Kff", d, e, e.P.Fmap, e.P.Fmap, ni, nj, e.P.Kff, Fbtmp, Yold, tol, verb, restore)
		}
	}
}

func testing_ana_num_K(tst *testing.T, label string, d *Domain, e Elem, Imap, Jmap []int, nI, nJ int, Kana [][]float64, Fbtmp, Yold []float64, tol float64, verb bool, restore func()) {
	var imap, jmap []int
	if nI <= len(Imap) {
		imap = Imap[:nI]
	}
	if nJ <= len(Jmap) {
		jmap = Jmap[:nJ]
	}
	//derivfcn := num.DerivFwd
	//derivfcn := num.DerivBwd
	derivfcn := num.DerivCen
	var tmp float64
	for i, I := range imap {
		for j, J := range jmap {
			dnum := derivfcn(func(x float64, args ...interface{}) (res float64) {
				tmp, d.Sol.Y[J] = d.Sol.Y[J], x
				for k := 0; k < d.Ny; k++ {
					Fbtmp[k] = 0
					d.Sol.ΔY[k] = d.Sol.Y[k] - Yold[k]
				}
				restore()
				e.Update(d.Sol)
				e.AddToRhs(Fbtmp, d.Sol)
				res = -Fbtmp[I]
				d.Sol.Y[J] = tmp
				return res
			}, d.Sol.Y[J])
			chk.AnaNum(tst, io.Sf(label+"%3d%3d", i, j), tol, Kana[i][j], dnum, verb)
		}
	}
}

func testing_get_aux_vectors(d *Domain) (Fbtmp, ΔYbkp, Yold []float64) {
	Fbtmp = make([]float64, d.Ny)
	Yold = make([]float64, d.Ny)
	ΔYbkp = make([]float64, d.Ny)
	for i := 0; i < d.Ny; i++ {
		Yold[i] = d.Sol.Y[i] - d.Sol.ΔY[i]
		ΔYbkp[i] = d.Sol.ΔY[i]
	}
	return
}

func testing_up_restore_states(it, nip int, ele *ElemUP, u_states, u_statesBkp []*msolid.State, p_states, p_statesBkp []*mporous.State) {
	if it == 0 {
		for k := 0; k < nip; k++ {
			ele.U.States[k].Set(u_states[k])
			ele.P.States[k].Set(p_states[k])
		}
		return
	}
	for k := 0; k < nip; k++ {
		ele.U.States[k].Set(u_statesBkp[k])
		ele.P.States[k].Set(p_statesBkp[k])
	}
}
