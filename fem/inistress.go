// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/utl"
)

// SetIniStress sets the initial state with initial stresses
func (o *Domain) SetIniStress(stg *inp.Stage) (ok bool) {

	// check for hydrost data
	dat := stg.IniStress
	if dat == nil {
		return true
	}

	// set elements with homogeneous stress state
	if dat.Hom {

		// isotropic state
		if dat.Iso {
			for _, e := range o.ElemIntvars {

				// get element's integration points data
				ele := e.(Elem)
				_, d := ele.OutIpsData()
				nip := len(d)

				// build map with isotropic and homogeneus state
				v := utl.DblVals(nip, dat.S0)
				ivs := map[string][]float64{"sx": v, "sy": v, "sz": v}

				// set element's states
				if LogErrCond(!e.SetIvs(ivs), "homogeneous/isotropic: element's internal values setting failed") {
					return
				}
			}
			log.Printf("dom: initial homogeneous/isotropic state set with σ0 = %g", dat.S0)
			return true
		}

		// plane-strain state
		if dat.Psa {
			sz := dat.Nu * (dat.Sh + dat.Sv)
			for _, e := range o.ElemIntvars {

				// get element's integration points data
				ele := e.(Elem)
				_, d := ele.OutIpsData()
				nip := len(d)

				// build map with plane-strain and homogeneus state
				vx := utl.DblVals(nip, dat.Sh)
				vy := utl.DblVals(nip, dat.Sv)
				vz := utl.DblVals(nip, sz)
				ivs := map[string][]float64{"sx": vx, "sy": vy, "sz": vz}

				// set element's states
				if LogErrCond(!e.SetIvs(ivs), "homogeneous/plane-strain: element's internal values setting failed") {
					return
				}
			}
			log.Printf("dom: initial homogeneous/plane-strain state set with sx=%g sy=%g sz=%g", dat.Sh, dat.Sv, sz)
			return true
		}
	}
	return true
}
