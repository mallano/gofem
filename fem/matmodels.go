// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/mconduct"
	"github.com/cpmech/gofem/mporous"
	"github.com/cpmech/gofem/mreten"
	"github.com/cpmech/gofem/msolid"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// GetAndInitPorousModel get porous model from material name
// It returns nil on errors, after logging
func GetAndInitPorousModel(matname string) *mporous.Model {

	// materials
	cndmat, lrmmat, pormat, err := Global.Sim.Mdb.GroupGet3(matname, "c", "l", "p")
	if LogErr(err, io.Sf("materials database failed on getting %q (porous) group\n", matname)) {
		return nil
	}

	// conductivity models
	simfnk := Global.Sim.Data.FnameKey
	getnew := false
	cnd := mconduct.GetModel(simfnk, cndmat.Name, cndmat.Model, getnew)
	if LogErrCond(cnd == nil, "cannot allocate conductivity models with name=%q", cndmat.Model) {
		return nil
	}

	// retention model
	lrm := mreten.GetModel(simfnk, lrmmat.Name, lrmmat.Model, getnew)
	if LogErrCond(lrm == nil, "cannot allocate liquid retention model with name=%q", lrmmat.Model) {
		return nil
	}

	// porous model
	mdl := mporous.GetModel(simfnk, pormat.Name, getnew)
	if LogErrCond(mdl == nil, "cannot allocate model for porous medium with name=%q", pormat.Name) {
		return nil
	}

	// initialise all models
	if LogErr(cnd.Init(cndmat.Prms), "cannot initialise conductivity model") {
		return nil
	}
	if LogErr(lrm.Init(lrmmat.Prms), "cannot initialise liquid retention model") {
		return nil
	}
	if LogErr(mdl.Init(pormat.Prms, cnd, lrm), "cannot initialise porous model") {
		return nil
	}

	// results
	return mdl
}

func GetAndInitSolidModel(matname string, ndim int) (msolid.Model, fun.Prms) {

	// material name
	matdata := Global.Sim.Mdb.Get(matname)
	if LogErrCond(matdata == nil, "materials database failed on getting %q (solid) material\n", matname) {
		return nil, nil
	}
	mdlname := matdata.Model

	// handle groups
	if mdlname == "group" {
		if s_matname, found := io.Keycode(matdata.Extra, "s"); found {
			matname = s_matname
			matdata = Global.Sim.Mdb.Get(matname)
			if LogErrCond(matdata == nil, "materials database failed on getting %q (solid/sub) material\n", matname) {
				return nil, nil
			}
			mdlname = matdata.Model
		} else {
			LogErrCond(true, "cannot find solid model in grouped material data. 's' subkey needed in Extra field")
			return nil, nil
		}
	}

	// initialise model
	mdl := msolid.GetModel(Global.Sim.Data.FnameKey, matname, mdlname, false)
	if LogErrCond(mdl == nil, "cannot find solid model named %q", mdlname) {
		return nil, nil
	}
	if LogErr(mdl.Init(ndim, Global.Sim.Data.Pstress, matdata.Prms), "solid model initialisation failed") {
		return nil, nil
	}

	// results
	return mdl, matdata.Prms
}
