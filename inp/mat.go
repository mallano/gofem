// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"log"

	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// Material holds material data
type Material struct {
	Name  string   `json:"name"`  // name of material
	Desc  string   `json:"desc"`  // description of material
	Model string   `json:"model"` // name of model ex: 'dp', 'vm', 'elast', etc.
	Extra string   `json:"extra"` // extra information about this material
	Prms  fun.Prms `json:"prms"`  // prms holds all model parameters for this material
}

// Mats holds materials
type MatsData []*Material

// MatDb implements a database of materials
type MatDb struct {
	Functions FuncsData `json:"functions"` // all functions
	Materials MatsData  `json:"materials"` // all materials
}

// ReadMat reads all materials data from a .mat JSON file
//  Note: returns nil on errors
func ReadMat(fn string) *MatDb {

	// new mat
	var o MatDb

	// read file
	b, err := utl.ReadFile(fn)
	if LogErr(err, "mat: cannot open materials file "+fn+"\n") {
		return nil
	}

	// decode
	if LogErr(json.Unmarshal(b, &o), "mat: cannot unmarshal materials file "+fn+"\n") {
		return nil
	}

	// log
	log.Printf("mat: fn=%s nfunctions=%d nmaterials=%d\n", fn, len(o.Functions), len(o.Materials))
	return &o
}

// Get returns a material
//  Note: returns nil if not found
func (o *MatDb) Get(name string) *Material {
	for _, mat := range o.Materials {
		if mat.Name == name {
			return mat
		}
	}
	return nil
}

// String prints one function
func (o *Material) String() string {
	fun.G_extraindent = "  "
	fun.G_openbrackets = false
	return utl.Sf("    {\n      \"name\"  : %q,\n      \"desc\"  : %q,\n      \"model\" : %q,\n      \"prms\"  : [\n%v\n    }", o.Name, o.Desc, o.Model, o.Prms)
}

// String prints materials
func (o MatsData) String() string {
	l := "  \"materials\" : [\n"
	for i, m := range o {
		if i > 0 {
			l += ",\n"
		}
		l += utl.Sf("%v", m)
	}
	l += "\n  ]"
	return l
}

// String outputs all materials
func (o MatDb) String() string {
	return utl.Sf("{\n%v,\n%v\n}", o.Functions, o.Materials)
}

// MatfileOld2New converts an old mat file to new mat file format
//  convertsymbols -- convert symbols with Greek characters to ANSI
func MatfileOld2New(dirout string, fnnew, fnold string, convertsymbols bool) {

	// oldmatdata implements the old data format
	type oldmatdata struct {
		Name  string
		Desc  string
		Model string
		Prms  []string
		Vals  []float64
		Units []string
		Extra string
	}

	// data holder
	var mats_old []oldmatdata

	// read file
	b, err := utl.ReadFile(fnold)
	if err != nil {
		utl.PfRed("cannot open file: %v", err.Error())
		return
	}

	// decode
	err = json.Unmarshal(b, &mats_old)
	if err != nil {
		utl.PfRed("cannot unmarshal file: %v", err.Error())
		return
	}

	// new data holder
	var mats_new MatDb
	for _, m := range mats_old {
		prms := []*fun.Prm{}
		has_units := len(m.Units) == len(m.Prms)
		for i, prm := range m.Prms {
			name := prm
			if convertsymbols {
				if n, ok := conversiontable[prm]; ok {
					name = n
				}
			}
			if has_units {
				prms = append(prms, &fun.Prm{N: name, V: m.Vals[i], U: m.Units[i]})
			} else {
				prms = append(prms, &fun.Prm{N: name, V: m.Vals[i]})
			}
		}
		mat_new := Material{
			Name:  m.Name,
			Desc:  m.Desc,
			Model: m.Model,
			Extra: m.Extra,
			Prms:  prms,
		}
		mats_new.Materials = append(mats_new.Materials, &mat_new)
	}

	// save file
	if dirout == "" {
		utl.WriteFileS(fnnew, mats_new.String())
		return
	}
	utl.WriteFileSD(dirout, fnnew, mats_new.String())
}

// convert greek to ansi
var conversiontable = map[string]string{
	"α":   "alp",
	"αl":  "alpL",
	"β":   "bet",
	"βl":  "betaL",
	"βg":  "betaG",
	"βd":  "betD",
	"βw":  "betW",
	"β1":  "bet1",
	"β2":  "bet2",
	"ν":   "nu",
	"φ":   "phi",
	"λ":   "lam",
	"λd":  "lamD",
	"λw":  "lamW",
	"λ0l": "lam0L",
	"λ1l": "lam1L",
	"λ0g": "lam0G",
	"λ1g": "lam1G",
	"ρL":  "RhoL",
	"ρG":  "RhoG",
	"ρS":  "RhoS",
	"RΘg": "RthG",
}
