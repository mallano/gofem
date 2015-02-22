// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"log"
	"path/filepath"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
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
func ReadMat(dir, fn string) *MatDb {

	// new mat
	var o MatDb

	// read file
	b, err := io.ReadFile(filepath.Join(dir, fn))
	if LogErr(err, "mat: cannot open materials file "+dir+"/"+fn) {
		return nil
	}

	// decode
	if LogErr(json.Unmarshal(b, &o), "mat: cannot unmarshal materials file "+fn) {
		return nil
	}

	// log
	log.Printf("mat: fn=%s nfunctions=%d nmaterials=%d\n", fn, len(o.Functions), len(o.Materials))
	return &o
}

// Get returns a material
//  Note: returns nil if not found
func (o MatDb) Get(name string) *Material {
	for _, mat := range o.Materials {
		if mat.Name == name {
			return mat
		}
	}
	return nil
}

// GroupGet parses group data
//  Note: returns nil on failure
func (o MatDb) GroupGet(matname, key string) *Material {
	mat := o.Get(matname)
	if mat == nil {
		return nil
	}
	if submatname, found := io.Keycode(mat.Extra, key); found {
		return o.Get(submatname)
	}
	return nil
}

// GroupGet3 parses group data
func (o MatDb) GroupGet3(matname, key1, key2, key3 string) (m1, m2, m3 *Material, err error) {
	mat := o.Get(matname)
	if mat == nil {
		err = chk.Err("cannot find material named %q", matname)
		return
	}
	if submat1, found := io.Keycode(mat.Extra, key1); found {
		m1 = o.Get(submat1)
	} else {
		err = chk.Err("cannot find key %q in grouped material data %q", key1, mat.Extra)
		return
	}
	if submat2, found := io.Keycode(mat.Extra, key2); found {
		m2 = o.Get(submat2)
	} else {
		err = chk.Err("cannot find key %q in grouped material data %q", key2, mat.Extra)
		return
	}
	if submat3, found := io.Keycode(mat.Extra, key3); found {
		m3 = o.Get(submat3)
	} else {
		err = chk.Err("cannot find key %q in grouped material data %q", key3, mat.Extra)
		return
	}
	if m1 == nil || m2 == nil || m3 == nil {
		err = chk.Err("material data in grouped materials cannot be parsed")
	}
	return
}

// String prints one function
func (o *Material) String() string {
	fun.G_extraindent = "  "
	fun.G_openbrackets = false
	return io.Sf("    {\n      \"name\"  : %q,\n      \"desc\"  : %q,\n      \"model\" : %q,\n      \"prms\"  : [\n%v\n    }", o.Name, o.Desc, o.Model, o.Prms)
}

// String prints materials
func (o MatsData) String() string {
	l := "  \"materials\" : [\n"
	for i, m := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("%v", m)
	}
	l += "\n  ]"
	return l
}

// String outputs all materials
func (o MatDb) String() string {
	return io.Sf("{\n%v,\n%v\n}", o.Functions, o.Materials)
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
	b, err := io.ReadFile(fnold)
	if err != nil {
		io.PfRed("cannot open file: %v", err.Error())
		return
	}

	// decode
	err = json.Unmarshal(b, &mats_old)
	if err != nil {
		io.PfRed("cannot unmarshal file: %v", err.Error())
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
		io.WriteFileS(fnnew, mats_new.String())
		return
	}
	io.WriteFileSD(dirout, fnnew, mats_new.String())
}

// MatfileNew2Old converts a new mat file to the old mat file format
//  convertsymbols -- convert back symbols with Greek characters to UTF-8
func MatfileNew2Old(dirout string, fnold, fnnew string, convertsymbols bool) {

	// read file
	b, err := io.ReadFile(fnnew)
	if err != nil {
		io.PfRed("cannot open file: %v", err.Error())
		return
	}

	// decode
	var mats_new MatDb
	err = json.Unmarshal(b, &mats_new)
	if err != nil {
		io.PfRed("cannot unmarshal file: %v", err.Error())
		return
	}

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

	// convert
	var mats_old []oldmatdata
	for _, m := range mats_new.Materials {
		var oldmat oldmatdata
		oldmat.Name = m.Name
		oldmat.Desc = m.Desc
		oldmat.Model = m.Model
		oldmat.Extra = m.Extra
		for _, prm := range m.Prms {
			name := prm.N
			if convertsymbols {
				if n, ok := invertedconversiontable[prm.N]; ok {
					name = n
				}
			}
			oldmat.Prms = append(oldmat.Prms, name)
			oldmat.Vals = append(oldmat.Vals, prm.V)
			if len(prm.U) > 0 {
				oldmat.Units = append(oldmat.Units, prm.U)
			}
		}
		mats_old = append(mats_old, oldmat)
	}

	// encode
	buf, err := json.MarshalIndent(mats_old, "", "  ")
	if err != nil {
		return
	}

	// save file
	if dirout == "" {
		io.WriteFileS(fnold, string(buf))
		return
	}
	io.WriteFileSD(dirout, fnold, string(buf))
}

// convert greek to ansi
var conversiontable = map[string]string{
	"α":   "alp",
	"αl":  "alpl",
	"β":   "bet",
	"βl":  "betl",
	"βg":  "betg",
	"βd":  "betd",
	"βw":  "betw",
	"β1":  "bet1",
	"β2":  "bet2",
	"μ":   "mu",
	"ν":   "nu",
	"φ":   "phi",
	"λ":   "lam",
	"λd":  "lamd",
	"λw":  "lamw",
	"λ0l": "lam0l",
	"λ1l": "lam1l",
	"λ0g": "lam0g",
	"λ1g": "lam1g",
	"ρL":  "RhoL",
	"ρG":  "RhoG",
	"ρ":   "rho",
	"ρS":  "rho",
	"RΘg": "RthG",
	"τy0": "tauy0",
}

var invertedconversiontable = map[string]string{
	"alp":   "α",
	"alpL":  "αl",
	"bet":   "β",
	"betaL": "βl",
	"betaG": "βg",
	"betD":  "βd",
	"betW":  "βw",
	"bet1":  "β1",
	"bet2":  "β2",
	"nu":    "ν",
	"phi":   "φ",
	"lam":   "λ",
	"lamD":  "λd",
	"lamW":  "λw",
	"lam0L": "λ0l",
	"lam1L": "λ1l",
	"lam0G": "λ0g",
	"lam1G": "λ1g",
	"Rho":   "ρ",
	"RhoL":  "ρL",
	"RhoG":  "ρG",
	"RhoS":  "ρS",
	"RthG":  "RΘg",
}
