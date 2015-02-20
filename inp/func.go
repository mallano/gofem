// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
)

// FuncData holds function definition
type FuncData struct {
	Name string   `json:"name"` // name of function. ex: zero, load, myfunction1, etc.
	Type string   `json:"type"` // type of function. ex: cte, rmp
	Prms fun.Prms `json:"prms"` // parameters
}

// Funcs holds functions
type FuncsData []*FuncData

// Get returns function by name
//  Note: returns nil if not found
func (o FuncsData) Get(name string) fun.Func {
	if name == "zero" {
		return &fun.Zero
	}
	for _, f := range o {
		if f.Name == name {
			fcn, err := fun.New(f.Type, f.Prms)
			if LogErr(err, "FuncsData.Get") {
				return nil
			}
			return fcn
		}
	}
	return nil
}

// String prints one function
func (o FuncData) String() string {
	fun.G_extraindent = "  "
	fun.G_openbrackets = false
	return io.Sf("    {\n      \"name\":%q, \"type\":%q, \"prms\" : [\n%v\n    }", o.Name, o.Type, o.Prms)
}

// String prints functions
func (o FuncsData) String() string {
	if len(o) == 0 {
		return "  \"functions\" : []"
	}
	l := "  \"functions\" : [\n"
	for i, f := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("%v", f)
	}
	l += "\n  ]"
	return l
}
