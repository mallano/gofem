// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/utl"
)

// FuncData holds function definition
type FuncData struct {
	Name string   `json:"name"` // name of function. ex: zero, load, myfunction1, etc.
	Type string   `json:"type"` // type of function. ex: cte, rmp
	Prms fun.Prms `json:"prms"` // parameters
}

// Funcs holds functions
type FuncsData []*FuncData

// GetOrPanic returns function or panic
func (o FuncsData) GetOrPanic(name string) fun.Func {
	if name == "zero" {
		return &fun.Zero
	}
	for _, f := range o {
		if f.Name == name {
			return fun.New(f.Type, f.Prms)
		}
	}
	utl.Panic("cannot find function named %s", name)
	return nil
}

// String prints one function
func (o FuncData) String() string {
	fun.G_extraindent = "  "
	fun.G_openbrackets = false
	return utl.Sf("    {\n      \"name\":%q, \"type\":%q, \"prms\" : [\n%v\n    }", o.Name, o.Type, o.Prms)
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
		l += utl.Sf("%v", f)
	}
	l += "\n  ]"
	return l
}
