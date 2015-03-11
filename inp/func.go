// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"github.com/cpmech/gosl/fun"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// PlotFdata holds information to plot functions
type PlotFdata struct {
	Ti    float64  `json:"ti"`    // initial time
	Tf    float64  `json:"tf"`    // final time
	Np    int      `json:"np"`    // number of points
	Skip  []string `json:"skip"`  // skip functions
	WithG bool     `json:"withg"` // with dF/dt
	WithH bool     `json:"withh"` // with d²F/dt²
	Eps   bool     `json:"eps"`   // save eps instead of png
}

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

// PlotAll plot all functions
func (o FuncsData) PlotAll(pd *PlotFdata, dirout, fnkey string) {
	ext := "png"
	if pd.Eps {
		ext = "eps"
	}
	for _, f := range o {
		if utl.StrIndexSmall(pd.Skip, f.Name) >= 0 {
			continue
		}
		fn := io.Sf("fcn-%s-%s.%s", fnkey, f.Name, ext)
		ff := o.Get(f.Name)
		if ff != nil {
			fun.PlotT(ff, dirout, fn, pd.Ti, pd.Tf, nil, pd.Np, "", pd.WithG, pd.WithH, true, false, nil)
		}
	}
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////////

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
