// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/io"
)

// Layer holds information of one soil layer
type Layer struct {
	Tags  []int   // tags of cells within this layer
	Zmin  float64 // coordinate (elevation) at bottom of layer
	Zmax  float64 // coordinate (elevation) at top of layer
	DsigV float64 // absolute value of vertical stress incremented added by this layer
}

// StressInc1 returns the (absolute) increment of vertical stress from z to the top of this layer.
func (o Layer) StressInc(z float64) (ΔσVabs float64) {
	return
}

func (o Layer) State(σ0abs, z float64) (sx, sy, sz, pl, pg float64, err error) {
	return
}

// Layers is a set of Layer
type Layers []*Layer

// Len the length of Layers
func (o Layers) Len() int {
	return len(o)
}

// Swap swaps two Layers
func (o Layers) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Layers: sort from top to bottom
func (o Layers) Less(i, j int) bool {
	return o[i].Zmin > o[j].Zmin
}

// SetGeoSt sets the initial state to a hydrostatic condition
func (o *Domain) SetGeoSt(stg *inp.Stage) (ok bool) {
	return
}

// String prints a json formatted string with Layers' content
func (o Layers) String() string {
	if len(o) == 0 {
		return "[]"
	}
	l := "[\n"
	for i, lay := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("  { \"Tags\":%v, \"Zmin\":%g, \"Zmax\":%g, \"DsigV\":%g }", lay.Tags, lay.Zmin, lay.Zmax, lay.DsigV)
	}
	l += "\n]"
	return l
}
