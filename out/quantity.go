// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

// Quantity holds node or ip quantity
type Quantity struct {
	Value float64 // value of, e.g., "ux"
	Dist  float64 // distance from reference point (if along line)
}

// Quantities is a set of Quantity
type Quantities []*Quantity

// Len the length of Quantities
func (o Quantities) Len() int {
	return len(o)
}

// Swap swaps two Quantities
func (o Quantities) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Quantities considering Dist
func (o Quantities) Less(i, j int) bool {
	return o[i].Dist < o[j].Dist
}
