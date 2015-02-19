// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "github.com/cpmech/gosl/io"

// Quantity holds node or ip quantity
type Quantity struct {
	Value *float64  // value of, e.g., "ux"
	Dist  float64   // distance from reference point (if along line)
	X     []float64 // coordinates
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

// parse_key parses key like "duxdt" returning "ux" and time derivative number
//  Output: {key, number-of-time-derivatives}
//  Examples:  "ux"      => "ux", 0
//             "duxdt"   => "ux", 1
//             "d2uxdt2" => "ux", 2
func parse_key(key string) (string, int) {
	if len(key) > 3 {
		n := len(key)
		if key[:1] == "d" && key[n-2:] == "dt" {
			return key[1 : n-2], 1
		}
		if len(key) > 5 {
			if key[:2] == "d2" && key[n-3:] == "dt2" {
				return key[2 : n-3], 2
			}
		}
	}
	return key, 0
}

func get_nod_quantity(key string, nid int, dist float64) *Quantity {
	key, ntderiv := parse_key(key)
	nod := Dom.Nodes[nid]
	dof := nod.GetDof(key)
	if dof != nil {
		switch ntderiv {
		case 0:
			return &Quantity{&Dom.Sol.Y[dof.Eq], dist, nod.Vert.C}
		case 1:
			return &Quantity{&Dom.Sol.Dydt[dof.Eq], dist, nod.Vert.C}
		case 2:
			return &Quantity{&Dom.Sol.D2ydt2[dof.Eq], dist, nod.Vert.C}
		}
	}
	return nil
}

func get_ip_quantity(key string, ipid int, dist float64) *Quantity {
	ip := Ipoints[ipid]
	if ip.P != nil {
		switch key {
		case "pl":
			return &Quantity{&ip.P.Pl, dist, ip.X}
		case "sl":
			return &Quantity{&ip.P.Sl, dist, ip.X}
		}
	}
	if ip.U != nil {
		switch key {
		case "sx":
			return &Quantity{&ip.U.Sig[0], dist, ip.X}
		case "sy":
			return &Quantity{&ip.U.Sig[1], dist, ip.X}
		case "sz":
			return &Quantity{&ip.U.Sig[2], dist, ip.X}
		case "sxy":
			return &Quantity{&ip.U.Sig[3], dist, ip.X}
		case "syz":
			if len(ip.U.Sig) == 6 {
				return &Quantity{&ip.U.Sig[4], dist, ip.X}
			}
		case "szx":
			if len(ip.U.Sig) == 6 {
				return &Quantity{&ip.U.Sig[5], dist, ip.X}
			}
		}
	}
	return nil
}

func (o Quantity) String() string {
	l := io.Sf("{\"value\":%g, \"dist\":%g, \"x\":[%g,%g", *o.Value, o.Dist, o.X[0], o.X[1])
	if len(o.X) == 3 {
		l += io.Sf(",%g", o.X[2])
	}
	l += "]}"
	return l
}

func (o Quantities) String() string {
	l := "{\n"
	for i, q := range o {
		if i > 0 {
			l += ",\n"
		}
		l += io.Sf("  %v", q)
	}
	l += "\n}"
	return l
}
