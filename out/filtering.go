// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/gm"
)

// Quantity holds node or ip quantity
type Quantity struct {
	Name  string
	Value float64
}

// Locator defines interface for locating space positions
type Locator interface {
	Get(key string) []*Quantity // returns value @ X, Y, Z or Id (returs nil if not found)
}

// At implements locator at point
type At struct {
	X, Y, Z float64
}

// Along implements locator along line
type Along struct {
	A []float64 // first point on line
	B []float64 // second point on line
}

// IdsOrTags implements locator a specific nodes
type IdsOrTags []int

func (o At) Get(key string) []*Quantity {

	// coordinates
	c := []float64{o.X, o.Y, o.Z}

	// node quantity
	entry := NodBins.Find(c)
	if entry != nil {
		q := get_nod_quantity(key, entry)
		if q != nil {
			return []*Quantity{q}
		}
	}

	// integration point quantity
	entry = IpsBins.Find(c)
	if entry != nil {
		q := get_ip_quantity(key, entry)
		if q != nil {
			return []*Quantity{q}
		}
	}
	return nil
}

func (o Along) Get(key string) []*Quantity {
	return nil
}
func (o IdsOrTags) Get(key string) []*Quantity {
	return nil
}

func (o IdsOrTags) GetNodes() (nodes []*fem.Node) {
	for _, idortag := range o {
		if idortag < 0 {
			tag := idortag
			verts := Dom.Msh.VertTag2verts[tag]
			for _, v := range verts {
				nodes = append(nodes, Dom.Nodes[v.Id])
			}
		} else {
			id := idortag
			nodes = append(nodes, Dom.Nodes[id])
		}
	}
	return
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

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

func get_nod_quantity(key string, entry *gm.BinEntry) *Quantity {
	if entry != nil {
		key, ntderiv := parse_key(key)
		nod := Dom.Nodes[entry.Id]
		dof := nod.GetDof(key)
		if dof != nil {
			switch ntderiv {
			case 0:
				return &Quantity{key, Dom.Sol.Y[dof.Eq]}
			case 1:
				return &Quantity{key, Dom.Sol.Dydt[dof.Eq]}
			case 2:
				return &Quantity{key, Dom.Sol.D2ydt2[dof.Eq]}
			}
		}
	}
	return nil
}

func get_ip_quantity(key string, entry *gm.BinEntry) *Quantity {
	if entry != nil {
		o := Ipoints[entry.Id]
		if o.P != nil {
			switch key {
			case "pl":
				return &Quantity{key, o.P.Pl}
			case "sl":
				return &Quantity{key, o.P.Sl}
			}
		}
		if o.U != nil {
			switch key {
			case "sx":
				return &Quantity{key, o.U.Sig[0]}
			case "sy":
				return &Quantity{key, o.U.Sig[1]}
			case "sz":
				return &Quantity{key, o.U.Sig[2]}
			}
		}
	}
	return nil
}
