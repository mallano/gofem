// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

// PointLocator defines interface for locating space positions
type PointLocator interface {
	AtPoint(key string) Quantities // returns value @ X, Y, Z or Id (returs nil if not found)
}

// LineLocator defines interface for locating space positions
type LineLocator interface {
	AlongLine(key string) Quantities
}

// At implements locator at point => PointLocator
type At []float64

// IdsOrTags implements locator a specific nodes => PointLocator
type IdsOrTags []int

// Along implements locator along line => LineLocator
type Along struct {
	A []float64 // first point on line
	B []float64 // second point on line
}

// AlongX implements LineLocator
type AlongX []float64

// AlongY implements LineLocator
type AlongY []float64

// AlongZ implements LineLocator
type AlongZ []float64

// AtPoint returns quantity at point
func (o At) AtPoint(key string) Quantities {

	// node quantity
	entry := NodBins.Find(o)
	if entry != nil {
		q := get_nod_quantity(key, entry.Id, 0)
		if q != nil {
			return Quantities{q}
		}
	}

	// integration point quantity
	entry = IpsBins.Find(o)
	if entry != nil {
		q := get_ip_quantity(key, entry.Id, 0)
		if q != nil {
			return Quantities{q}
		}
	}
	return nil
}

// AtPoint returns quantity at point
func (o IdsOrTags) AtPoint(key string) (res Quantities) {
	for _, idortag := range o {
		if idortag < 0 {
			tag := idortag
			verts := Dom.Msh.VertTag2verts[tag]
			for _, v := range verts {
				q := get_nod_quantity(key, v.Id, 0)
				if q != nil {
					res = append(res, q)
				}
			}
		} else {
			id := idortag
			q := get_nod_quantity(key, id, 0)
			if q != nil {
				res = append(res, q)
			}
		}
	}
	return
}

// Along returns quantity along line
func (o Along) AlongLine(key string) Quantities {
	return nil
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

func get_nod_quantity(key string, nid int, dist float64) *Quantity {
	key, ntderiv := parse_key(key)
	nod := Dom.Nodes[nid]
	dof := nod.GetDof(key)
	if dof != nil {
		switch ntderiv {
		case 0:
			return &Quantity{Dom.Sol.Y[dof.Eq], dist}
		case 1:
			return &Quantity{Dom.Sol.Dydt[dof.Eq], dist}
		case 2:
			return &Quantity{Dom.Sol.D2ydt2[dof.Eq], dist}
		}
	}
	return nil
}

func get_ip_quantity(key string, ipid int, dist float64) *Quantity {
	o := Ipoints[ipid]
	if o.P != nil {
		switch key {
		case "pl":
			return &Quantity{o.P.Pl, dist}
		case "sl":
			return &Quantity{o.P.Sl, dist}
		}
	}
	if o.U != nil {
		switch key {
		case "sx":
			return &Quantity{o.U.Sig[0], dist}
		case "sy":
			return &Quantity{o.U.Sig[1], dist}
		case "sz":
			return &Quantity{o.U.Sig[2], dist}
		}
	}
	return nil
}
