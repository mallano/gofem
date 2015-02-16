// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

// PointLocator defines interface for locating space positions
type PointLocator interface {
	AtPoint(key string) Quantities
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
	id := NodBins.Find(o)
	if id >= 0 {
		q := get_nod_quantity(key, id, 0)
		if q != nil {
			return Quantities{q}
		}
	}

	// integration point quantity
	id = IpsBins.Find(o)
	if id >= 0 {
		q := get_ip_quantity(key, id, 0)
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
