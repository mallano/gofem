// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

// PointLocator defines interface for locating space positions
type Locator interface {
	Locate(key string) Quantities
}

// At implements locator at point => PointLocator
type At []float64

// Verts implements locator at point => PointLocator
// Ids or tags of vertices can be stored in Verts
type Verts []int

// Cells implements locator at point => PointLocator
// Pairs of ids or tags of cells and integration points indices can be stored in Cells
type Cells [][]int

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
func (o At) Locate(key string) Quantities {

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
func (o Verts) Locate(key string) (res Quantities) {
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

// AtPoint returns quantity at point
func (o Cells) Locate(key string) (res Quantities) {
	ncells := len(o)
	for i := 0; i < ncells; i++ {
		if len(o[i]) != 2 {
			continue
		}
		idortag := o[i][0]
		if idortag < 0 {
			tag := idortag
			cells := Dom.Msh.CellTag2cells[tag]
			for _, c := range cells {
				cid := c.Id
				idx := o[i][1]
				ipid := Cid2ips[cid][idx]
				q := get_ip_quantity(key, ipid, 0)
				if q != nil {
					res = append(res, q)
				}
			}
		} else {
			cid := idortag
			idx := o[i][1]
			ipid := Cid2ips[cid][idx]
			q := get_ip_quantity(key, ipid, 0)
			if q != nil {
				res = append(res, q)
			}
		}
	}
	return
}

// Along returns quantity along line
func (o Along) Locate(key string) (res Quantities) {

	// node quantities
	ids := NodBins.FindAlongLine(o.A, o.B, TolC)
	for _, id := range ids {
		q := get_nod_quantity(key, id, 0)
		if q != nil {
			res = append(res, q)
		}
	}

	// integration point quantitites
	ids = IpsBins.FindAlongLine(o.A, o.B, TolC)
	for _, id := range ids {
		q := get_ip_quantity(key, id, 0)
		if q != nil {
			res = append(res, q)
		}
	}
	return
}

// AllCells returns all cell/ip indices
func AllCells() Cells {
	cells := make([][]int, 0)
	for i, ips := range Cid2ips {
		for j, _ := range ips {
			cells = append(cells, []int{i, j})
		}
	}
	return cells
}
