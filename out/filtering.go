// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

// PointLocator defines interface for locating space positions
type Locator interface {
	Locate() Points
}

// At implements locator at point => PointLocator
type At []float64

// N implements node locator
// Ids or tags of vertices can be stored in Verts
type N []int

// P implements [element][integrationPoint] locator
// Pairs of ids or tags of cells and integration points indices can be stored in Cells
type P [][]int

// Along implements locator along line => LineLocator
//  Example: with 2 points in 3D: {{0,0,0}, {1,1,1}}
type Along [][]float64

// AlongX implements LineLocator
type AlongX []float64

// AlongY implements LineLocator
type AlongY []float64

// AlongZ implements LineLocator
type AlongZ []float64

// AtPoint returns quantity at point
func (o At) Locate() Points {

	// node quantity
	id := NodBins.Find(o)
	if id >= 0 {
		q := get_nod_point(id, 0)
		if q != nil {
			return Points{q}
		}
	}

	// integration point quantity
	id = IpsBins.Find(o)
	if id >= 0 {
		q := get_ip_point(id, 0)
		if q != nil {
			return Points{q}
		}
	}
	return nil
}

// AtPoint returns quantity at point
func (o N) Locate() (res Points) {
	for _, idortag := range o {
		if idortag < 0 {
			tag := idortag
			verts := Dom.Msh.VertTag2verts[tag]
			for _, v := range verts {
				q := get_nod_point(v.Id, 0)
				if q != nil {
					res = append(res, q)
				}
			}
		} else {
			id := idortag
			q := get_nod_point(id, 0)
			if q != nil {
				res = append(res, q)
			}
		}
	}
	return
}

// AtPoint returns quantity at point
func (o P) Locate() (res Points) {
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
				q := get_ip_point(ipid, 0)
				if q != nil {
					res = append(res, q)
				}
			}
		} else {
			cid := idortag
			idx := o[i][1]
			ipid := Cid2ips[cid][idx]
			q := get_ip_point(ipid, 0)
			if q != nil {
				res = append(res, q)
			}
		}
	}
	return
}

// Along returns quantity along line
func (o Along) Locate() (res Points) {

	// check if there are two points
	if len(o) != 2 {
		return
	}
	A := o[0]
	B := o[1]

	// node quantities
	ids := NodBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		q := get_nod_point(id, 0)
		if q != nil {
			res = append(res, q)
		}
	}

	// integration point quantitites
	ids = IpsBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		q := get_ip_point(id, 0)
		if q != nil {
			res = append(res, q)
		}
	}
	return
}

// AllIps returns all cell/ip indices
func AllIps() P {
	p := make([][]int, 0)
	for i, ips := range Cid2ips {
		for j, _ := range ips {
			p = append(p, []int{i, j})
		}
	}
	return p
}
