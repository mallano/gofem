// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "sort"

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

// Locate finds points
func (o At) Locate() Points {

	// node
	id := NodBins.Find(o)
	if id >= 0 {
		q := get_nod_point(id, nil)
		if q != nil {
			return Points{q}
		}
	}

	// integration point
	id = IpsBins.Find(o)
	if id >= 0 {
		q := get_ip_point(id, nil)
		if q != nil {
			return Points{q}
		}
	}
	return nil
}

// Locate finds points
func (o N) Locate() (res Points) {
	var A []float64 // reference point
	for _, idortag := range o {
		if idortag < 0 {
			tag := idortag
			verts := Dom.Msh.VertTag2verts[tag]
			for _, v := range verts {
				q := get_nod_point(v.Id, A)
				if q != nil {
					res = append(res, q)
					if A == nil {
						A = q.X
					}
				}
			}
		} else {
			id := idortag
			q := get_nod_point(id, A)
			if q != nil {
				res = append(res, q)
				if A == nil {
					A = q.X
				}
			}
		}
	}
	return
}

// Locate finds points
func (o P) Locate() (res Points) {
	var A []float64 // reference point
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
				q := get_ip_point(ipid, A)
				if q != nil {
					res = append(res, q)
					if A == nil {
						A = q.X
					}
				}
			}
		} else {
			cid := idortag
			idx := o[i][1]
			ipid := Cid2ips[cid][idx]
			q := get_ip_point(ipid, A)
			if q != nil {
				res = append(res, q)
				if A == nil {
					A = q.X
				}
			}
		}
	}
	return
}

// Along finds points
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
		q := get_nod_point(id, A)
		if q != nil {
			res = append(res, q)
		}
	}

	// integration point quantitites
	ids = IpsBins.FindAlongLine(A, B, TolC)
	for _, id := range ids {
		q := get_ip_point(id, A)
		if q != nil {
			res = append(res, q)
		}
	}
	sort.Sort(res)
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
