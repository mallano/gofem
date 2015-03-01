// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "math"

// Point holds information about one specific node xor one integration point
type Point struct {
	Vid  int                  // vertex id; -1 => not node
	IpId int                  // ip id; -1 => not ip
	X    []float64            // coordinates
	Dist float64              // distance from reference point (if along line)
	Vals map[string][]float64 // combines keys to results along time; e.g. "ux" : {0,0.0001,0,0002,...}
}

// Points is a set of Point
type Points []*Point

// Len the length of Points
func (o Points) Len() int {
	return len(o)
}

// Swap swaps two Points
func (o Points) Swap(i, j int) {
	o[i], o[j] = o[j], o[i]
}

// Less compares Points considering Dist
func (o Points) Less(i, j int) bool {
	return o[i].Dist < o[j].Dist
}

func get_nod_point(vid int, A []float64) *Point {
	nod := Dom.Vid2node[vid]
	if nod != nil {
		var dist float64
		if A != nil {
			dist = dist_point_point(nod.Vert.C, A)
		}
		return &Point{vid, -1, nod.Vert.C, dist, make(map[string][]float64)}
	}
	return nil
}

func get_ip_point(ipid int, A []float64) *Point {
	ip := Ipoints[ipid]
	if ip != nil {
		var dist float64
		if A != nil {
			dist = dist_point_point(ip.X, A)
		}
		return &Point{-1, ipid, ip.X, dist, make(map[string][]float64)}
	}
	return nil
}

// dist_point_point returns the distance between two points
func dist_point_point(a, b []float64) float64 {
	if len(a) == 2 {
		return math.Sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
	}
	return math.Sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]))
}
