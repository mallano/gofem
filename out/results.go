// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"strings"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

// ResultsMap maps aliases to points
type ResultsMap map[string]Points

// Define defines labels
//  labels -- an alias to a group of points, an individual point, or labels to many points
//            example: "A", "left-column" or "a b c". If the number of points found is different
//            than the number of labels, a group is created.
//  Note:
//    To use spaces in labels, prefix the label with an exclamation mark; e.g "!right column"
func Define(labels string, loc Locator) {

	// check
	if len(labels) < 1 {
		chk.Panic("labels must have at least one character. %q is invalid", labels)
	}

	// locate points
	pts := loc.Locate()
	if len(pts) < 1 {
		chk.Panic("cannot define entities with alias/labels=%q and locator=%v", labels, loc)
	}

	// set results map
	if labels[0] == '!' {
		R[labels[1:]] = pts
		return
	}
	lbls := strings.Fields(labels)
	if len(lbls) == len(pts) {
		for i, l := range lbls {
			R[l] = []*Point{pts[i]}
		}
		return
	}
	R[labels] = pts
}

// LoadResults loads all results after points are defined
//  times -- specified selected output times
//           use nil to indicate that all times are required
func LoadResults(times []float64) {

	// selected output times and indices
	if times == nil {
		times = Sum.Times
	}
	I, T = utl.GetITout(Sum.Times, times, TolT)

	// for each selected output time
	for _, tidx := range I {

		// input results into domain
		if !Dom.In(tidx) {
			chk.Panic("cannot load results into domain; please check log file")
		}

		// for each quantity
		for _, pts := range R {
			for _, p := range pts {

				// node or integration point id
				nid := p.Nid
				pid := p.IpId

				// handle node
				if nid >= 0 {
					nod := Dom.Nodes[nid]
					for _, dof := range nod.Dofs {
						if dof != nil {
							utl.StrDblsMapAppend(&p.Vals, dof.Key, Dom.Sol.Y[dof.Eq])
						}
					}
				}

				// handle integration point
				if pid >= 0 {
					dat := Ipoints[pid]
					for key, val := range dat.V {
						utl.StrDblsMapAppend(&p.Vals, key, *val)
					}
				}
			}
		}
	}
}

// GetX gets a time series corresponding to a given label
//  tidx -- the time-output-index if label corresponds to a set of points;
//          thus the resulting slice is a spatial series @ time T[tidx].
//          If label corresponts to a single point, a time series is returned
//          and tidx is ignored.
func GetX(key, label string, tidx int) []float64 {
	if tidx < 0 {
		tidx = I[len(I)-1]
	}
	if pts, ok := R[label]; ok {
		if len(pts) == 1 {
			for k, v := range pts[0].Vals {
				if k == key {
					return v
				}
			}
		} else {
			res := make([]float64, len(pts))
			for i, p := range pts {
				for k, v := range p.Vals {
					if k == key {
						res[i] = v[tidx]
					}
				}
			}
			return res
		}
	}
	chk.Panic("cannot get %q at %q", key, label)
	return nil
}

// GetCoords returns the coordinates of a single point
func GetCoords(label string) []float64 {
	if pts, ok := R[label]; ok {
		if len(pts) == 1 {
			return pts[0].X
		}
	}
	chk.Panic("cannot get coordinates of point with label %q (make sure this label corresponds to a single point)", label)
	return nil
}

// GetDist returns the distance from a reference point on the given line with selected points
func GetDist(label string) []float64 {
	if pts, ok := R[label]; ok {
		res := make([]float64, len(pts))
		for i, p := range pts {
			res[i] = p.Dist
		}
		return res
	}
	chk.Panic("cannot get distance with label %q", label)
	return nil
}
