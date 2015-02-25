// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"strings"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/utl"
)

// ResultsMap maps aliases to points
type ResultsMap map[string]Points

// Define defines aliases
//  alias -- an alias to a group of points, an individual point, or to a set of points.
//           Example: "A", "left-column" or "a b c". If the number of points found is different
//           than the number of aliases, a group is created.
//  Note:
//    To use spaces in aliases, prefix the alias with an exclamation mark; e.g "!right column"
func Define(alias string, loc Locator) {

	// check
	if len(alias) < 1 {
		chk.Panic("alias must have at least one character. %q is invalid", alias)
	}

	// locate points
	pts := loc.Locate()
	if len(pts) < 1 {
		chk.Panic("cannot define entities with alias=%q and locator=%v", alias, loc)
	}

	// set results map
	if alias[0] == '!' {
		R[alias[1:]] = pts
		return
	}
	lbls := strings.Fields(alias)
	if len(lbls) == len(pts) {
		for i, l := range lbls {
			R[l] = []*Point{pts[i]}
		}
		return
	}
	R[alias] = pts
}

// LoadResults loads all results after points are defined
//  times -- specified selected output times
//           use nil to indicate that all times are required
func LoadResults(times []float64) {

	// selected output times and indices
	if times == nil {
		times = fem.Global.Sum.OutTimes
	}
	I, T = utl.GetITout(fem.Global.Sum.OutTimes, times, TolT)

	// for each selected output time
	for _, tidx := range I {

		// input results into domain
		if !Dom.In(tidx) {
			chk.Panic("cannot load results into domain; please check log file")
		}

		// for each point
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

// GetRes gets results as a time or space series corresponding to a given alias
// for a single point or set of points.
//  idxI -- index in I slice corresponding to selected output time; use -1 for the last item.
//          If alias defines a single point, the whole time series is returned and idxI is ignored.
func GetRes(key, alias string, idxI int) []float64 {
	if idxI < 0 {
		idxI = len(I) - 1
	}
	if pts, ok := R[alias]; ok {
		if len(pts) == 1 {
			for k, v := range pts[0].Vals {
				if k == key {
					return v
				}
			}
		} else {
			var res []float64
			for _, p := range pts {
				for k, v := range p.Vals {
					if k == key {
						res = append(res, v[idxI])
					}
				}
			}
			return res
		}
	}
	chk.Panic("cannot get %q at %q", key, alias)
	return nil
}

// GetCoords returns the coordinates of a single point
func GetCoords(alias string) []float64 {
	if pts, ok := R[alias]; ok {
		if len(pts) == 1 {
			return pts[0].X
		}
	}
	chk.Panic("cannot get coordinates of point with alias %q (make sure this alias corresponds to a single point)", alias)
	return nil
}

// GetDist returns the distance from a reference point on the given line with selected points
// if they contain a given key
//  key -- use any to get coordinates of points with any key such as "ux", "pl", etc.
func GetDist(key, alias string) (dist []float64) {
	any := key == "any"
	if pts, ok := R[alias]; ok {
		for _, p := range pts {
			for k, _ := range p.Vals {
				if k == key || any {
					dist = append(dist, p.Dist)
					break
				}
			}
		}
		return
	}
	chk.Panic("cannot get distance with key %q and alias %q", key, alias)
	return
}

// GetXYZ returns the x-y-z coordinates of selected points that have a specified key
//  key -- use any to get coordinates of points with any key such as "ux", "pl", etc.
func GetXYZ(key, alias string) (x, y, z []float64) {
	any := key == "any"
	if pts, ok := R[alias]; ok {
		for _, p := range pts {
			for k, _ := range p.Vals {
				if k == key || any {
					x = append(x, p.X[0])
					y = append(y, p.X[1])
					if len(p.X) == 3 {
						z = append(z, p.X[2])
					}
					break
				}
			}
		}
		return
	}
	chk.Panic("cannot get x-y-z coordinates with key %q and alias %q", key, alias)
	return
}
