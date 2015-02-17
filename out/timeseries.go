// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"sort"

	"github.com/cpmech/gosl/utl"
)

// TseriesDat holds information of one key to be collected along time; e.g. "pl"
type TseriesDat struct {
	Qts Quantities // [npoints] points where values must be analysed/plotted
	Sty Styles     // [npoints] styles => one item can have many points
}

// TseriesClear clears Tplot data
func TseriesClear() {
	TseriesKeys = make([]string, 0)
	TseriesData = make([]*TseriesDat, 0)
}

// TseriesStart initialise structures containing all time-series results
func TseriesStart() {
	nkeys := len(TseriesKeys)
	TseriesR = make([][][]float64, nkeys)
	for i, dat := range TseriesData {
		TseriesR[i] = make([][]float64, len(dat.Qts))
		for j, _ := range dat.Qts {
			TseriesR[i][j] = make([]float64, Sum.NumTidx)
		}
		sort.Sort(dat.Qts)
	}
	return
}

// Tseries specifies a variable at a particular point to be plotted along time
func Tseries(key string, loc PointLocator, sty Styles) {
	qts := loc.AtPoint(key)
	if qts == nil {
		return
	}
	n := len(qts)
	if len(sty) != n {
		sty = GetDefaultStyles(n)
	}
	idx := utl.StrIndexSmall(TseriesKeys, key)
	if idx < 0 {
		TseriesKeys = append(TseriesKeys, key)
		TseriesData = append(TseriesData, &TseriesDat{Qts: qts, Sty: sty})
		return
	}
	TseriesData[idx].Qts = append(TseriesData[idx].Qts, qts...)
	TseriesData[idx].Sty = append(TseriesData[idx].Sty, sty...)
}
