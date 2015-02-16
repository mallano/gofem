// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"sort"

	"github.com/cpmech/gosl/utl"
)

// TplotDat holds information of one key to be ploted along time; e.g. "pl"
type TplotDat struct {
	Qts Quantities // [npoints] points where values must be plotted
	Sty Styles     // [npoints] styles => one item can have many points
}

// Global variables
var (
	TplotKeys []string    // [nkeys] all keys
	TplotData []*TplotDat // [nkeys] all items
)

// TplotClear clears Tplot data
func TplotClear() {
	TplotKeys = make([]string, 0)
	TplotData = make([]*TplotDat, 0)
}

// TplotStart initialise structures for Time Plots
//  Returns results: R[nkeys][nqts...?][ntimes]
func TplotStart() {
	nkeys := len(TplotKeys)
	R = make([][][]float64, nkeys)
	for i, dat := range TplotData {
		R[i] = make([][]float64, len(dat.Qts))
		for j, _ := range dat.Qts {
			R[i][j] = make([]float64, Sum.NumTidx)
		}
		sort.Sort(dat.Qts)
	}
	return
}

// Tplot specifies a variable at a particular point to be plotted along time
func Tplot(key string, loc PointLocator, sty Styles) {
	qts := loc.AtPoint(key)
	if qts == nil {
		return
	}
	n := len(qts)
	if len(sty) != n {
		sty = GetDefaultStyles(n)
	}
	idx := utl.StrIndexSmall(TplotKeys, key)
	if idx < 0 {
		TplotKeys = append(TplotKeys, key)
		TplotData = append(TplotData, &TplotDat{Qts: qts, Sty: sty})
		return
	}
	TplotData[idx].Qts = append(TplotData[idx].Qts, qts...)
	TplotData[idx].Sty = append(TplotData[idx].Sty, sty...)
}
