// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"sort"

	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

type Styles []*plt.LineData

func GetDefaultStyles(n int) Styles {
	return make([]*plt.LineData, n)
}

// TplotDat holds information of one key to be ploted along time; e.g. "pl"
type TplotDat struct {
	Qts Quantities // [npoints] points where values must be plotted
	Sty Styles     // [npoints] styles => one item can have many points
}

var (
	TplotKeys []string    // [nkeys] all keys
	TplotData []*TplotDat // [nkeys] all items
)

func TplotClear() {
	TplotKeys = make([]string, 0)
	TplotData = make([]*TplotDat, 0)
}

// R[nkeys][nqts...?][ntimes]
func TplotStart() (R [][][]float64) {
	nkeys := len(TplotKeys)
	R = make([][][]float64, nkeys)
	for i, dat := range TplotData {
		R[i] = make([][]float64, len(dat.Qts))
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

/*
func get_tplot_quantities() (T []float64, V map[string][]float64, err error) {
	utl.Pforan("Sum = %v\n", Sum)
	T = make([]float64, Sum.NumTidx)
	V = make(map[string][]float64)
	for tidx := 0; tidx < Sum.NumTidx; tidx++ {
		if !Dom.ReadSol(tidx) {
			return nil, nil, utl.Err("ReadSol failed. See log files\n")
		}
		utl.Pforan("tidx = %v\n", tidx)
		T[tidx] = Dom.Sol.T
		for _, key := range TplotKeys {
			for _, item := range TplotData[key] {
				Q := item.Loc.AtPoint(key)
				for _, q := range Q {
					utl.StrDblsMapAppend(&V, key, q.Value)
				}
			}
		}
	}
	return
}
*/
