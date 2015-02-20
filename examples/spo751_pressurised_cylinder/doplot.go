// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func main() {
	//utl.Tsilent = false
	defer func() {
		if err := recover(); err != nil {
			utl.Pf("[1;31mERROR: %v [0m\n", err)
		}
	}()

	ok := out.Start("spo751.sim", 0, 0)
	if !ok {
		utl.Panic("Start failed")
	}
	defer out.End()

	// commands for reading time-series
	out.Tseries("ux", &out.Verts{20}, out.Styles{{Label: "Outer", Marker: "o"}})

	// apply commands
	err := out.Apply()
	if err != nil {
		utl.Panic("test failed: %v\n", err)
	}

	ux := out.TseriesRes[out.TseriesKey2idx["ux"]][0]
	p := make([]float64, len(out.TseriesTimes))
	for i := 0; i < len(p); i++ {
		p[i] = out.TseriesTimes[i] * 0.2 / 0.96
	}

	plt.Plot(ux, p, "'ro-'")
	plt.Show()

	// show figure
	//out.Show(nil)
}
