// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"testing"

	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

func Test_out01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("out01")

	stageIdx := 0
	regionIdx := 0
	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"

	if ok := With(datadir+"p01.sim", stageIdx, regionIdx, func() {

		// get first ip coordinates
		xip := Ipoints[0].X
		utl.Pfcyan("xip = %v\n", xip)

		// configure time-plots
		Tplot("pl", &At{2.5, 0}, nil)
		Tplot("pl", &At{2.5, 10}, nil)
		Tplot("sl", &At{xip[0], xip[1]}, nil)

		// check slices
		nnod := 27
		nele := 4
		nip := 4
		utl.IntAssert(len(Dom.Nodes), nnod)
		utl.IntAssert(len(Ipoints), nele*nip)
		utl.IntAssert(len(TplotKeys), 2)
		utl.IntAssert(len(TplotData), 2)
		utl.CompareStrs(tst, "TplotKeys", TplotKeys, []string{"pl", "sl"})

		// check quantities
		for i, dat := range TplotData {
			key := TplotKeys[i]
			utl.Pforan("key=%v => dat=%v\n", key, dat)
			if key == "pl" {
				utl.IntAssert(len(dat.Qts), 2)
				utl.IntAssert(len(dat.Sty), 2)
			}
			if key == "sl" {
				utl.IntAssert(len(dat.Qts), 1)
				utl.IntAssert(len(dat.Sty), 1)
			}
		}

		// do plot
		_, err := PlotAll()
		if err != nil {
			return
		}

		// show figure
		if !utl.Tsilent {
			Show()
		}

	}); !ok {
		tst.Errorf("test failed\n")
	}
}

func Test_out02(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("out02")

	stageIdx := 0
	regionIdx := 0
	datadir := "$GOPATH/src/github.com/cpmech/gofem/fem/data/"

	if ok := With(datadir+"p01.sim", stageIdx, regionIdx, func() {

		// configure time-plots
		Tplot("pl", &At{2.5, 0}, Styles{{Label: "A", Marker: "o"}})
		Tplot("pl", &At{2.5, 10}, Styles{{Label: "B"}})

		// do plot
		sp, err := PlotAll()
		if err != nil {
			return
		}

		// settings
		I := sp["pl"]
		plt.Subplot(I[0], I[1], I[2])
		plt.AxisYrange(-10, 110)
		plt.Gll("$t$", "$p_{\\ell}$", "")

		// show figure
		if !utl.Tsilent {
			Show()
		}

	}); !ok {
		tst.Errorf("test failed\n")
	}
}
