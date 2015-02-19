// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"math/rand"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

func Test_imap(tst *testing.T) {

	//utl.Tsilent = false
	chk.PrintTitle("Test imap")

	for name, shape := range factory {
		gndim := shape.Gndim
		if gndim == 1 {
			continue
		}

		io.Pfyel("--------------------------------- %-6s---------------------------------\n", name)

		// check inverse mapping
		tol := 1e-14
		noise := 0.01
		if name == "tri10" {
			tol = 1e-14
		}
		if shape.FaceNverts > 2 {
			noise = 0.0
		}
		nverts := shape.Nverts
		C := la.MatAlloc(gndim, nverts)
		s := []float64{rand.Float64(), rand.Float64(), rand.Float64()} // scale factors
		la.MatCopy(C, 1.0, shape.NatCoords)
		_ = tol
		io.Pf("nverts:%v\n", nverts)
		io.Pf("gndim:%v\n", gndim)
		for i := 0; i < gndim; i++ {
			for j := 0; j < nverts; j++ {
				C[i][j] *= s[i]
				C[i][j] += noise * rand.Float64() // noise
			}
		}

		r := make([]float64, 3)
		x := make([]float64, 3)
		R := la.MatAlloc(gndim, nverts)

		for j := 0; j < nverts; j++ {
			for i := 0; i < gndim; i++ {
				x[i] = C[i][j]
			}
			err := shape.InvMap(r, x, C)
			io.Pf("r:%v\n", r)
			_ = err
			for i := 0; i < gndim; i++ {
				R[i][j] = r[i]
			}
		}

		chk.Matrix(tst, "checking", tol, R, shape.NatCoords)

		io.PfGreen("OK\n")
	}
}
