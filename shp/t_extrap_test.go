// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

import (
	"math/rand"
	"strconv"
	"strings"
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/la"
)

func Test_extrap(tst *testing.T) {

	//utl.Tsilent = false
	chk.PrintTitle("Test extrapolation")

	GetIpNums := func(name string) []int {
		var nums []int
		for key, _ := range ipsfactory {
			name_n := strings.Split(key, "_")
			if name_n[0] == name {
				nip := io.Atoi(name_n[1])
				nums = append(nums, nip)
			}
		}
		return nums
	}

	for name, shape := range factory {

		gndim := shape.Gndim
		if gndim == 1 {
			continue
		}

		io.Pfyel("--------------------------------- %-6s---------------------------------\n", name)

		for _, nip := range GetIpNums(shape.Type) {
			if nip <= 1 {
				continue
			}
			io.Pfblue("nip = %v\n", nip)
			tol := 1.0e-13

			// create a N vector with nodal values
			nverts := shape.Nverts
			X := shape.NatCoords
			delta := rand.Float64()
			N := make([]float64, shape.Nverts)
			for i := 0; i < shape.Nverts; i++ {
				var x, y, z float64
				x = X[0][i]
				y = X[1][i]
				if shape.Gndim == 3 {
					z = X[2][i]
				}
				N[i] = x + y + z + delta
			}
			io.Pfblue("N    = %v\n", N)

			// calculate P vector with corresponding values at ips
			key := name + "_" + strconv.Itoa(nip)
			ips := ipsfactory[key]
			P := make([]float64, nip)
			Xip := la.MatAlloc(nip, 4) // ips local coordinates
			for i := 0; i < nip; i++ {
				var x, y, z float64
				x = ips[i].R
				y = ips[i].S
				z = ips[i].T
				Xip[i][0] = x
				Xip[i][1] = y
				Xip[i][2] = z
				P[i] = x + y + z + delta
			}

			// Allocate E matrix
			E := la.MatAlloc(nverts, nip)

			// Calculate extrapolator matrix
			shape.Extrapolator(E, Xip)

			// Recalculate nodal values NN = E*P
			NN := make([]float64, shape.Nverts)
			la.MatVecMul(NN, 1.0, E, P)
			io.Pfblue("N ext= %v\n", NN)

			// Compare vectors N and NN
			msg := name + "_" + strconv.Itoa(nip)
			chk.Vector(tst, msg, tol, NN, N)
		}

	}
}
