// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

func Test_state01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mERROR:", err, "[0m\n")
		}
	}()

	//utl.Tsilent = false
	chk.PrintTitle("state01")

	nsig, nalp, nphi, large := 4, 1, 3, false
	state0 := NewState(nsig, nalp, nphi, large)
	io.Pforan("state0 = %+v\n", state0)
	chk.Vector(tst, "sig", state0.Sig, []float64{0, 0, 0, 0})
	chk.Vector(tst, "alp", state0.Alp, []float64{0})
	chk.Vector(tst, "phi", state0.Phi, []float64{0, 0, 0})

	state0.Sig[0] = 10.0
	state0.Sig[1] = 11.0
	state0.Sig[2] = 12.0
	state0.Sig[3] = 13.0
	state0.Alp[0] = 20.0
	state0.Phi[0] = 30.0
	state0.Phi[1] = 31.0
	state0.Phi[2] = 32.0

	state1 := NewState(nsig, nalp, nphi, large)
	state1.Set(state0)
	io.Pforan("state1 = %+v\n", state1)
	chk.Vector(tst, "sig", state1.Sig, []float64{10, 11, 12, 13})
	chk.Vector(tst, "alp", state1.Alp, []float64{20})
	chk.Vector(tst, "phi", state1.Phi, []float64{30, 31, 32})

	state2 := state1.GetCopy()
	io.Pforan("state2 = %+v\n", state2)
	chk.Vector(tst, "sig", state2.Sig, []float64{10, 11, 12, 13})
	chk.Vector(tst, "alp", state2.Alp, []float64{20})
	chk.Vector(tst, "phi", state2.Phi, []float64{30, 31, 32})
}
