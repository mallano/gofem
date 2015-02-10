// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/utl"
)

func Test_state01(tst *testing.T) {

	prevTs := utl.Tsilent
	defer func() {
		utl.Tsilent = prevTs
		if err := recover(); err != nil {
			tst.Error("[1;31mSome error has happened:[0m\n", err)
		}
	}()

	//utl.Tsilent = false
	utl.TTitle("state01")

	nsig, nalp, nphi, large := 4, 1, 3, false
	state0 := NewState(nsig, nalp, nphi, large)
	utl.Pforan("state0 = %+v\n", state0)
	utl.CompareDbls(tst, "sig", state0.Sig, []float64{0, 0, 0, 0})
	utl.CompareDbls(tst, "alp", state0.Alp, []float64{0})
	utl.CompareDbls(tst, "phi", state0.Phi, []float64{0, 0, 0})

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
	utl.Pforan("state1 = %+v\n", state1)
	utl.CompareDbls(tst, "sig", state1.Sig, []float64{10, 11, 12, 13})
	utl.CompareDbls(tst, "alp", state1.Alp, []float64{20})
	utl.CompareDbls(tst, "phi", state1.Phi, []float64{30, 31, 32})

	state2 := state1.GetCopy()
	utl.Pforan("state2 = %+v\n", state2)
	utl.CompareDbls(tst, "sig", state2.Sig, []float64{10, 11, 12, 13})
	utl.CompareDbls(tst, "alp", state2.Alp, []float64{20})
	utl.CompareDbls(tst, "phi", state2.Phi, []float64{30, 31, 32})
}
