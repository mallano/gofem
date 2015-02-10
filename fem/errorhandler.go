// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

// Stop decides whether a serial or parallel execution has to be stopped or not
func Stop(err error, msg string) bool {

	// serial run
	if !global.Distr {
		if err != nil {
			utl.Pf("\n")
			utl.CallerInfo(3)
			utl.CallerInfo(2)
			utl.PfMag("simulation failed on %s with %v\n", msg, err)
			return true // stop
		}
		return false // continue
	}

	// parallel run
	for i := 0; i < global.Nproc; i++ {
		global.WspcStop[i] = 0 // all continue
	}
	if err != nil {
		utl.Pf("\n")
		utl.CallerInfo(3)
		utl.CallerInfo(2)
		utl.PfMag("simulation failed in proc # %d on %s with %v\n", global.Rank, msg, err)
		global.WspcStop[global.Rank] = 1 // this processor wants to stop
	}
	mpi.IntAllReduceMax(global.WspcStop, global.WspcInum)
	for i := 0; i < global.Nproc; i++ {
		if global.WspcStop[i] > 0 {
			return true // stop
		}
	}
	return false // continue
}

// PanicOrNot decides to panic if any processor wants to panic
func PanicOrNot(dopanic bool, msg string, prm ...interface{}) {

	// serial run
	if !global.Distr {
		if dopanic {
			utl.Pf("\n")
			panic(utl.Sf(msg, prm...))
		}
		return
	}

	// parallel run
	for i := 0; i < global.Nproc; i++ {
		global.WspcStop[i] = 0 // all ok
	}
	if dopanic {
		global.WspcStop[global.Rank] = 1 // this processor wants to panic
	}
	mpi.IntAllReduceMax(global.WspcStop, global.WspcInum)
	for i := 0; i < global.Nproc; i++ {
		if global.WspcStop[i] > 0 {
			if global.Root {
				utl.Pf("\n")
				utl.CallerInfo(3)
				utl.CallerInfo(2)
				panic(utl.Sf(msg, prm...))
			} else {
				panic(utl.Sf(msg, prm...))
			}
		}
	}
}
