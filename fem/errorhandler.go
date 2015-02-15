// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"log"

	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

func LogErr(err error, msg string) (stop bool) {
	if err != nil {
		fullmsg := "ERROR: " + msg + " : " + err.Error()
		log.Printf(fullmsg)
		Global.WspcStop[Global.Rank] = 1
		return true
	}
	return
}

func LogErrCond(condition bool, msg string, prm ...interface{}) (stop bool) {
	if condition {
		fullmsg := "ERROR: " + utl.Sf(msg, prm...)
		log.Printf(fullmsg)
		Global.WspcStop[Global.Rank] = 1
		return true
	}
	return
}

func Stop() bool {
	if !Global.Distr {
		if Global.WspcStop[Global.Rank] > 0 {
			utl.CallerInfo(3)
			utl.CallerInfo(2)
			utl.PfRed("simulation stopped due to errors. see log files\n")
			return true
		}
		return false
	}
	mpi.IntAllReduceMax(Global.WspcStop, Global.WspcInum)
	for i := 0; i < Global.Nproc; i++ {
		if Global.WspcStop[i] > 0 {
			if Global.Root {
				utl.CallerInfo(3)
				utl.CallerInfo(2)
				utl.PfRed("simulation stopped due to errors. see log files\n")
			}
			return true
		}
	}
	return false
}
