// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"log"
	"os"

	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

// errFile holds a handle to errors logger file
var logFile os.File

// InitLogFile initialises logger
func InitLogFile(dirout, fnamekey string) (err error) {

	// create log file
	var rank int
	if mpi.IsOn() {
		rank = mpi.Rank()
	}
	logFile, err := os.Create(utl.Sf("%s/%s_p%d.log", dirout, fnamekey, rank))
	if err != nil {
		return
	}

	// connect logger to output file
	log.SetOutput(logFile)
	return
}

// FlusLog saves log (flushes to disk)
func FlushLog() {
	logFile.Close()
}

// LogErr logs error and returs stop flag
func LogErr(err error, msg string) (stop bool) {
	if err != nil {
		fullmsg := "ERROR: " + msg + " : " + err.Error()
		log.Printf(fullmsg)
		return true
	}
	return false
}

// LogErr logs error using condition (==true) to stop and returs stop flag
func LogErrCond(condition bool, msg string, prm ...interface{}) (stop bool) {
	if condition {
		fullmsg := "ERROR: " + utl.Sf(msg, prm...)
		log.Printf(fullmsg)
		return true
	}
	return false
}
