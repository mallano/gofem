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
func InitLogFile(dirout, fnamekey string) {

	// create log file
	var rank int
	if mpi.IsOn() {
		rank = mpi.Rank()
	}
	logFile, err := os.Create(utl.Sf("%s/%s_p%d.log", dirout, fnamekey, rank))
	if err != nil {
		utl.Panic("cannot create log file: %v", err)
	}

	// connect logger to output file
	log.SetOutput(logFile)
}

// FlusLog saves log (flushes to disk)
func FlushLog() {
	logFile.Close()
}
