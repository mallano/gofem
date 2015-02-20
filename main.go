// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/mpi"
	"github.com/cpmech/gosl/utl"
)

func main() {

	// catch errors
	defer func() {
		if mpi.Rank() == 0 {
			if err := recover(); err != nil {
				io.PfRed("ERROR: %v\n", err)
			}
		}
		mpi.Stop(false)
	}()
	mpi.Start(false)

	// message
	io.PfWhite("\nGofem v3 -- Go Finite Element Method\n\n")
	io.Pf("Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.\n")
	io.Pf("Use of this source code is governed by a BSD-style\n")
	io.Pf("license that can be found in the LICENSE file.\n\n")

	// simulation filenamepath
	flag.Parse()
	var fnamepath string
	if len(flag.Args()) > 0 {
		fnamepath = flag.Arg(0)
	} else {
		chk.Panic("Please, provide a filename. Ex.: cylinder.sim")
	}

	// other options
	erasefiles := true
	verbose := true
	if len(flag.Args()) > 1 {
		erasefiles = io.Atob(flag.Arg(1))
	}
	if len(flag.Args()) > 2 {
		verbose = io.Atob(flag.Arg(2))
	}

	// profiling?
	defer utl.DoProf(false)()

	// start global variables and log
	if !fem.Start(fnamepath, erasefiles, verbose) {
		chk.Panic("Start failed\n")
		return
	}

	// make sure to flush log
	defer fem.End()

	// run simulation
	if !fem.Run() {
		chk.Panic("Run failed\n")
		return
	}
}
