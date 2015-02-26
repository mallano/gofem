// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"os"
	"path"

	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

// Summary records summary of outputs
type Summary struct {

	// essential
	OutTimes []float64   // [nOutTimes] output times
	Resids   utl.DblList // [nTimes][nIter] residuals (if Stat is on; includes all stages)

	// auxiliary
	cinc int // current time-step increment (used)
}

// SaveSums saves summary to disc
func (o Summary) Save() (ok bool) {

	// skip if not root
	if !Global.Root {
		return true
	}

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// encode summary
	if LogErr(enc.Encode(o), "SaveSum") {
		return
	}

	// save file
	return save_file("SaveSum", "summary", out_sum_path(Global.Rank), &buf)
}

// ReadSum reads summary back
func (o *Summary) Read() (ok bool) {

	// open file
	fil, err := os.Open(out_sum_path(0)) // read always from proc # 0
	if LogErr(err, "ReadSum") {
		return
	}
	defer func() {
		LogErr(fil.Close(), "ReadSum: cannot close file")
	}()

	// decode summary
	dec := GetDecoder(fil)
	err = dec.Decode(&o)
	if LogErr(err, "ReadSum") {
		return
	}
	return true
}

// AddResid adds the residual value for a given iteration
func (o Summary) AddResid(iter int, resid float64) {
	// update time-step counter
	if iter == 0 && o.cinc != 0 {
		o.cinc += 1
	}

	// add residual
	o.Resids.Append(o.cinc, resid)
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_sum_path(proc int) string {
	return path.Join(Global.Sim.Data.DirOut, io.Sf("%s_p%d_sum.%s", Global.Sim.Data.FnameKey, proc, Global.Sim.Data.Encoder))
}
