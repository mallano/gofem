// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"os"
	"path"

	"github.com/cpmech/gosl/io"
)

// Summary records summary of outputs
type Summary struct {
	Times   []float64 // [nOutTimes] output times
	StgTidx []int     // [nstg] first stage's time-output-index
}

// SaveSums saves summary to disc
func SaveSum(sum *Summary) (ok bool) {

	// skip if not root
	if !Global.Root {
		return true
	}

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// encode summary
	if LogErr(enc.Encode(sum), "SaveSum") {
		return
	}

	// save file
	return save_file("SaveSum", "summary", out_sum_path(Global.Rank), &buf)
}

// ReadSum reads summary back
func ReadSum() *Summary {

	// open file
	fil, err := os.Open(out_sum_path(0)) // read always from proc # 0
	if LogErr(err, "ReadSum") {
		return nil
	}
	defer func() {
		if LogErr(fil.Close(), "ReadSum: cannot close file") {
			return
		}
	}()

	// decode summary
	var sum Summary
	dec := GetDecoder(fil)
	err = dec.Decode(&sum)
	if LogErr(err, "ReadSum") {
		return nil
	}
	return &sum
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_sum_path(proc int) string {
	return path.Join(Global.Sim.Data.DirOut, io.Sf("%s_p%d_sum.%s", Global.Sim.Data.FnameKey, proc, Global.Sim.Data.Encoder))
}
