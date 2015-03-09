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
	Nproc    int         // number of processors used in last last run; equal to 1 if not distributed
	OutTimes []float64   // [nOutTimes] output times
	Resids   utl.DblList // [nTimes][nIter] residuals (if Stat is on; includes all stages)
}

// SaveSums saves summary to disc
func (o Summary) Save() (ok bool) {

	// set flags before saving
	o.Nproc = Global.Nproc

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
	fn := out_sum_path(Global.Dirout, Global.Fnkey, Global.Rank)
	return save_file("SaveSum", "summary", fn, &buf)
}

// ReadSum reads summary back
//  Note: returns nil on errors
func ReadSum(dirout, fnamekey string) (o *Summary) {

	// open file
	fn := out_sum_path(dirout, fnamekey, 0) // reading always from proc # 0
	fil, err := os.Open(fn)
	if LogErr(err, "ReadSum") {
		return nil
	}
	defer func() {
		LogErr(fil.Close(), "ReadSum: cannot close file")
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

func out_sum_path(dirout, fnamekey string, proc int) string {
	return path.Join(dirout, io.Sf("%s_p%d_sum.%s", fnamekey, proc, Global.Enc))
}
