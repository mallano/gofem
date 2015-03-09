// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"encoding/gob"
	"encoding/json"
	goio "io"
	"os"
	"path"

	"github.com/cpmech/gosl/io"
)

// Encoder defines encoders; e.g. gob or json
type Encoder interface {
	Encode(e interface{}) error
}

// Decoder defines decoders; e.g. gob or json
type Decoder interface {
	Decode(e interface{}) error
}

// GetEncoder returns a new encoder
func GetEncoder(w goio.Writer) Encoder {
	if Global.Enc == "json" {
		return json.NewEncoder(w)
	}
	return gob.NewEncoder(w)
}

// GetDecoder returns a new decoder
func GetDecoder(r goio.Reader) Decoder {
	if Global.Enc == "json" {
		return json.NewDecoder(r)
	}
	return gob.NewDecoder(r)
}

// SaveSol saves Solution to a file which name is set with tidx (time output index)
func (o Domain) SaveSol(tidx int) (ok bool) {

	// skip if not root
	if !Global.Root {
		return true
	}

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// encode Sol
	if LogErr(enc.Encode(o.Sol), "SaveSol") {
		return
	}

	// save file
	fn := out_nod_path(Global.Dirout, Global.Fnkey, tidx, Global.Rank)
	return save_file("SaveSol", "solution", fn, &buf)
}

// ReadSol reads Solution from a file which name is set with tidx (time output index)
func (o *Domain) ReadSol(dir, fnkey string, tidx int) (ok bool) {

	// open file
	fn := out_nod_path(dir, fnkey, tidx, 0) // reading always from proc # 0
	fil, err := os.Open(fn)
	if LogErr(err, "ReadSol") {
		return
	}
	defer func() {
		LogErr(fil.Close(), "ReadSol: cannot close file")
	}()

	// get decoder
	dec := GetDecoder(fil)

	// decode Sol
	if LogErr(dec.Decode(&o.Sol), "ReadSol") {
		return
	}
	return true
}

// SaveIvs saves elements's internal values to a file which name is set with tidx (time output index)
func (o Domain) SaveIvs(tidx int) (ok bool) {

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// elements that go to file
	enc.Encode(o.MyCids)

	// encode internal variables
	for _, e := range o.Elems {
		if !e.Encode(enc) {
			return
		}
	}

	// save file
	fn := out_ele_path(Global.Dirout, Global.Fnkey, tidx, Global.Rank)
	return save_file("SaveIvs", "internal values", fn, &buf)
}

// ReadIvs reads elements's internal values from a file which name is set with tidx (time output index)
func (o *Domain) ReadIvs(dir, fnkey string, tidx, proc int) (ok bool) {

	// open file
	fn := out_ele_path(dir, fnkey, tidx, proc)
	fil, err := os.Open(fn)
	if LogErr(err, "ReadIvs") {
		return
	}
	defer func() {
		LogErr(fil.Close(), "ReadIvs: cannot close file")
	}()

	// decoder
	dec := GetDecoder(fil)

	// elements that are in file
	dec.Decode(&o.MyCids)

	// decode internal variables
	for _, cid := range o.MyCids {
		elem := o.Cid2elem[cid]
		if LogErrCond(elem == nil, "ReadIvs: cannot find element with cid=%d", cid) {
			return
		}
		if !elem.Decode(dec) {
			return
		}
	}
	return true
}

// Out performs output of Solution and Internal values to files
func (o *Domain) Out(tidx int) (ok bool) {
	if !o.SaveSol(tidx) {
		return
	}
	return o.SaveIvs(tidx)
}

// In performes the inverse operation from Out
func (o *Domain) In(sum *Summary, tidx int) (ok bool) {
	for i := 0; i < sum.Nproc; i++ {
		if !o.ReadIvs(sum.Dirout, sum.Fnkey, tidx, i) {
			return
		}
	}
	return o.ReadSol(sum.Dirout, sum.Fnkey, tidx)
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_nod_path(dir, fnkey string, tidx, proc int) string {
	return path.Join(dir, io.Sf("%s_p%d_nod_%010d.%s", fnkey, proc, tidx, Global.Enc))
}

func out_ele_path(dir, fnkey string, tidx, proc int) string {
	return path.Join(dir, io.Sf("%s_p%d_ele_%010d.%s", fnkey, proc, tidx, Global.Enc))
}

func save_file(function, category, filename string, buf *bytes.Buffer) (ok bool) {
	fil, err := os.Create(filename)
	if LogErr(err, function) {
		return
	}
	defer func() {
		LogErr(fil.Close(), io.Sf("cannot close %s file", category))
	}()
	_, err = fil.Write(buf.Bytes())
	if LogErr(err, io.Sf("cannot write to %s file", category)) {
		return
	}
	return true
}
