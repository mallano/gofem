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

	// encode t
	if LogErr(enc.Encode(o.Sol.T), "SaveSol") {
		return
	}

	// encode Y
	if LogErr(enc.Encode(o.Sol.Y), "SaveSol") {
		return
	}

	// encode L
	if LogErr(enc.Encode(o.Sol.L), "SaveSol") {
		return
	}

	// transient simulations
	if !Global.Sim.Data.Steady {

		// encode Dydt
		if LogErr(enc.Encode(o.Sol.Dydt), "SaveSol") {
			return
		}

		// encode D2ydt2
		if LogErr(enc.Encode(o.Sol.D2ydt2), "SaveSol") {
			return
		}
	}

	// save file
	return save_file("SaveSol", "solution", out_nod_path(tidx, Global.Rank), &buf)
}

// ReadSol reads Solution from a file which name is set with tidx (time output index)
func (o *Domain) ReadSol(tidx int) (ok bool) {

	// open file
	fil, err := os.Open(out_nod_path(tidx, 0)) // read always from proc # 0
	if LogErr(err, "ReadSol") {
		return
	}
	defer func() {
		LogErr(fil.Close(), "ReadSol: cannot close file")
	}()

	// get decoder
	dec := GetDecoder(fil)

	// decode t
	if LogErr(dec.Decode(&o.Sol.T), "ReadSol") {
		return
	}

	// decode Y
	if LogErr(dec.Decode(&o.Sol.Y), "ReadSol") {
		return
	}

	// decode L
	if LogErr(dec.Decode(&o.Sol.L), "ReadSol") {
		return
	}

	// transient simulations
	if !Global.Sim.Data.Steady {

		// decode Dydt
		if LogErr(dec.Decode(&o.Sol.Dydt), "ReadSol") {
			return
		}

		// decode D2ydt2
		if LogErr(dec.Decode(&o.Sol.D2ydt2), "ReadSol") {
			return
		}
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
	return save_file("SaveIvs", "internal values", out_ele_path(tidx, Global.Rank), &buf)
}

// ReadIvs reads elements's internal values from a file which name is set with tidx (time output index)
func (o *Domain) ReadIvs(tidx, proc int) (ok bool) {

	// open file
	fil, err := os.Open(out_ele_path(tidx, proc))
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
		if !o.ReadIvs(tidx, i) {
			return
		}
	}
	return o.ReadSol(tidx)
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_nod_path(tidx, proc int) string {
	return path.Join(Global.Dirout, io.Sf("%s_p%d_nod_%010d.%s", Global.Fnkey, proc, tidx, Global.Enc))
}

func out_ele_path(tidx, proc int) string {
	return path.Join(Global.Dirout, io.Sf("%s_p%d_ele_%010d.%s", Global.Fnkey, proc, tidx, Global.Enc))
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
