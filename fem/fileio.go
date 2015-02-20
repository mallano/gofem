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
	if Global.Sim.Data.Encoder == "json" {
		return json.NewEncoder(w)
	}
	return gob.NewEncoder(w)
}

// GetDecoder returns a new decoder
func GetDecoder(r goio.Reader) Decoder {
	if Global.Sim.Data.Encoder == "json" {
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
	fil, err := os.Create(out_nod_path(tidx, Global.Rank))
	if LogErr(err, "SaveSol") {
		return
	}
	defer fil.Close()
	fil.Write(buf.Bytes())
	return true
}

// ReadSol reads Solution from a file which name is set with tidx (time output index)
func (o *Domain) ReadSol(tidx int) (ok bool) {

	// open file
	fil, err := os.Open(out_nod_path(tidx, 0)) // read always from proc # 0
	if LogErr(err, "ReadSol") {
		return
	}
	defer fil.Close()

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

	// encode internal variables
	for _, e := range o.Elems {
		if !e.Encode(enc) {
			return
		}
	}

	// save file
	fil, err := os.Create(out_ele_path(tidx))
	if LogErr(err, "SaveIvs") {
		return
	}
	defer fil.Close()
	fil.Write(buf.Bytes())
	return true
}

// ReadIvs reads elements's internal values from a file which name is set with tidx (time output index)
func (o *Domain) ReadIvs(tidx int) (ok bool) {

	// open file
	fil, err := os.Open(out_ele_path(tidx))
	if LogErr(err, "ReadIvs") {
		return
	}
	defer fil.Close()

	// decode internal variables
	dec := GetDecoder(fil)
	for _, e := range o.Elems {
		if !e.Decode(dec) {
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
func (o *Domain) In(tidx int) (ok bool) {
	if !o.ReadIvs(tidx) {
		return
	}
	return o.ReadSol(tidx)
}

// Summary records summary of outputs
type Summary struct {
	TidxIni []int // [nstg] first stage's time-output-index
	NumTidx int   // number of tidx
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
	fil, err := os.Create(out_sum_path(Global.Rank))
	if LogErr(err, "SaveSum") {
		return
	}
	defer fil.Close()
	fil.Write(buf.Bytes())
	return true
}

// ReadSum reads summary back
func ReadSum() *Summary {

	// open file
	fil, err := os.Open(out_sum_path(0)) // read always from proc # 0
	if LogErr(err, "ReadSum") {
		return nil
	}
	defer fil.Close()

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

func out_nod_path(tidx, proc int) string {
	return path.Join(Global.Sim.Data.DirOut, io.Sf("%s_p%d_nod_%010d.%s", Global.Sim.Data.FnameKey, proc, tidx, Global.Sim.Data.Encoder))
}

func out_ele_path(tidx int) string {
	return path.Join(Global.Sim.Data.DirOut, io.Sf("%s_p%d_ele_%010d.%s", Global.Sim.Data.FnameKey, Global.Rank, tidx, Global.Sim.Data.Encoder))
}

func out_sum_path(proc int) string {
	return path.Join(Global.Sim.Data.DirOut, io.Sf("%s_p%d_sum.%s", Global.Sim.Data.FnameKey, proc, Global.Sim.Data.Encoder))
}
