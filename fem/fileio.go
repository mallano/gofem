// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"bytes"
	"encoding/gob"
	"encoding/json"
	"io"
	"os"
	"path"

	"github.com/cpmech/gosl/utl"
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
func GetEncoder(w io.Writer) Encoder {
	if global.Sim.Data.Encoder == "gob" {
		return gob.NewEncoder(w)
	}
	return json.NewEncoder(w)
}

// GetDecoder returns a new decoder
func GetDecoder(r io.Reader) Decoder {
	if global.Sim.Data.Encoder == "gob" {
		return gob.NewDecoder(r)
	}
	return json.NewDecoder(r)
}

// SaveSol saves Solution to a file which name is set with tidx (time output index)
func (o Domain) SaveSol(tidx int) (err error) {

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// encode t
	err = enc.Encode(o.Sol.T)
	if err != nil {
		return
	}

	// encode Y
	err = enc.Encode(o.Sol.Y)
	if err != nil {
		return
	}

	// transient simulations
	if !global.Sim.Data.Steady {

		// encode Dydt
		err = enc.Encode(o.Sol.Dydt)
		if err != nil {
			return
		}

		// encode D2ydt2
		err = enc.Encode(o.Sol.D2ydt2)
		if err != nil {
			return
		}
	}

	// save file
	fil, err := os.Create(out_nod_path(tidx))
	if err != nil {
		return
	}
	defer fil.Close()
	fil.Write(buf.Bytes())
	return
}

// ReadSol reads Solution from a file which name is set with tidx (time output index)
func (o *Domain) ReadSol(tidx int) (err error) {

	// open file
	fil, err := os.Open(out_nod_path(tidx))
	if err != nil {
		return
	}
	defer fil.Close()

	// get decoder
	dec := GetDecoder(fil)

	// decode t
	err = dec.Decode(&o.Sol.T)
	if err != nil {
		return
	}

	// decode Y
	err = dec.Decode(&o.Sol.Y)
	if err != nil {
		return
	}

	// transient simulations
	if !global.Sim.Data.Steady {

		// decode Dydt
		err = dec.Decode(&o.Sol.Dydt)
		if err != nil {
			return
		}

		// decode D2ydt2
		err = dec.Decode(&o.Sol.D2ydt2)
		if err != nil {
			return
		}
	}
	return
}

// SaveIvs saves elements's internal values to a file which name is set with tidx (time output index)
func (o Domain) SaveIvs(tidx int) (err error) {

	// buffer and encoder
	var buf bytes.Buffer
	enc := GetEncoder(&buf)

	// encode internal variables
	for _, e := range o.ElemWriters {
		err = e.Encode(enc)
		if err != nil {
			return
		}
	}

	// save file
	fil, err := os.Create(out_ele_path(tidx))
	if err != nil {
		return
	}
	defer fil.Close()
	fil.Write(buf.Bytes())
	return
}

// ReadSol reads elements's internal values from a file which name is set with tidx (time output index)
func (o *Domain) ReadIvs(tidx int) (err error) {

	// open file
	fil, err := os.Open(out_ele_path(tidx))
	if err != nil {
		return
	}
	defer fil.Close()

	// decode internal variables
	dec := GetDecoder(fil)
	for _, e := range o.ElemWriters {
		err = e.Decode(dec)
		if err != nil {
			return
		}
	}
	return
}

// Out performs output of Solution and Internal values to files
func (o *Domain) Out(tidx int) (err error) {
	err = o.SaveSol(tidx)
	if err != nil {
		return
	}
	err = o.SaveIvs(tidx)
	return
}

// In performes the inverse operation from Out
func (o *Domain) In(tidx int) (err error) {
	err = o.ReadIvs(tidx)
	if err != nil {
		return
	}
	err = o.ReadSol(tidx)
	return
}

// auxiliary ///////////////////////////////////////////////////////////////////////////////////////

func out_nod_path(tidx int) string {
	return path.Join(global.Sim.Data.DirOut, utl.Sf("%s_nod_%010d.%s", global.Sim.Data.FnameKey, tidx, global.Sim.Data.Encoder))
}

func out_ele_path(tidx int) string {
	return path.Join(global.Sim.Data.DirOut, utl.Sf("%s_ele_%010d_p%d.%s", global.Sim.Data.FnameKey, tidx, global.Rank, global.Sim.Data.Encoder))
}
