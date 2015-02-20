// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package msolid

import (
	"testing"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

func Test_path01(tst *testing.T) {

	//verbose()
	chk.PrintTitle("path01")

	ndim := 2
	var pth Path
	err := pth.ReadJson(ndim, "data/path01.json")
	if err != nil {
		tst.Errorf("test failed: %v\n", err)
		return
	}
	io.Pforan("pth = %+v\n", pth)
	chk.Vector(tst, "sx", 1e-17, pth.Sx, []float64{1, 1, 1, 0})
}
