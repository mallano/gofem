// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/io"
)

func GetIntegrationPoints(extra, cellType string) (ipsElem, ipsFace []*shp.Ipoint) {

	// get number of ips of element
	var nip int
	if s_nip, found := io.Keycode(extra, "nip"); found {
		nip = io.Atoi(s_nip)
	}

	// get number of ips of face
	var nipf int
	if s_nipf, found := io.Keycode(extra, "nipf"); found {
		nipf = io.Atoi(s_nipf)
	}

	// get integration points of element
	var err error
	ipsElem, err = shp.GetIps(cellType, nip)
	if LogErr(err, io.Sf("cannot get integration points for element with shape type=%q and nip=%d", cellType, nip)) {
		return nil, nil
	}

	// get integration points of face
	faceType := shp.GetFaceType(cellType)
	ipsFace, err = shp.GetIps(faceType, nipf)
	if LogErr(err, io.Sf("cannot get integration points for face with face-shape type=%q and nip=%d", faceType, nip)) {
		return nil, nil
	}
	return
}

func GetSolidFlags(extra string) (useB, debug bool, thickness float64) {

	// defaults
	useB = false
	debug = false
	thickness = 1.0

	// flag: use B matrix
	if s_useB, found := io.Keycode(extra, "useB"); found {
		useB = io.Atob(s_useB)
	}

	// fix useB flag in case of axisymmetric simulation
	if Global.Sim.Data.Axisym {
		useB = true
	}

	// flag: thickess => plane-stress
	if s_thick, found := io.Keycode(extra, "thick"); found {
		thickness = io.Atof(s_thick)
	}

	// fix thickness flag
	if !Global.Sim.Data.Pstress {
		thickness = 1.0
	}

	// flag: debug
	if s_debug, found := io.Keycode(extra, "debug"); found {
		debug = io.Atob(s_debug)
	}
	return
}

func GetSeepFaceFlags(extra string) (Macaulay bool, BetRamp, Kappa float64) {

	// defaults
	Macaulay = false
	BetRamp = math.Ln2 / 0.01
	Kappa = 1.0

	// use macaulay function ?
	if s_mac, found := io.Keycode(extra, "mac"); found {
		Macaulay = io.Atob(s_mac)
	}

	// coefficient for smooth ramp function
	if s_bet, found := io.Keycode(extra, "bet"); found {
		BetRamp = io.Atof(s_bet)
	}

	// κ coefficient
	if s_kap, found := io.Keycode(extra, "kap"); found {
		Kappa = io.Atof(s_kap)
	}
	return
}
