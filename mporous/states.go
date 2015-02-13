// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mporous

// StateLG holds state variables for porous media with liquid and gas
type StateLG struct {
	Sl   float64 // sl: liquid saturation
	RhoL float64 // ρL: real (intrinsic) density of liquid
	RhoG float64 // ρG: real (intrinsic) density of gas
	Ns0  float64 // ns0: initial partial fraction of solids
	Wet  bool    // wetting flag
	Dpc  float64 // Δpc: step increment of capillary pressure
}
