// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "github.com/cpmech/gosl/utl"

// errFunc returns a callback function to print error message
func errFunc(err error) func() {
	return func() {
		utl.Pfred("ERROR: %v\n", err)
	}
}

// max returns the max between two floats
func max(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

// min returns the min between two floats
func min(a, b float64) float64 {
	if a < b {
		return a
	}
	return b
}
