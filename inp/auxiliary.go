// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

// imax returns the max between two integers
func imax(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// imin returns the min between two integers
func imin(a, b int) int {
	if a < b {
		return a
	}
	return b
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
