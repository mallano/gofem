// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

// ParseKey parses key like "duxdt" returning "ux" and time derivative number
//  Output: {key, number-of-time-derivatives}
//  Examples:  "ux"      => "ux", 0
//             "duxdt"   => "ux", 1
//             "d2uxdt2" => "ux", 2
func ParseKey(key string) (string, int) {
	if len(key) > 3 {
		n := len(key)
		if key[:1] == "d" && key[n-2:] == "dt" {
			return key[1 : n-2], 1
		}
		if len(key) > 5 {
			if key[:2] == "d2" && key[n-3:] == "dt2" {
				return key[2 : n-3], 2
			}
		}
	}
	return key, 0
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
