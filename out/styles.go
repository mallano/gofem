// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import (
	"github.com/cpmech/gosl/plt"
	"github.com/cpmech/gosl/utl"
)

// Styles
type Styles []plt.LineData

func GetDefaultStyles(qts Quantities) Styles {
	sty := make([]plt.LineData, len(qts))
	for i, q := range qts {
		sty[i].Label = utl.Sf("x=%v", q.X)
	}
	return sty
}

func GetTexLabel(key, unit string) string {
	l := "$"
	switch key {
	case "time":
		l += "t"
	case "ux":
		l += "u_x"
	case "uy":
		l += "u_y"
	case "uz":
		l += "u_z"
	case "pl":
		l += "p_{\\ell}"
	case "pg":
		l += "p_g"
	case "sx":
		l += "\\sigma_x"
	case "sy":
		l += "\\sigma_y"
	case "sz":
		l += "\\sigma_z"
	case "sxy":
		l += "\\sigma_{xy}"
	case "syz":
		l += "\\sigma_{yz}"
	case "szx":
		l += "\\sigma_{zx}"
	default:
		l += "UNKNOWN_LABEL"
	}
	if unit != "" {
		l += "\\;" + unit
	}
	l += "$"
	return l
}
