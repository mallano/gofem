// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "github.com/cpmech/gosl/plt"

// Styles
type Styles []plt.LineData

func GetDefaultStyles(n int) Styles {
	return make([]plt.LineData, n)
}
