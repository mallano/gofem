// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package fem

import "github.com/cpmech/gosl/fun"

// NaturalBc holds information on natural boundary conditioins such as
// distributed loads or fluxes acting on surfaces
type NaturalBc struct {
	Key     string   // key such as qn, qn0, ql, seepH, seepP, etc...
	IdxFace int      // local index of face
	Fcn     fun.Func // function callback
	Extra   string   // extra information
}
