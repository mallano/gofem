// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package mporous implements models for porous media based on the Theory of Porous Media
package mporous

import "github.com/cpmech/gosl/fun"

// LGcond defines liquid-gas conductivity models
type LGcond interface {
	Init(prms fun.Prms) error   // Init initialises this structure
	GetPrms() fun.Prms          // gets (an example) of parameters
	Klr(sl float64) float64     // Klr returns klr
	Kgr(sl float64) float64     // Kgr returns kgr
	DklrDsl(sl float64) float64 // DklrDsl returns ∂klr/∂sl
	DkgrDsl(sl float64) float64 // DkgrDsl returns ∂kgr/∂sl
}
