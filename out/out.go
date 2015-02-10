// Copyright 2015 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package out

import "github.com/cpmech/gofem/fem"

type Extractor interface {
}

type TimeSlicer interface {
}

type SpaceSlicer interface {
}

type Node struct {
	fem.Node
}

type Ipoit struct {
}
