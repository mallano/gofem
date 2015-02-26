// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

func GetNverts(cellType string) int {
	shape, ok := factory[cellType]
	if !ok {
		return -1
	}
	return shape.Nverts
}

func GetFaceLocalVerts(cellType string, fidx int) []int {
	shape, ok := factory[cellType]
	if !ok {
		return nil
	}
	return shape.FaceLocalV[fidx]
}

func GetFaceType(cellType string) string {
	shape, ok := factory[cellType]
	if !ok {
		return ""
	}
	return shape.FaceType
}

func GetBasicType(cellType string) string {
	shape, ok := factory[cellType]
	if !ok {
		return ""
	}
	return shape.BasicType
}
