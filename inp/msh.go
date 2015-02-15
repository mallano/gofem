// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package inp

import (
	"encoding/json"
	"log"
	"math"

	"github.com/cpmech/gofem/shp"

	"github.com/cpmech/gosl/utl"
)

// constants
const Ztol = 1e-7

// Vert holds vertex data
type Vert struct {
	Id  int       // id
	Tag int       // tag
	C   []float64 // coordinates (size==2 or 3)
}

// Cell holds cell data
type Cell struct {
	Id     int    // id
	Tag    int    // tag
	Geo    int    // geometry type (gemlab code)
	Type   string // geometry type (string)
	Part   int    // partition id
	Verts  []int  // vertices
	FTags  []int  // edge (2D) or face (3D) tags
	STags  []int  // seam tags (for 3D only; it is actually a 3D edge tag)
	JlinId int    // joint line id
	JsldId int    // joint solid id

	// specific problems data
	SeepVerts map[int]bool // local vertices ids of vertices on seepage faces

	// derived
	Shp *shp.Shape // shape structure
}

// CellFaceId structure
type CellFaceId struct {
	C   *Cell // cell
	Fid int   // face id
}

// CellSeamId structure
type CellSeamId struct {
	C   *Cell // cell
	Sid int   // seam id
}

// Mesh holds a mesh for FE analyses
type Mesh struct {

	// from JSON
	Verts []*Vert // vertices
	Cells []*Cell // cells

	// derived
	Ndim       int     // space dimension
	Xmin, Xmax float64 // min and max x-coordinate
	Ymin, Ymax float64 // min and max x-coordinate
	Zmin, Zmax float64 // min and max x-coordinate

	// derived: maps
	VertTag2verts map[int][]*Vert      // vertex tag => set of vertices
	CellTag2cells map[int][]*Cell      // cell tag => set of cells
	FaceTag2cells map[int][]CellFaceId // face tag => set of cells
	SeamTag2cells map[int][]CellSeamId // seam tag => set of cells
	Ctype2cells   map[string][]*Cell   // cell type => set of cells
	Part2cells    map[int][]*Cell      // partition number => set of cells
}

// ReadMsh reads a mesh for FE analyses
//  Note: returns nil on errors
func ReadMsh(fn string) *Mesh {

	// new mesh
	var o Mesh

	// read file
	b, err := utl.ReadFile(fn)
	if LogErr(err, "msh: cannot open mesh file "+fn+"\n") {
		return nil
	}

	// decode
	if LogErr(json.Unmarshal(b, &o), "msh: cannot unmarshal mesh file "+fn+"\n") {
		return nil
	}

	// check
	if LogErrCond(len(o.Verts) < 2, "msh: mesh must have at least 2 vertices and 1 cell") {
		return nil
	}
	if LogErrCond(len(o.Cells) < 1, "msh: mesh must have at least 2 vertices and 1 cell") {
		return nil
	}

	// vertex related derived data
	o.Ndim = 2
	o.Xmin = o.Verts[0].C[0]
	o.Ymin = o.Verts[0].C[1]
	if len(o.Verts[0].C) > 2 {
		o.Zmin = o.Verts[0].C[2]
	}
	o.Xmax = o.Xmin
	o.Ymax = o.Ymin
	o.Zmax = o.Zmin
	o.VertTag2verts = make(map[int][]*Vert)
	for i, v := range o.Verts {

		// check vertex id
		if LogErrCond(v.Id != i, "msh: vertices must be sequentially numbered. %d != %d\n", v.Id, i) {
			return nil
		}

		// ndim
		nd := len(v.C)
		if LogErrCond(nd < 2 || nd > 4, "msh: ndim must be 2 or 3\n") {
			return nil
		}
		if nd == 3 {
			if math.Abs(v.C[2]) > Ztol {
				o.Ndim = 3
			}
		}

		// tags
		if v.Tag < 0 {
			verts := o.VertTag2verts[v.Tag]
			o.VertTag2verts[v.Tag] = append(verts, v)
		}

		// limits
		o.Xmin = min(o.Xmin, v.C[0])
		o.Xmax = max(o.Xmax, v.C[0])
		o.Ymin = min(o.Ymin, v.C[1])
		o.Ymax = max(o.Ymax, v.C[1])
		if nd > 2 {
			o.Zmin = min(o.Zmin, v.C[2])
			o.Zmax = max(o.Zmax, v.C[2])
		}
	}

	// derived data
	o.CellTag2cells = make(map[int][]*Cell)
	o.FaceTag2cells = make(map[int][]CellFaceId)
	o.SeamTag2cells = make(map[int][]CellSeamId)
	o.Ctype2cells = make(map[string][]*Cell)
	o.Part2cells = make(map[int][]*Cell)
	for i, c := range o.Cells {

		// shape
		c.Shp = shp.Get(c.Type)

		// check id and tag
		if LogErrCond(c.Id != i, "msh: cells must be sequentially numbered. %d != %d\n", c.Id, i) {
			return nil
		}
		if LogErrCond(c.Tag >= 0, "msh: cell tags must be negative\n") {
			return nil
		}

		// face tags
		cells := o.CellTag2cells[c.Tag]
		o.CellTag2cells[c.Tag] = append(cells, c)
		for i, ftag := range c.FTags {
			if ftag < 0 {
				pairs := o.FaceTag2cells[ftag]
				o.FaceTag2cells[ftag] = append(pairs, CellFaceId{c, i})
			}
		}

		// seam tags
		if c.Shp.Gndim == 3 {
			for i, stag := range c.STags {
				if stag < 0 {
					pairs := o.SeamTag2cells[stag]
					o.SeamTag2cells[stag] = append(pairs, CellSeamId{c, i})
				}
			}
		}

		// cell type => cells
		cells = o.Ctype2cells[c.Type]
		o.Ctype2cells[c.Type] = append(cells, c)

		// partition => cells
		cells = o.Part2cells[c.Part]
		o.Part2cells[c.Part] = append(cells, c)
	}

	// log
	log.Printf("msh: fn=%s nverts=%d ncells=%d ncelltags=%d nfacetags=%d nseamtags=%d nverttags=%d ncelltypes=%d npart=%d\n", fn, len(o.Verts), len(o.Cells), len(o.CellTag2cells), len(o.FaceTag2cells), len(o.SeamTag2cells), len(o.VertTag2verts), len(o.Ctype2cells), len(o.Part2cells))
	return &o
}

// String returns a JSON representation of *Vert
func (o *Vert) String() string {
	l := utl.Sf("{\"id\":%4d, \"tag\":%6d, \"c\":[", o.Id, o.Tag)
	for i, x := range o.C {
		if i > 0 {
			l += ", "
		}
		l += utl.Sf("%23.15e", x)
	}
	l += "] }"
	return l
}

// String returns a JSON representation of *Cell
func (o *Cell) String() string {
	l := utl.Sf("{\"id\":%d, \"tag\":%d, \"type\":%q, \"part\":%d, \"verts\":[", o.Id, o.Tag, o.Type, o.Part)
	for i, x := range o.Verts {
		if i > 0 {
			l += ", "
		}
		l += utl.Sf("%d", x)
	}
	l += "], \"ftags\":["
	for i, x := range o.FTags {
		if i > 0 {
			l += ", "
		}
		l += utl.Sf("%d", x)
	}
	l += "] }"
	return l
}

// String returns a JSON representation of *Mesh
func (o Mesh) String() string {
	l := "{\n  \"verts\" : [\n"
	for i, x := range o.Verts {
		if i > 0 {
			l += ",\n"
		}
		l += utl.Sf("    %v", x)
	}
	l += "\n  ],\n  \"cells\" : [\n"
	for i, x := range o.Cells {
		if i > 0 {
			l += ",\n"
		}
		l += utl.Sf("    %v", x)
	}
	l += "\n  ]\n}"
	return l
}
