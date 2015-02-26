// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"
	"flag"

	"github.com/cpmech/gofem/fem"
	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
)

const IP_TAG_INI = 70000

var (
	ndim  int         // space dimension
	verts []*inp.Vert // all vertices
	cells []*inp.Cell // all cells
	nodes []*fem.Node // active/allocated nodes
	elems []fem.Elem  // active/allocated elements

	dirout  string // directory for output
	fnkey   string // filename key
	steady  bool   // steady simulation
	withips bool   // output integration points

	ukeys   = []string{"ux", "uy", "uz"}                      // displacement keys
	skeys   = []string{"sx", "sy", "sz", "sxy", "syz", "szx"} // stress keys
	rwlkeys = []string{"rwlx", "rwly", "rwlz"}                // ρl * wl keys
	plkeys  = []string{"pl"}                                  // liquid pressure keys
	pgkeys  = []string{"pg"}                                  // gas pressure keys

	is_sig     map[string]bool     // is sigma key? "sx" => true
	is_rwl     map[string]bool     // is rwl key? "rwlx" => true
	label2keys map[string][]string // maps, e.g., "u" => ukeys
)

func init() {
	is_sig = map[string]bool{"sx": true, "sy": true, "sz": true, "sxy": true, "syz": true, "szx": true}
	is_rwl = map[string]bool{"rwlx": true, "rwly": true, "rwlz": true}
	label2keys = map[string][]string{"u": ukeys, "sig": skeys, "rwl": rwlkeys, "pl": plkeys, "pg": pgkeys}
}

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// input data
	simfn := "data/twoqua4.sim"
	withips = true
	stgidx := 0

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		withips = io.Atob(flag.Arg(1))
	}
	if len(flag.Args()) > 2 {
		stgidx = io.Atoi(flag.Arg(2))
	}

	// check extension
	if io.FnExt(simfn) == "" {
		simfn += ".sim"
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  simfn   = %20s // simulation filename\n", simfn)
	io.Pf("  withips = %20v // with integration points\n", withips)
	io.Pf("  stgidx  = %20v // stage index\n", stgidx)
	io.Pf("\n")

	// start analysis process
	out.Start(simfn, stgidx, 0)

	// global variables
	ndim = out.Dom.Msh.Ndim
	verts = out.Dom.Msh.Verts
	cells = out.Dom.Msh.Cells
	nodes = out.Dom.Nodes
	elems = out.Dom.Elems
	dirout = fem.Global.Sim.Data.DirOut
	fnkey = fem.Global.Sim.Data.FnameKey
	steady = fem.Global.Sim.Data.Steady

	// flags
	has_sig := out.Ipkeys["sx"]
	has_rwl := out.Ipkeys["rwlx"]

	// buffers
	pvd := make(map[string]*bytes.Buffer)
	geo := make(map[string]*bytes.Buffer)
	vtu := make(map[string]*bytes.Buffer)
	if _, ok := out.Dom.YandC["ux"]; ok {
		pvd["u"] = new(bytes.Buffer)
		geo["u"] = new(bytes.Buffer)
		vtu["u"] = new(bytes.Buffer)
	}
	if _, ok := out.Dom.YandC["pl"]; ok {
		pvd["pl"] = new(bytes.Buffer)
		geo["pl"] = new(bytes.Buffer)
		vtu["pl"] = new(bytes.Buffer)
	}
	if _, ok := out.Dom.YandC["pg"]; ok {
		pvd["pg"] = new(bytes.Buffer)
		geo["pg"] = new(bytes.Buffer)
		vtu["pg"] = new(bytes.Buffer)
	}
	if len(out.Ipkeys) > 0 {
		pvd["ips"] = new(bytes.Buffer)
		geo["ips"] = new(bytes.Buffer)
		vtu["ips"] = new(bytes.Buffer)
	}

	// headers
	for _, b := range pvd {
		pvd_header(b)
	}

	// process results
	for tidx, t := range fem.Global.Sum.OutTimes {

		// input results into domain
		if !out.Dom.In(tidx) {
			chk.Panic("cannot load results into domain; please check log file")
		}

		// generate topology
		if tidx == 0 {
			for label, b := range geo {
				topology(b, label == "ips", label == "pl" || label == "pg")
			}
		}

		// for each data buffer
		for label, b := range vtu {

			// reset buffer
			b.Reset()

			// points data
			if label == "ips" {
				pdata_open(b)
				if has_sig {
					pdata_write(b, "sig", skeys, true)
				}
				if has_rwl {
					pdata_write(b, "rwl", rwlkeys, true)
				}
				for key, _ := range out.Ipkeys {
					if !is_sig[key] && !is_rwl[key] {
						pdata_write(b, key, []string{key}, true)
					}
				}
				pdata_close(b)
			} else {
				pdata_open(b)
				pdata_write(b, label, label2keys[label], false)
				pdata_close(b)
			}

			// cells data
			cdata_write(b, label == "ips")
		}

		/*
			// node's points-data: displacements
			if has_ux {
				pdata_write(b_dat_ge, "u", ukeys, out.Dom.Sol.Y, false, true)
				if !steady {
					pdata_write(b_dat_ge, "v", ukeys, out.Dom.Sol.Dydt, false, true)
					pdata_write(b_dat_ge, "a", ukeys, out.Dom.Sol.D2ydt2, false, true)
				}
			}
			// node's points-data: pressure
			if has_pl {
				pdata_write(b_dat_pl, "pl", []string{"pl"}, out.Dom.Sol.Y, false, true)
			}
			if has_pg {
				pdata_write(b_dat_pg, "pg", []string{"pl"}, out.Dom.Sol.Y, false, true)
			}
			// node's points-data: general scalars
			for key, _ := range out.Dom.YandC {
				if _, isnongen := nongeneral[key]; !isnongen {
					pdata_write(b_dat_ge, key, []string{key}, out.Dom.Sol.Y, false, true)
				}
			}
			// element's points-data: stresses
			if has_sx {
				pdata_write(b_dat_ge, "sig", skeys, out.Dom.Sol.Y, true, false)
			}
			// element's points-data: ρl * wl
			if has_rwlx {
				pdata_write(b_dat_ge, "rwl", rwlkeys, out.Dom.Sol.Y, true, false)
			}
			// element's points-data: general
			for key, _ := range out.Ipkeys {
				if _, isnongen := nongeneral[key]; !isnongen {
					pdata_write(b_dat_ge, key, []string{key}, out.Dom.Sol.Y, true, false)
				}
			}
			// close points-data section
			close_pdata(b_dat_ge)
			close_pdata(b_dat_pl)
			close_pdata(b_dat_pg)
		*/

		// vtu
		//vtu_write(b_top_pg, b_dat_pg, tidx, "_pg")

		// pvd
		for label, b := range pvd {
			pvd_line(b, tidx, t, label)
		}
	}

	// write pvd files
	for label, b := range pvd {
		pvd_write(b, label)
	}
}

// headers and footers ///////////////////////////////////////////////////////////////////////////////

func pvd_header(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n")
}

func pvd_line(buf *bytes.Buffer, tidx int, time float64, label string) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<DataSet timestep=\"%23.15e\" file=\"%s_%06d_%s.vtu\" />\n", time, fnkey, tidx, label)
}

func pvd_write(buf *bytes.Buffer, label string) {
	if buf == nil {
		return
	}
	io.Ff(buf, "</Collection>\n</VTKFile>")
	io.WriteFileV(io.Sf("%s/%s_%s.pvd", dirout, fnkey, label), buf)
}

func vtu_write(top, dat *bytes.Buffer, tidx int, suffix string) {
	if top == nil || dat == nil {
		return
	}
	nv := len(nodes)
	nc := len(elems)
	if withips {
		nv += len(out.Ipoints)
		nc += len(out.Ipoints)
	}
	var hdr, foo bytes.Buffer
	io.Ff(&hdr, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n")
	io.Ff(&hdr, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nv, nc)
	io.Ff(&foo, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")
	io.WriteFile(io.Sf("%s/%s_%06d%s.vtu", dirout, fnkey, tidx, suffix), &hdr, top, dat, &foo)
}

// topology ////////////////////////////////////////////////////////////////////////////////////////

func topology(buf *bytes.Buffer, ips, lbb bool) {
	if buf == nil {
		return
	}

	// coordinates
	io.Ff(buf, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	var z float64
	if ips {
		for _, p := range out.Ipoints {
			if ndim == 3 {
				z = p.X[2]
			}
			io.Ff(buf, "%23.15e %23.15e %23.15e ", p.X[0], p.X[1], z)
		}
	} else {
		for _, v := range verts {
			if ndim == 3 {
				z = v.C[2]
			}
			io.Ff(buf, "%23.15e %23.15e %23.15e ", v.C[0], v.C[1], z)
		}
	}
	io.Ff(buf, "\n</DataArray>\n</Points>\n")

	// connectivities
	io.Ff(buf, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	if ips {
		for i, _ := range out.Ipoints {
			io.Ff(buf, "%d ", i)
		}
	} else {
		for _, e := range elems {
			cell := cells[e.Id()]
			verts := cell.Verts
			nverts := len(verts)
			if lbb {
				nverts = shp.GetNverts(shp.GetBasicType(cell.Type))
			}
			for j := 0; j < nverts; j++ {
				io.Ff(buf, "%d ", verts[j])
			}
		}
	}

	// offsets of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	var offset int
	if ips {
		for range out.Ipoints {
			offset += 1
			io.Ff(buf, "%d ", offset)
		}
	} else {
		for _, e := range elems {
			cell := cells[e.Id()]
			nverts := len(cell.Verts)
			if lbb {
				nverts = shp.GetNverts(shp.GetBasicType(cell.Type))
			}
			if cell.Type == "qua9" {
				nverts = 8
			}
			offset += nverts
			io.Ff(buf, "%d ", offset)
		}
	}

	// types of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
	if ips {
		for range out.Ipoints {
			io.Ff(buf, "%d ", shp.VTK_VERTEX)
		}
	} else {
		for _, e := range elems {
			cell := cells[e.Id()]
			vtk := shp.VTK_QUADRATIC_QUAD // TODO: shp.GetVtk(cell.Type)
			if lbb {
				// TODO: vtk = shp.GetVtk(shp.GetBasicType(cell.Type))
				vtk = shp.VTK_QUADRATIC_QUAD // TODO: <<<< REMOVE THIS
			}
			if cell.Type == "qua9" {
				vtk = shp.VTK_QUADRATIC_QUAD
			}
			if vtk < 0 {
				chk.Panic("cannot handle cell type %q", cell.Type)
			}
			io.Ff(buf, "%d ", vtk)
		}
	}
	io.Ff(buf, "\n</DataArray>\n</Cells>\n")
	return
}

// points data /////////////////////////////////////////////////////////////////////////////////////

func pdata_open(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<PointData Scalars=\"TheScalars\">\n")
}

func pdata_close(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "</PointData>\n")
}

func get_zeros(n int) (l string) {
	for i := 0; i < n; i++ {
		l += "0 "
	}
	return
}

func pdata_write(buf *bytes.Buffer, label string, keys []string, ips bool) {
	if buf == nil {
		return
	}

	nkeys := len(keys)
	zeros := get_zeros(nkeys)

	io.Ff(buf, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", label, nkeys)
	if ips {
		// loop over integration points
		for _, p := range out.Ipoints {
			l := ""
			for _, key := range keys {
				if v, ok := p.V[key]; ok {
					io.Ff(buf, "%23.15e ", *v)
				} else {
					l += "0 "
				}
			}
			if l == "" {
				l = zeros
			}
			io.Ff(buf, l)
		}
	} else {
		// loop over vertices
		Y := out.Dom.Sol.Y
		switch label {
		case "v":
			Y = out.Dom.Sol.Dydt
		case "a":
			Y = out.Dom.Sol.D2ydt2
		}
		for _, v := range verts {
			n := out.Dom.Vid2node[v.Id]
			l := zeros
			if n != nil {
				l = ""
				for _, key := range keys {
					eq := n.GetEq(key)
					if eq >= 0 {
						l += io.Sf("%23.15e ", Y[eq])
					} else {
						l += "0 "
					}
				}
			}
			io.Ff(buf, l)
		}
	}
	io.Ff(buf, "\n</DataArray>\n")
}

func cdata_write(buf *bytes.Buffer, ips bool) {
	if buf == nil {
		return
	}

	// open
	io.Ff(buf, "<CellData Scalars=\"TheScalars\">\n")

	// ids
	io.Ff(buf, "<DataArray type=\"Float64\" Name=\"eid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	if ips {
		for _, p := range out.Ipoints {
			io.Ff(buf, "%d ", p.Eid)
		}
	} else {
		for _, e := range elems {
			io.Ff(buf, "%d ", e.Id())
		}
	}

	// cells positive tags
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Float64\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	if ips {
		for _, p := range out.Ipoints {
			ptag := IP_TAG_INI + iabs(cells[p.Eid].Tag)
			io.Ff(buf, "%d ", ptag)
		}
	} else {
		for _, e := range elems {
			ptag := iabs(cells[e.Id()].Tag)
			io.Ff(buf, "%d ", ptag)
		}
	}

	// close
	io.Ff(buf, "\n</DataArray>\n</CellData>\n")
}

func iabs(val int) int {
	if val < 0 {
		return -val
	}
	return val
}
