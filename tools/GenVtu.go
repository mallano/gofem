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
)

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

	// build map of non-general keys
	nongeneral := fem.GetIsEssenKeyMap()
	for _, key := range ukeys {
		nongeneral[key] = true
	}
	for _, key := range skeys {
		nongeneral[key] = true
	}
	for _, key := range rwlkeys {
		nongeneral[key] = true
	}
	nongeneral["pl"] = true
	nongeneral["pg"] = true

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

	// check for keys
	_, has_ux := out.Dom.YandC["ux"]
	_, has_pl := out.Dom.YandC["pl"]
	_, has_pg := out.Dom.YandC["pg"]
	_, has_sx := out.Ipkeys["sx"]
	_, has_rwlx := out.Ipkeys["rwlx"]

	// buffers for general scalars
	b_top_ge := new(bytes.Buffer) // points/general
	b_pvd_ge := new(bytes.Buffer) // pvd file/general

	// buffers for "pl"
	var b_top_pl *bytes.Buffer // points for "pl" file
	var b_pvd_pl *bytes.Buffer // pvd file for "pl" keys
	if has_pl {
		b_top_pl = new(bytes.Buffer)
		b_pvd_pl = new(bytes.Buffer)
	}

	// buffers for "pg"
	var b_top_pg *bytes.Buffer // points for "pg" file
	var b_pvd_pg *bytes.Buffer // pvd file for "pg" keys
	if has_pg {
		b_top_pg = new(bytes.Buffer)
		b_pvd_pg = new(bytes.Buffer)
	}

	// headers
	pvd_header(b_pvd_ge)
	pvd_header(b_pvd_pl)
	pvd_header(b_pvd_pg)

	// process results
	for tidx, t := range out.Sum.Times {

		// input results into domain
		if !out.Dom.In(tidx) {
			chk.Panic("cannot load results into domain; please check log file")
		}

		// generate topology
		if tidx == 0 {
			topology(b_top_ge)
			topology(b_top_pl)
			topology(b_top_pg)
		}

		// buffers for data
		var b_dat_ge, b_dat_pl, b_dat_pg *bytes.Buffer
		b_dat_ge = new(bytes.Buffer)
		if has_pl {
			b_dat_pl = new(bytes.Buffer)
		}
		if has_pg {
			b_dat_pg = new(bytes.Buffer)
		}

		// open points-data section
		open_pdata(b_dat_ge)
		open_pdata(b_dat_pl)
		open_pdata(b_dat_pg)

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

		// tags
		cdata_ids_and_tags(b_dat_ge, false, false)
		cdata_ids_and_tags(b_dat_pl, false, false)
		cdata_ids_and_tags(b_dat_pg, false, false)

		// vtu
		vtu_write(b_top_ge, b_dat_ge, tidx, "")
		vtu_write(b_top_pl, b_dat_pl, tidx, "_pl")
		vtu_write(b_top_pg, b_dat_pg, tidx, "_pg")

		// pvd
		pvd_line(b_pvd_ge, tidx, t, "")
		pvd_line(b_pvd_pl, tidx, t, "_pl")
		pvd_line(b_pvd_pg, tidx, t, "_pg")
	}

	// write pvd file
	pvd_write(b_pvd_ge, "")
	pvd_write(b_pvd_pl, "_pl")
	pvd_write(b_pvd_pg, "_pg")
}

// headers and footers ///////////////////////////////////////////////////////////////////////////////

func pvd_header(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n")
}

func pvd_line(buf *bytes.Buffer, tidx int, time float64, suffix string) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<DataSet timestep=\"%23.15e\" file=\"%s_%06d%s.vtu\" />\n", time, fnkey, tidx, suffix)
}

func pvd_write(buf *bytes.Buffer, suffix string) {
	if buf == nil {
		return
	}
	io.Ff(buf, "</Collection>\n</VTKFile>")
	io.WriteFileV(io.Sf("%s/%s%s.pvd", dirout, fnkey, suffix), buf)
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

func topology(buf *bytes.Buffer) {
	if buf == nil {
		return
	}

	// all nodes coordinates
	io.Ff(buf, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	for _, v := range verts {
		if ndim == 2 {
			io.Ff(buf, "%23.15e %23.15e 0 ", v.C[0], v.C[1])
		} else {
			io.Ff(buf, "%23.15e %23.15e %23.15e ", v.C[0], v.C[1], v.C[2])
		}
	}
	if withips {
		for _, p := range out.Ipoints {
			if ndim == 2 {
				io.Ff(buf, "%23.15e %23.15e 0 ", p.X[0], p.X[1])
			} else {
				io.Ff(buf, "%23.15e %23.15e %23.15e ", p.X[0], p.X[1], p.X[2])
			}
		}
	}
	io.Ff(buf, "\n</DataArray>\n</Points>\n")

	// connectivity of active elements
	io.Ff(buf, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	for _, e := range elems {
		for _, n := range cells[e.Id()].Verts {
			io.Ff(buf, "%d ", n)
		}
	}
	if withips {
		npt := len(verts)
		for range out.Ipoints {
			io.Ff(buf, "%d ", npt)
			npt += 1
		}
	}

	// offsets of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	var offset int
	for _, e := range elems {
		offset += len(cells[e.Id()].Verts)
		io.Ff(buf, "%d ", offset)
	}
	if withips {
		for range out.Ipoints {
			offset += 1
			io.Ff(buf, "%d ", offset)
		}
	}

	// types of active elements
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
	for _, e := range elems {
		c := cells[e.Id()]
		vtk := c.Shp.VtkCode
		if vtk == shp.VTK_POLY_VERTEX {
			switch c.Type {
			case "qua9":
				vtk = shp.VTK_QUADRATIC_QUAD
			default:
				chk.Panic("cannot handle cell type %q", c.Type)
			}
		}
		io.Ff(buf, "%d ", vtk)
	}
	if withips {
		for range out.Ipoints {
			io.Ff(buf, "%d ", shp.VTK_VERTEX)
		}
	}
	io.Ff(buf, "\n</DataArray>\n</Cells>\n")
	return
}

// points data /////////////////////////////////////////////////////////////////////////////////////

func open_pdata(buf *bytes.Buffer) {
	if buf == nil {
		return
	}
	io.Ff(buf, "<PointData Scalars=\"TheScalars\">\n")
}

func close_pdata(buf *bytes.Buffer) {
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

func pdata_write(buf *bytes.Buffer, label string, keys []string, Y []float64, skipNodes, skipIps bool) {
	if buf == nil {
		return
	}

	nkeys := len(keys)
	zeros := get_zeros(nkeys)

	io.Ff(buf, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"ascii\">\n", label, nkeys)

	// loop over nodes
	for _, v := range verts {
		n := out.Dom.Vid2node[v.Id]
		l := zeros
		if n != nil && !skipNodes {
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

	// loop over integration points
	if withips {
		for _, p := range out.Ipoints {
			l := zeros
			if !skipIps {
				l = ""
				for _, key := range keys {
					if v, ok := p.V[key]; ok {
						io.Ff(buf, "%23.15e ", *v)
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

func cdata_ids_and_tags(buf *bytes.Buffer, skipNodes, skipIps bool) {
	if buf == nil {
		return
	}

	io.Ff(buf, "<CellData Scalars=\"TheScalars\">\n")

	// ids
	io.Ff(buf, "<DataArray type=\"Float64\" Name=\"eid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, e := range elems {
		io.Ff(buf, "%d ", e.Id())
	}
	if withips {
		for _, p := range out.Ipoints {
			io.Ff(buf, "%d ", p.Eid)
		}
	}

	// cells positive tags
	io.Ff(buf, "\n</DataArray>\n<DataArray type=\"Float64\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, e := range elems {
		ptag := iabs(cells[e.Id()].Tag)
		io.Ff(buf, "%d ", ptag)
	}
	if withips {
		for _, p := range out.Ipoints {
			ptag := IP_TAG_INI + iabs(cells[p.Eid].Tag)
			io.Ff(buf, "%d ", ptag)
		}
	}
	io.Ff(buf, "\n</DataArray>\n</CellData>\n")
}

func iabs(val int) int {
	if val < 0 {
		return -val
	}
	return val
}
