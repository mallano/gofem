// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build ignore

package main

import (
	"bytes"
	"flag"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gofem/out"
	"github.com/cpmech/gofem/shp"
	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/io"
	"github.com/cpmech/gosl/utl"
)

const IP_TAG0 = 70000

func main() {

	// finalise analysis process and catch errors
	defer out.End()

	// input data
	simfn := "data/twoqua4.sim"
	withIps := false

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		simfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		withIps = utl.Atob(flag.Arg(1))
	}

	// print input data
	io.Pf("\nInput data\n")
	io.Pf("==========\n")
	io.Pf("  matOld   = %20s // old material filename\n", matOld)
	io.Pf("  matNew   = %20s // new material filename\n", matNew)
	io.Pf("  convSymb = %20v // do convert symbols\n", convSymb)
	io.Pf("\n")

	// start analysis process
	out.Start(simfn, 0, 0)

	// filename key
	fnk := utl.FnKey(simfn)

	// buffers
	var pts, cls bytes.Buffer            // points and cells buffers
	var pts_pl, cls_pl bytes.Buffer      // points and cells for "pl" file
	var pts_pg, cls_pg bytes.Buffer      // points and cells for "pg" file
	var pvd, pvd_pl, pvd_pg bytes.Buffer // pvd files

	// total number of integration points
	var nipsAll int

	// process results
	outidx := 0
	for tidx, t := range out.Sum.Times {

		// generate topology
		if tidx == 0 {
			pts.Reset()
			cls.Reset()
			nipsAll = topology(&pts, &cls, g.M, r.A2e, withIps, r.E)
			if r.Lbb {
				pts_p.Reset()
				cls_p.Reset()
				topology_lbb(&pts_p, &cls_p, g.M, r.A2e)
			}
		}

		// header
		var hdr bytes.Buffer
		vtu_header(&hdr, g.Nv+nipsAll, len(r.A2e)+nipsAll)

		// open points-data section
		var dat, dat_p bytes.Buffer
		open_pdata(&dat, &dat_p, r.Lbb)

		// points-data: displacements
		if r.Has_ux {
			pdata_vector(&dat, "", "u", func(n int) int { return r.N[n].Ux_start }, g.M.Ndim, r.N, false, r.A2e, withIps, false, r.E)
		}

		// points-data: pressure
		if r.Lbb {
			pdata_scalar(&dat_p, "", "pl", r.N, r.A2e, withIps, false, r.E)
		}

		// points-data: scalars
		for _, key := range r.ScalarKeys {
			if r.Lbb && key == "pl" {
				continue
			}
			//io.Ffpink("scalar field = %v\n", key)
			pdata_scalar(&dat, "", key, r.N, r.A2e, withIps, false, r.E)
		}

		// points-data: @ integration points
		if withIps {
			allKeys := make(map[string]bool)
			for _, e := range r.A2e {
				for _, key := range r.E[e].Keys {
					allKeys[key] = true
				}
			}
			if _, ok := allKeys["rwlx"]; ok {
				pdata_vector(&dat, "ip_", "rwl", nil, g.M.Ndim, r.N, true, r.A2e, withIps, true, r.E)
			}
			ak := utl.StrBoolMapSort(allKeys)
			for _, key := range ak {
				if key != "rwlx" && key != "rwly" && key != "rwlz" {
					pdata_scalar(&dat, "ip_", key, r.N, r.A2e, withIps, true, r.E)
				}
			}
		}

		// close points-data section
		close_pdata(&dat, &dat_p, r.Lbb)

		// cells-data section
		gen_cells_data(&dat, g.M, r.A2e, withIps, r.E)

		// footer
		var foo bytes.Buffer
		vtu_footer(&foo)

		// write vtu file
		utl.WriteFile(io.Sf("vtufiles/%s_%06d.vtu", fnk, outidx), &hdr, &pts, &cls, &dat, &foo)
		if r.Lbb {
			utl.WriteFile(io.Sf("vtufiles/%s_%06d_p.vtu", fnk, outidx), &hdr, &pts_p, &cls_p, &dat_p, &foo)
			haslbb = true
		}

		// pvd line
		pvd_line(&pvd, &pvd_pl, r.Lbb, time, fnk, outidx)

		// next outidx
		outidx += 1
		return false
	}

	// write pvd file
	pvd_footer(&pvd, &pvd_pl, haslbb)
	utl.WriteFileV(io.Sf("vtufiles/%s.pvd", fnk), &pvd)
	if haslbb {
		utl.WriteFileV(io.Sf("vtufiles/%s_p.pvd", fnk), &pvd_pl)
	}
}

// headers and footers ///////////////////////////////////////////////////////////////////////////////

func pvd_header(pvd *bytes.Buffer) {
	io.Ff(pvd, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n")
}

func pvd_line(pvd *bytes.Buffer, time float64, fnk string, outidx int) {
	io.Ff(pvd, "<DataSet timestep=\"%23.15e\" file=\"%s_%06d.vtu\" />\n", time, fnk, outidx)
}

func pvd_footer(pvd *bytes.Buffer) {
	io.Ff(pvd, "</Collection>\n</VTKFile>")
}

func vtu_header(hdr *bytes.Buffer, nv, nc int) {
	io.Ff(hdr, "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n")
	io.Ff(hdr, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nv, nc)
}

func vtu_footer(foo *bytes.Buffer) {
	io.Ff(foo, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n")
}

// topology //////////////////////////////////////////////////////////////////////////////////////////

func topology(pts, cls *bytes.Buffer, withIps bool) (totalNips int) {

	// auxiliary
	msh := out.Dom.Msh
	verts := msh.Verts
	cells := msh.Cells
	elems := out.Dom.Elems

	// all nodes coordinates
	io.Ff(pts, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	for _, v := range verts {
		if ndim == 2 {
			io.Ff(pts, "%23.15e %23.15e 0 ", v.C[0], v.C[1])
		} else {
			io.Ff(pts, "%23.15e %23.15e %23.15e ", v.C[0], v.C[1], v.C[2])
		}
	}
	if withIps {
		for _, e := range elems {
			ipids := out.Cid2ips[e.Id()]
			for _, ipid := range ipids {
				ip := out.Ipoints[ipid]
				if ndim == 2 {
					io.Ff(pts, "%23.15e %23.15e 0 ", ip.X[0], ip.X[1])
				} else {
					io.Ff(pts, "%23.15e %23.15e %23.15e ", ip.X[0], ip.X[1], ip.X[2])
				}
				totalNips += 1
			}
		}
	}
	io.Ff(pts, "\n</DataArray>\n</Points>\n")

	// connectivity of active elements
	io.Ff(cls, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	for _, e := range elems {
		for _, n := range cells[e.Id()].Verts {
			io.Ff(cls, "%d ", n)
		}
	}
	if withIps {
		npt := len(verts)
		for _, e := range elems {
			ipids := out.Cid2ips[e.Id()]
			for range ipids {
				io.Ff(cls, "%d ", npt)
				npt += 1
			}
		}
	}

	// offsets of active elements
	io.Ff(cls, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	var offset int
	for _, e := range elems {
		offset += len(cells[e.Id()].Verts)
		io.Ff(cls, "%d ", offset)
	}
	if withIps {
		for _, e := range elems {
			ipids := out.Cid2ips[e.Id()]
			for range ipids {
				offset += 1
				io.Ff(cls, "%d ", offset)
			}
		}
	}

	// types of active elements
	io.Ff(cls, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
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
		io.Ff(cls, "%d ", vtk)
	}
	if withIps {
		for _, e := range elems {
			ipids := out.Cid2ips[e.Id()]
			for range ipids {
				io.Ff(cls, "%d ", shp.VTK_VERTEX)
			}
		}
	}
	io.Ff(cls, "\n</DataArray>\n</Cells>\n")
	return
}

// points data ///////////////////////////////////////////////////////////////////////////////////////

func open_pdata(dat *bytes.Buffer) {
	io.Ff(dat, "<PointData Scalars=\"TheScalars\">\n")
}

func close_pdata(dat *bytes.Buffer) {
	io.Ff(dat, "</PointData>\n")
}

func pdata_scalar(dat *bytes.Buffer, prefix, key string, withIps, skipNodes bool) {

	// auxiliary
	msh := out.Dom.Msh
	verts := msh.Verts
	cells := msh.Cells
	nodes := out.Dom.Vid2node
	elems := out.Dom.Elems

	// solution, dof-key (ukey) and number of time derivatives
	ukey, nderiv := out.ParseKey(key)
	Y := out.Dom.Sol.Y
	switch nderiv {
	case 2:
		Y = out.Dom.Sol.D2ydt2
	case 1:
		Y = out.Dom.Sol.Dydt
	}

	io.Ff(dat, "<DataArray type=\"Float64\" Name=\"%s%s\" NumberOfComponents=\"1\" format=\"ascii\">\n", prefix, key)

	// loop over nodes
	for _, n := range nodes {
		l := "0 "
		if n != nil && !skipNodes { // active and not to skip
			eq := n.GetEq(ukey)
			if eq >= 0 {
				l = io.Sf("%23.15e ", Y[eq])
			}
		}
		io.Ff(dat, "l")
	}

	// loop over integration points
	if withIps {
		for _, e := range elems {
			ipids := out.Cid2ips[e.Id()]
			if len(ipids) > 0 {
				if v, ok := out.Ipoints[ipids[0]].V[key]; ok {
					io.Ff(dat, "%23.15e ", *v)
				}
			} else {
				for i := 0; i < len(E[e].IpsC); i++ {
					io.Ff(dat, "0 ")
				}
			}
		}
	}
	io.Ff(dat, "\n</DataArray>\n")
}

func pdata_vector(dat *bytes.Buffer, prefix, key string, start func(n int) int, ndim int, nodes []out.Node, is_extrap bool, withIps, skipNodes bool, E []out.Elem) {
	io.Ff(dat, "<DataArray type=\"Float64\" Name=\"%s%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", prefix, key)
	v := []float64{0, 0, 0}
	for n, _ := range nodes {
		if skipNodes {
			io.Ff(dat, "0 0 0  ")
			continue
		}
		if nodes[n].Active {
			for i := 0; i < ndim; i++ {
				v[i] = nodes[n].Vals[start(n)+i]
				if is_extrap {
					v[i] /= nodes[n].Count[start(n)+i]
				}
			}
		}
		io.Ff(dat, "%23.15e %23.15e %23.15e  ", v[0], v[1], v[2])
	}
	if withIps {
		if key == "rwl" {
			if ndim == 2 {
				for _, e := range elems {
					i0, has0 := E[e].Key2idx["rwlx"]
					i1, has1 := E[e].Key2idx["rwly"]
					if has0 && has1 {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "%23.15e %23.15e 0  ", E[e].Vals[i][i0], E[e].Vals[i][i1])
						}
					} else {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "0 0 0  ")
						}
					}
				}
			} else {
				for _, e := range elems {
					i0, has0 := E[e].Key2idx["rwlx"]
					i1, has1 := E[e].Key2idx["rwly"]
					i2, has2 := E[e].Key2idx["rwlz"]
					if has0 && has1 && has2 {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "%23.15e %23.15e %23.15e  ", E[e].Vals[i][i0], E[e].Vals[i][i1], E[e].Vals[i][i2])
						}
					} else {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "0 0 0  ")
						}
					}
				}
			}
		} else {
			for _, e := range elems {
				for i := 0; i < len(E[e].IpsC); i++ {
					io.Ff(dat, "0 0 0  ")
				}
			}
		}
	}
	io.Ff(dat, "\n</DataArray>\n")
}

func pdata_tensor(dat *bytes.Buffer, prefix, key string, start func(n int) int, ndim int, nodes []out.Node, is_extrap bool, withIps, skipNodes bool, E []out.Elem) {
	io.Ff(dat, "<DataArray type=\"Float64\" Name=\"%s%s\" NumberOfComponents=\"6\" format=\"ascii\">\n", prefix, key)
	v := []float64{0, 0, 0, 0, 0, 0}
	for n, _ := range nodes {
		if skipNodes {
			io.Ff(dat, "0 0 0 0 0 0  ")
			continue
		}
		if nodes[n].Active {
			for i := 0; i < 2*ndim; i++ {
				v[i] = nodes[n].Vals[start(n)+i]
				if is_extrap {
					v[i] /= nodes[n].Count[start(n)+i]
				}
			}
		}
		io.Ff(dat, "%23.15e %23.15e %23.15e %23.15e %23.15e %23.15e  ", v[0], v[1], v[2], v[3], v[4], v[5])
	}
	if withIps {
		if key == "sE" {
			if ndim == 2 {
				for _, e := range elems {
					i0, has0 := E[e].Key2idx["sxE"]
					i1, has1 := E[e].Key2idx["syE"]
					i2, has2 := E[e].Key2idx["szE"]
					i3, has3 := E[e].Key2idx["sxyE"]
					if has0 && has1 && has2 && has3 {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "%23.15e %23.15e %23.15e %23.15e 0 0  ", E[e].Vals[i][i0], E[e].Vals[i][i1], E[e].Vals[i][i2], E[e].Vals[i][i3])
						}
					} else {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "0 0 0 0 0 0  ")
						}
					}
				}
			} else {
				for _, e := range elems {
					i0, has0 := E[e].Key2idx["sxE"]
					i1, has1 := E[e].Key2idx["syE"]
					i2, has2 := E[e].Key2idx["szE"]
					i3, has3 := E[e].Key2idx["sxyE"]
					i4, has4 := E[e].Key2idx["syzE"]
					i5, has5 := E[e].Key2idx["szxE"]
					if has0 && has1 && has2 && has3 && has4 && has5 {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "%23.15e %23.15e %23.15e %23.15e %23.15e %23.15e  ", E[e].Vals[i][i0], E[e].Vals[i][i1], E[e].Vals[i][i2], E[e].Vals[i][i3], E[e].Vals[i][i4], E[e].Vals[i][i5])
						}
					} else {
						for i := 0; i < len(E[e].IpsC); i++ {
							io.Ff(dat, "0 0 0 0 0 0  ")
						}
					}
				}
			}
		} else {
			for _, e := range elems {
				for i := 0; i < len(E[e].IpsC); i++ {
					io.Ff(dat, "0 0 0 0 0 0  ")
				}
			}
		}
	}
	io.Ff(dat, "\n</DataArray>\n")
}

// cells data ////////////////////////////////////////////////////////////////////////////////////////

func gen_cells_data(dat *bytes.Buffer, msh *inp.Mesh, withIps bool, E []out.Elem) {
	// cells positive tags and geo types
	io.Ff(dat, "<CellData Scalars=\"TheScalars\">\n")
	io.Ff(dat, "<DataArray type=\"Float64\" Name=\"tag\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, e := range elems {
		ptag := 0
		if cells[e].Tag < 0 {
			ptag = -cells[e].Tag
		}
		io.Ff(dat, "%d ", ptag)
	}
	if withIps {
		for _, e := range elems {
			for i := 0; i < len(E[e].IpsC); i++ {
				ptag := 0
				if cells[e].Tag < 0 {
					ptag = -cells[e].Tag
				}
				io.Ff(dat, "%d ", IP_TAG0+ptag)
			}
		}
	}
	io.Ff(dat, "\n</DataArray>\n<DataArray type=\"Float64\" Name=\"geo\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, e := range elems {
		io.Ff(dat, "%d ", cells[e].Geo)
	}
	if withIps {
		for _, e := range elems {
			for i := 0; i < len(E[e].IpsC); i++ {
				io.Ff(dat, "-1 ")
			}
		}
	}
	io.Ff(dat, "\n</DataArray>\n<DataArray type=\"Float64\" Name=\"eid\" NumberOfComponents=\"1\" format=\"ascii\">\n")
	for _, e := range elems {
		io.Ff(dat, "%d ", cells[e].Id)
	}
	if withIps {
		ncl := len(cells)
		for _, e := range elems {
			for i := 0; i < len(E[e].IpsC); i++ {
				io.Ff(dat, "%d ", ncl)
				ncl += 1
			}
		}
	}
	io.Ff(dat, "\n</DataArray>\n</CellData>\n")
}

func topology_lbb(pts, cls *bytes.Buffer, msh *inp.Mesh, elems []int) {

	// all nodes coordinates
	io.Ff(pts, "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n")
	for n, _ := range verts {
		if ndim == 2 {
			io.Ff(pts, "%23.15e %23.15e 0 ", verts[n].C[0], verts[n].C[1])
		} else {
			io.Ff(pts, "%23.15e %23.15e %23.15e ", verts[n].C[0], verts[n].C[1], verts[n].C[2])
		}
	}
	io.Ff(pts, "\n</DataArray>\n</Points>\n")

	// connectivity of active elements; but only getting LBB nodes
	io.Ff(cls, "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
	for _, e := range elems {
		for i := 0; i < inp.Geo2nvc[inp.Geo2geob[cells[e].Geo]]; i++ {
			io.Ff(cls, "%d ", cells[e].Verts[i])
		}
	}

	// offsets of active elements
	io.Ff(cls, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
	var offset int
	for _, e := range elems {
		offset += inp.Geo2nvc[inp.Geo2geob[cells[e].Geo]]
		io.Ff(cls, "%d ", offset)
	}

	// types of active elements
	io.Ff(cls, "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
	for _, e := range elems {
		io.Ff(cls, "%d ", inp.Geo2vtk[inp.Geo2geob[cells[e].Geo]])
	}
	io.Ff(cls, "\n</DataArray>\n</Cells>\n")
	return
}
