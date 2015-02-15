// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// package shp implements shape structures/routines
package shp

import "github.com/cpmech/gosl/la"

// constants
const MINDET = 1.0e-14 // minimum determinant allowed for dxdR

// ShpFunc is the shape functions callback function
type ShpFunc func(S []float64, dSdR [][]float64, r, s, t float64, derivs bool)

// Shape holds geometry data
type Shape struct {

	// geometry
	Type       string      // name; e.g. "lin2"
	Func       ShpFunc     // shape/derivs function callback function
	FaceFunc   ShpFunc     // face shape/derivs function callback function
	BasicType  string      // geometry of basic element; e.g. "qua8" => "qua4"
	FaceType   string      // geometry of face; e.g. "qua8" => "lin3"
	Gndim      int         // geometry of shape; e.g. "lin3" => gnd == 1 (even in 3D simulations)
	Nverts     int         // number of vertices in cell; e.g. "qua8" => 8
	VtkCode    int         // VTK code
	FaceNverts int         // number of vertices on face
	FaceLocalV [][]int     // face local vertices [nfaces][FaceNverts]
	NatCoords  [][]float64 // natural coordinates [gndim][nverts]

	// geometry: for seams (3D-edges)
	SeamType   int     // geometry of seam (3D-edge); e.g. "hex8" => "lin2"
	SeamLocalV [][]int // seam (3d-edge) local vertices [nseams][nVertsOnSeam]

	// scratchpad: volume
	S    []float64   // [nverts] shape functions
	G    [][]float64 // [nverts][gndim] G == dSdx. derivative of shape function
	J    float64     // Jacobian: determinant of dxdr
	dSdR [][]float64 // [nverts][gndim] derivatives of S w.r.t natural coordinates
	dxdR [][]float64 // [gndim][gndim] derivatives of real coordinates w.r.t natural coordinates
	dRdx [][]float64 // [gndim][gndim] dRdx == inverse(dxdR)

	// scratchpad: line
	Jvec3d []float64 // Jacobian: norm of dxdr for line elements (size==3)
	Gvec   []float64 // [nverts] G == dSdx. derivative of shape function

	// scratchpad: face
	Sf     []float64   // [facenverts] shape functions values
	Fnvec  []float64   // [gndim] face normal vector multiplied by Jf
	dSfdRf [][]float64 // [facenverts][gndim-1] derivatives of Sf w.r.t natural coordinates
	dxfdRf [][]float64 // [gndim][gndim-1] derivatives of real coordinates w.r.t natural coordinates
}

// factory holds all Shapes available
var factory = make(map[string]*Shape)

// Get returns an existent Shape structure
//  Note: returns nil on errors
func Get(geoType string) *Shape {
	s, ok := factory[geoType]
	if !ok {
		return nil
	}
	return s
}

// CalcAtIp calculates volume data such as S and G at natural coordinate r
//  Input:
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//   ip                -- integration point
//  Output:
//   S, DSdR, DxdR, DRdx, G, and J
func (o *Shape) CalcAtIp(x [][]float64, ip *Ipoint, derivs bool) (err error) {

	// S and dSdR
	o.Func(o.S, o.dSdR, ip.R, ip.S, ip.T, derivs)
	if !derivs {
		return
	}

	if o.Gndim == 1 {
		// calculate Jvec3d == dxdR
		for i := 0; i < len(x); i++ {
			o.Jvec3d[i] = 0.0
			for m := 0; m < o.Nverts; m++ {
				o.Jvec3d[i] += x[i][m] * o.dSdR[m][0] // dxdR := x * dSdR
			}
		}

		// calculate J = norm of Jvec3d
		o.J = la.VecNorm(o.Jvec3d)

		// calculate G
		for m := 0; m < o.Nverts; m++ {
			o.Gvec[m] = o.dSdR[m][0] / o.J
		}

		return
	}

	// dxdR := sum_n x * dSdR   =>  dx_i/dR_j := sum_n x^n_i * dS^n/dR_j
	for i := 0; i < len(x); i++ {
		for j := 0; j < o.Gndim; j++ {
			o.dxdR[i][j] = 0.0
			for n := 0; n < o.Nverts; n++ {
				o.dxdR[i][j] += x[i][n] * o.dSdR[n][j]
			}
		}
	}

	// dRdx := inv(dxdR)
	o.J, err = la.MatInv(o.dRdx, o.dxdR, MINDET)
	if err != nil {
		return
	}

	// G == dSdx := dSdR * dRdx  =>  dS^m/dR_i := sum_i dS^m/dR_i * dR_i/dx_j
	la.MatMul(o.G, 1, o.dSdR, o.dRdx)
	return
}

// CalcAtR calculates volume data such as S and G at natural coordinate r
//  Input:
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//   R                 -- local/natural coordinates
//  Output:
//   S, DSdR, DxdR, DRdx, G, and J
func (o *Shape) CalcAtR(x [][]float64, R []float64, derivs bool) (err error) {

	r := R[0]
	s := R[1]
	t := 0.0
	if len(R) == 3 {
		t = R[2]
	}

	// S and dSdR
	o.Func(o.S, o.dSdR, r, s, t, derivs)
	if !derivs {
		return
	}

	if o.Gndim == 1 {
		// calculate Jvec3d == dxdR
		for i := 0; i < len(x); i++ {
			o.Jvec3d[i] = 0.0
			for m := 0; m < o.Nverts; m++ {
				o.Jvec3d[i] += x[i][m] * o.dSdR[m][0] // dxdR := x * dSdR
			}
		}

		// calculate J = norm of Jvec3d
		o.J = la.VecNorm(o.Jvec3d)

		// calculate G
		for m := 0; m < o.Nverts; m++ {
			o.Gvec[m] = o.dSdR[m][0] / o.J
		}

		return
	}

	// dxdR := sum_n x * dSdR   =>  dx_i/dR_j := sum_n x^n_i * dS^n/dR_j
	for i := 0; i < len(x); i++ {
		for j := 0; j < o.Gndim; j++ {
			o.dxdR[i][j] = 0.0
			for n := 0; n < o.Nverts; n++ {
				o.dxdR[i][j] += x[i][n] * o.dSdR[n][j]
			}
		}
	}

	// dRdx := inv(dxdR)
	o.J, err = la.MatInv(o.dRdx, o.dxdR, MINDET)
	if err != nil {
		return
	}

	// G == dSdx := dSdR * dRdx  =>  dS^m/dR_i := sum_i dS^m/dR_i * dR_i/dx_j
	la.MatMul(o.G, 1, o.dSdR, o.dRdx)
	return
}

// CalcAtFaceIp calculates face data such as Sf and Fnvec
//  Input:
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//   ipf               -- local/natural coordinates of face
//   idxface           -- local index of face
//  Output:
//   Sf and Fnvec
func (o *Shape) CalcAtFaceIp(x [][]float64, ipf *Ipoint, idxface int) (err error) {

	// skip 1D elements
	if o.Gndim == 1 {
		return
	}

	// Sf and dSfdR
	o.FaceFunc(o.Sf, o.dSfdRf, ipf.R, ipf.S, ipf.T, true)

	// dxfdRf := sum_n x * dSfdRf   =>  dxf_i/dRf_j := sum_n xf^n_i * dSf^n/dRf_j
	for i := 0; i < len(x); i++ {
		for j := 0; j < o.Gndim-1; j++ {
			o.dxfdRf[i][j] = 0.0
			for k, n := range o.FaceLocalV[idxface] {
				o.dxfdRf[i][j] += x[i][n] * o.dSfdRf[k][j]
			}
		}
	}

	// face normal vector
	if o.Gndim == 2 {
		o.Fnvec[0] = o.dxfdRf[1][0]
		o.Fnvec[1] = -o.dxfdRf[0][0]
		return
	}
	o.Fnvec[0] = o.dxfdRf[1][0]*o.dxfdRf[2][1] - o.dxfdRf[2][0]*o.dxfdRf[1][1]
	o.Fnvec[1] = o.dxfdRf[2][0]*o.dxfdRf[0][1] - o.dxfdRf[0][0]*o.dxfdRf[2][1]
	o.Fnvec[2] = o.dxfdRf[0][0]*o.dxfdRf[1][1] - o.dxfdRf[1][0]*o.dxfdRf[0][1]
	return
}

// AxisymGetRadius returns the x0 == radius for axisymmetric computations
//  Note: must be called after CalcAtIp
func (o *Shape) AxisymGetRadius(x [][]float64) (radius float64) {
	for m := 0; m < o.Nverts; m++ {
		radius += o.S[m] * x[0][m]
	}
	return
}

// AxisymGetRadiusF (face) returns the x0 == radius for axisymmetric computations
//  Note: must be called after CalcAtFaceIp
func (o *Shape) AxisymGetRadiusF(x [][]float64, idxface int) (radius float64) {
	for m := 0; m < o.FaceNverts; m++ {
		radius += o.Sf[m] * x[0][o.FaceLocalV[idxface][m]]
	}
	return
}

// init_scratchpad initialise volume data (scratchpad)
func (o *Shape) init_scratchpad() {

	// volume data
	o.S = make([]float64, o.Nverts)
	o.dSdR = la.MatAlloc(o.Nverts, o.Gndim)
	o.dxdR = la.MatAlloc(o.Gndim, o.Gndim)
	o.dRdx = la.MatAlloc(o.Gndim, o.Gndim)
	o.G = la.MatAlloc(o.Nverts, o.Gndim)

	// face data
	if o.Gndim > 1 {
		o.Sf = make([]float64, o.FaceNverts)
		o.dSfdRf = la.MatAlloc(o.FaceNverts, o.Gndim-1)
		o.dxfdRf = la.MatAlloc(o.Gndim, o.Gndim-1)
		o.Fnvec = make([]float64, o.Gndim)
	}

	// lin data
	if o.Gndim == 1 {
		o.Jvec3d = make([]float64, 3)
		o.Gvec = make([]float64, o.Nverts)
	}
}
