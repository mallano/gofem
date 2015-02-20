// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// +build !appengine,!heroku

package shp

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/la"
)

// constants
const (
	INVMAP_TOL = 1.0e-10 // tolerance for inverse mapping function
	INVMAP_NIT = 25      // maximum number of iterations for inverse mapping
)

// InvMap computes the natural coordinates r, given the real coordinate y
//  Input:
//   y[ndim]           -- are the 2D/3D point coordinates
//   x[ndim][nverts+?] -- coordinates matrix of solid element
//  Output:
//   r[3] -- are the natural coordinates of given point
func (o *Shape) InvMap(r, y []float64, x [][]float64) (err error) {

	// check
	if o.Gndim == 1 {
		return chk.Err("Inverse mapping is not implemented in 1D\n")
	}

	var δRnorm float64
	e := make([]float64, o.Gndim)  // residual
	δr := make([]float64, o.Gndim) // corrector
	r[0], r[1], r[2] = 0, 0, 0     // first trial
	it := 0
	derivs := true
	for it = 0; it < INVMAP_NIT; it++ {

		// shape functions and derivatives
		o.Func(o.S, o.dSdR, r[0], r[1], r[2], derivs)

		// residual: e = y - x * S
		for i := 0; i < o.Gndim; i++ {
			e[i] = y[i]
			for j := 0; j < o.Nverts; j++ {
				e[i] -= x[i][j] * o.S[j]
			}
		}

		// Jmat == dxdR = x * dSdR;
		for i := 0; i < len(x); i++ {
			for j := 0; j < o.Gndim; j++ {
				o.dxdR[i][j] = 0.0
				for k := 0; k < o.Nverts; k++ {
					o.dxdR[i][j] += x[i][k] * o.dSdR[k][j] // dxdR := x * dSdR
				}
			}
		}

		// Jimat == dRdx = Jmat.inverse();
		o.J, err = la.MatInv(o.dRdx, o.dxdR, MINDET)
		if err != nil {
			return
		}

		// corrector: dR = Jimat * e
		for i := 0; i < o.Gndim; i++ {
			δr[i] = 0.0
			for j := 0; j < o.Gndim; j++ {
				δr[i] += o.dRdx[i][j] * e[j]
			}
		}

		// converged?
		δRnorm = 0.0
		for i := 0; i < o.Gndim; i++ {
			r[i] += δr[i]
			δRnorm += δr[i] * δr[i]
			// fix r outside range
			if r[i] < -1.0 || r[i] > 1.0 {
				if math.Abs(r[i]-(-1.0)) < INVMAP_TOL {
					r[i] = -1.0
				}
				if math.Abs(r[i]-1.0) < INVMAP_TOL {
					r[i] = 1.0
				}
			}

		}
		if math.Sqrt(δRnorm) < INVMAP_TOL {
			break
		}
	}

	// check
	if it == INVMAP_NIT {
		return
	}
	return
}

// GetNodesNatCoordsMat returns the matrix (ξ) with natural coordinates of nodes,
// augmented by one column which is filled with ones [nverts][ndim+1]
func (o *Shape) GetNodesNatCoordsMat() (ξ [][]float64) {
	ξ = la.MatAlloc(o.Nverts, o.Gndim+1)
	for i := 0; i < o.Nverts; i++ {
		for j := 0; j < o.Gndim; j++ {
			ξ[i][j] = o.NatCoords[j][i]
		}
		ξ[i][o.Gndim] = 1.0
	}
	return
}

// GetIpsNatCoordsMat returns the matrix (\hat{ξ}) with natural coordinates of interation
// points, augmented by one column which is filled with ones [nip][ndim+1]
func (o *Shape) GetIpsNatCoordsMat(ips []*Ipoint) (ξh [][]float64) {
	nip := len(ips)
	ξh = la.MatAlloc(nip, o.Gndim+1)
	for i := 0; i < nip; i++ {
		ξh[i][0] = ips[i].R
		ξh[i][1] = ips[i].S
		if o.Gndim == 3 {
			ξh[i][2] = ips[i].T
		}
		ξh[i][o.Gndim] = 1.0
	}
	return
}

// GetShapeMatAtIps returns a matrix formed by computing the shape functions
// at all integration points [nip][nverts]
func (o *Shape) GetShapeMatAtIps(ips []*Ipoint) (N [][]float64) {
	nip := len(ips)
	N = la.MatAlloc(nip, o.Nverts)
	derivs := false
	for i := 0; i < nip; i++ {
		r := ips[i].R
		s := ips[i].S
		t := ips[i].T
		o.Func(o.S, o.dSdR, r, s, t, derivs)
		for j := 0; j < o.Nverts; j++ {
			N[i][j] = o.S[j]
		}
	}
	return
}

// Extrapolator computes the extrapolation matrix for this Shape with a combination of integration points 'ips'
//  Note: E[nverts][nip] must be pre-allocated
func (o *Shape) Extrapolator(E [][]float64, ips []*Ipoint) (err error) {
	la.MatFill(E, 0)
	nip := len(ips)
	N := o.GetShapeMatAtIps(ips)
	if nip < o.Nverts {
		ξ := o.GetNodesNatCoordsMat()
		ξh := o.GetIpsNatCoordsMat(ips)
		ξhi := la.MatAlloc(o.Gndim+1, nip)
		Ni := la.MatAlloc(o.Nverts, nip)
		err = la.MatInvG(Ni, N, 1e-10)
		if err != nil {
			return
		}
		err = la.MatInvG(ξhi, ξh, 1e-10)
		if err != nil {
			return
		}
		ξhξhI := la.MatAlloc(nip, nip) // ξh * inv(ξh)
		for k := 0; k < o.Gndim+1; k++ {
			for j := 0; j < nip; j++ {
				for i := 0; i < nip; i++ {
					ξhξhI[i][j] += ξh[i][k] * ξhi[k][j]
				}
				for i := 0; i < o.Nverts; i++ {
					E[i][j] += ξ[i][k] * ξhi[k][j] // ξ * inv(ξh)
				}
			}
		}
		for i := 0; i < o.Nverts; i++ {
			for j := 0; j < nip; j++ {
				for k := 0; k < nip; k++ {
					I_kj := 0.0
					if j == k {
						I_kj = 1.0
					}
					E[i][j] += Ni[i][k] * (I_kj - ξhξhI[k][j])
				}
			}
		}
	} else {
		err = la.MatInvG(E, N, 1e-10)
		if err != nil {
			return
		}
	}
	return
}
