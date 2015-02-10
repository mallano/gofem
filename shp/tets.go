// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

// shapes
var tet4, tet10 Shape

// register shapes
func init() {

	// tet4
	tet4.Type = "tet4"
	tet4.Func = Tet4
	tet4.FaceFunc = Tri3
	tet4.BasicType = "tet4"
	tet4.FaceType = "tri3"
	tet4.Gndim = 3
	tet4.Nverts = 4
	tet4.VtkCode = VTK_TETRA
	tet4.FaceNverts = 3
	tet4.FaceLocalV = [][]int{{0, 3, 2}, {0, 1, 3}, {0, 2, 1}, {1, 2, 3}}
	tet4.NatCoords = [][]float64{
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	}
	tet4.init_scratchpad()
	factory["tet4"] = &tet4
	ipsfactory["tet4_0"] = ips_tet_1
	ipsfactory["tet4_1"] = ips_tet_1
	ipsfactory["tet4_4"] = ips_tet_4
	ipsfactory["tet4_5"] = ips_tet_5

	// tet10
	tet10.Type = "tet10"
	tet10.Func = Tet10
	tet10.FaceFunc = Tri6
	tet10.BasicType = "tet4"
	tet10.FaceType = "tri6"
	tet10.Gndim = 3
	tet10.Nverts = 10
	tet10.VtkCode = VTK_QUADRATIC_TETRA
	tet10.FaceNverts = 6
	tet10.FaceLocalV = [][]int{{0, 3, 2, 7, 9, 6}, {0, 1, 3, 4, 8, 7}, {0, 2, 1, 6, 5, 4}, {1, 2, 3, 5, 9, 8}}
	tet10.NatCoords = [][]float64{
		{0, 1, 0, 0, 0.5, 0.5, 0, 0, 0.5, 0},
		{0, 0, 1, 0, 0, 0.5, 0.5, 0, 0, 0.5},
		{0, 0, 0, 1, 0, 0, 0, 0.5, 0.5, 0.5},
	}
	tet10.init_scratchpad()
	factory["tet10"] = &tet10
	ipsfactory["tet10_0"] = ips_tet_4
	ipsfactory["tet10_4"] = ips_tet_4
	ipsfactory["tet10_5"] = ips_tet_5
}

// Tet4 calculates the shape functions (S) and derivatives of shape functions (dSdR) of tet4
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Tet4(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*                    t
	              |
	              3
	             /|`.
	             ||  `,
	            / |    ',
	            | |      \
	           /  |       `.
	           |  |         `,
	          /   |           `,
	          |   |             \
	         /    |              `.
	         |    |                ',
	        /     |                  \
	        |     0.,,_               `.
	       |     /     ``'-.,,__        `.
	       |    /              ``''-.,,_  ',
	      |    /                        `` 2 ,,s
	      |  ,'                       ,.-``
	     |  ,                    _,-'`
	     ' /                 ,.'`
	    | /             _.-``
	    '/          ,-'`
	   |/      ,.-``
	   /  _,-``
	  1 '`
	 /
	r
	*/
	S[0] = 1.0 - r - s - t
	S[1] = r
	S[2] = s
	S[3] = t

	if !derivs {
		return
	}

	dSdR[0][0] = -1.0
	dSdR[1][0] = 1.0
	dSdR[2][0] = 0.0
	dSdR[3][0] = 0.0

	dSdR[0][1] = -1.0
	dSdR[1][1] = 0.0
	dSdR[2][1] = 1.0
	dSdR[3][1] = 0.0

	dSdR[0][2] = -1.0
	dSdR[1][2] = 0.0
	dSdR[2][2] = 0.0
	dSdR[3][2] = 1.0
}

// Tet10 calculates the shape functions (S) and derivatives of shape functions (dSdR) of tet10
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Tet10(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*                    t
	              |
	              3
	             /|`.
	             ||  `,
	            / |    ',
	            | |      \
	           /  |       `.
	           |  |         `,
	          /   7            9
	          |   |             \
	         /    |              `.
	         |    |                ',
	        8     |                  \
	        |     0 ,,_               `.
	       |     /     ``'-., 6         `.
	       |    /               `''-.,,_  ',
	      |    /                        ``'2 ,,s
	      |   '                       ,.-``
	     |   4                   _,-'`
	     ' /                 ,.'`
	    | /             _ 5 `
	    '/          ,-'`
	   |/      ,.-``
	   /  _,-``
	  1 '`
	 /
	r
	*/
	u := 1.0 - r - s - t
	S[0] = u * (2.0*u - 1.0)
	S[1] = r * (2.0*r - 1.0)
	S[2] = s * (2.0*s - 1.0)
	S[3] = t * (2.0*t - 1.0)
	S[4] = 4.0 * u * r
	S[5] = 4.0 * r * s
	S[6] = 4.0 * s * u
	S[7] = 4.0 * u * t
	S[8] = 4.0 * r * t
	S[9] = 4.0 * s * t

	if !derivs {
		return
	}

	dSdR[0][0] = 4.0*(r+s+t) - 3.0
	dSdR[1][0] = 4.0*r - 1.0
	dSdR[2][0] = 0.0
	dSdR[3][0] = 0.0
	dSdR[4][0] = 4.0 - 8.0*r - 4.0*s - 4.0*t
	dSdR[5][0] = 4.0 * s
	dSdR[6][0] = -4.0 * s
	dSdR[7][0] = -4.0 * t
	dSdR[8][0] = 4.0 * t
	dSdR[9][0] = 0.0

	dSdR[0][1] = 4.0*(r+s+t) - 3.0
	dSdR[1][1] = 0.0
	dSdR[2][1] = 4.0*s - 1.0
	dSdR[3][1] = 0.0
	dSdR[4][1] = -4.0 * r
	dSdR[5][1] = 4.0 * r
	dSdR[6][1] = 4.0 - 4.0*r - 8.0*s - 4.0*t
	dSdR[7][1] = -4.0 * t
	dSdR[8][1] = 0.0
	dSdR[9][1] = 4.0 * t

	dSdR[0][2] = 4.0*(r+s+t) - 3.0
	dSdR[1][2] = 0.0
	dSdR[2][2] = 0.0
	dSdR[3][2] = 4.0*t - 1.0
	dSdR[4][2] = -4.0 * r
	dSdR[5][2] = 0.0
	dSdR[6][2] = -4.0 * s
	dSdR[7][2] = 4.0 - 4.0*r - 4.0*s - 8.0*t
	dSdR[8][2] = 4.0 * r
	dSdR[9][2] = 4.0 * s
}
