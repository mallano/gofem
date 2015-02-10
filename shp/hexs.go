// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

// shapes
var hex8, hex20 Shape

// register shapes
func init() {

	// hex8
	hex8.Type = "hex8"
	hex8.Func = Hex8
	hex8.FaceFunc = Qua4
	hex8.BasicType = "hex8"
	hex8.FaceType = "qua4"
	hex8.Gndim = 3
	hex8.Nverts = 8
	hex8.VtkCode = VTK_HEXAHEDRON
	hex8.FaceNverts = 4
	hex8.FaceLocalV = [][]int{{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 3, 2, 1}, {4, 5, 6, 7}}
	hex8.NatCoords = [][]float64{
		{-1, 1, 1, -1, -1, 1, 1, -1},
		{-1, -1, 1, 1, -1, -1, 1, 1},
		{-1, -1, -1, -1, 1, 1, 1, 1},
	}
	hex8.init_scratchpad()
	factory["hex8"] = &hex8
	ipsfactory["hex8_0"] = ips_hex_8
	ipsfactory["hex8_8"] = ips_hex_8
	ipsfactory["hex8_14"] = ips_hex_14
	ipsfactory["hex8_27"] = ips_hex_27

	// hex20
	hex20.Type = "hex20"
	hex20.Func = Hex20
	hex20.FaceFunc = Qua8
	hex20.BasicType = "hex8"
	hex20.FaceType = "qua8"
	hex20.Gndim = 3
	hex20.Nverts = 20
	hex20.VtkCode = VTK_QUADRATIC_HEXAHEDRON
	hex20.FaceNverts = 8
	hex20.FaceLocalV = [][]int{{0, 4, 7, 3, 16, 15, 19, 11}, {1, 2, 6, 5, 9, 18, 13, 17}, {0, 1, 5, 4, 8, 17, 12, 16}, {2, 3, 7, 6, 10, 19, 14, 18}, {0, 3, 2, 1, 11, 10, 9, 8}, {4, 5, 6, 7, 12, 13, 14, 15}}
	hex20.NatCoords = [][]float64{
		{-1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1, -1},
		{-1, -1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1},
		{-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0},
	}
	hex20.init_scratchpad()
	factory["hex20"] = &hex20
	ipsfactory["hex20_0"] = ips_hex_27
	ipsfactory["hex20_8"] = ips_hex_8
	ipsfactory["hex20_14"] = ips_hex_14
	ipsfactory["hex20_27"] = ips_hex_27

}

// Hex8 calculates the shape functions (S) and derivatives of shape functions (dSdR) of hex8
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Hex8(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*
	             4________________7
	           ,'|              ,'|
	         ,'  |            ,'  |
	       ,'    |          ,'    |
	     ,'      |        ,'      |
	   5'===============6'        |
	   |         |      |         |
	   |         |      |         |
	   |         0_____ | ________3
	   |       ,'       |       ,'
	   |     ,'         |     ,'
	   |   ,'           |   ,'
	   | ,'             | ,'
	   1________________2'
	*/
	S[0] = (1.0 - r - s + r*s - t + s*t + r*t - r*s*t) / 8.0
	S[1] = (1.0 + r - s - r*s - t + s*t - r*t + r*s*t) / 8.0
	S[2] = (1.0 + r + s + r*s - t - s*t - r*t - r*s*t) / 8.0
	S[3] = (1.0 - r + s - r*s - t - s*t + r*t + r*s*t) / 8.0
	S[4] = (1.0 - r - s + r*s + t - s*t - r*t + r*s*t) / 8.0
	S[5] = (1.0 + r - s - r*s + t - s*t + r*t - r*s*t) / 8.0
	S[6] = (1.0 + r + s + r*s + t + s*t + r*t + r*s*t) / 8.0
	S[7] = (1.0 - r + s - r*s + t + s*t - r*t - r*s*t) / 8.0

	if !derivs {
		return
	}

	dSdR[0][0] = (-1.0 + s + t - s*t) / 8.0
	dSdR[0][1] = (-1.0 + r + t - r*t) / 8.0
	dSdR[0][2] = (-1.0 + r + s - r*s) / 8.0

	dSdR[1][0] = (+1.0 - s - t + s*t) / 8.0
	dSdR[1][1] = (-1.0 - r + t + r*t) / 8.0
	dSdR[1][2] = (-1.0 - r + s + r*s) / 8.0

	dSdR[2][0] = (+1.0 + s - t - s*t) / 8.0
	dSdR[2][1] = (+1.0 + r - t - r*t) / 8.0
	dSdR[2][2] = (-1.0 - r - s - r*s) / 8.0

	dSdR[3][0] = (-1.0 - s + t + s*t) / 8.0
	dSdR[3][1] = (+1.0 - r - t + r*t) / 8.0
	dSdR[3][2] = (-1.0 + r - s + r*s) / 8.0

	dSdR[4][0] = (-1.0 + s - t + s*t) / 8.0
	dSdR[4][1] = (-1.0 + r - t + r*t) / 8.0
	dSdR[4][2] = (+1.0 - r - s + r*s) / 8.0

	dSdR[5][0] = (+1.0 - s + t - s*t) / 8.0
	dSdR[5][1] = (-1.0 - r - t - r*t) / 8.0
	dSdR[5][2] = (+1.0 + r - s - r*s) / 8.0

	dSdR[6][0] = (+1.0 + s + t + s*t) / 8.0
	dSdR[6][1] = (+1.0 + r + t + r*t) / 8.0
	dSdR[6][2] = (+1.0 + r + s + r*s) / 8.0

	dSdR[7][0] = (-1.0 - s - t - s*t) / 8.0
	dSdR[7][1] = (+1.0 - r + t - r*t) / 8.0
	dSdR[7][2] = (+1.0 - r + s - r*s) / 8.0
}

// Hex20 calculates the shape functions (S) and derivatives of shape functions (dSdR) of hex20
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Hex20(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*
	              4_______15_______7
	            ,'|              ,'|
	         12'  |            ,'  |
	        ,'    16         ,14   |
	      ,'      |        ,'      19
	    5'=====13========6'        |
	    |         |      |         |
	    |         |      |         |
	    |         0_____ | _11_____3
	   17       ,'       |       ,'
	    |     8'        18     ,'
	    |   ,'           |   ,10
	    | ,'             | ,'
	    1_______9________2'
	*/
	rp1 := 1.0 + r
	rm1 := 1.0 - r
	sp1 := 1.0 + s
	sm1 := 1.0 - s
	tp1 := 1.0 + t
	tm1 := 1.0 - t

	S[0] = rm1 * sm1 * tm1 * (-r - s - t - 2) / 8.0
	S[1] = rp1 * sm1 * tm1 * (r - s - t - 2) / 8.0
	S[2] = rp1 * sp1 * tm1 * (r + s - t - 2) / 8.0
	S[3] = rm1 * sp1 * tm1 * (-r + s - t - 2) / 8.0
	S[4] = rm1 * sm1 * tp1 * (-r - s + t - 2) / 8.0
	S[5] = rp1 * sm1 * tp1 * (r - s + t - 2) / 8.0
	S[6] = rp1 * sp1 * tp1 * (r + s + t - 2) / 8.0
	S[7] = rm1 * sp1 * tp1 * (-r + s + t - 2) / 8.0
	S[8] = (1.0 - r*r) * sm1 * tm1 / 4.0
	S[9] = rp1 * (1.0 - s*s) * tm1 / 4.0
	S[10] = (1.0 - r*r) * sp1 * tm1 / 4.0
	S[11] = rm1 * (1.0 - s*s) * tm1 / 4.0
	S[12] = (1.0 - r*r) * sm1 * tp1 / 4.0
	S[13] = rp1 * (1.0 - s*s) * tp1 / 4.0
	S[14] = (1.0 - r*r) * sp1 * tp1 / 4.0
	S[15] = rm1 * (1.0 - s*s) * tp1 / 4.0
	S[16] = rm1 * sm1 * (1.0 - t*t) / 4.0
	S[17] = rp1 * sm1 * (1.0 - t*t) / 4.0
	S[18] = rp1 * sp1 * (1.0 - t*t) / 4.0
	S[19] = rm1 * sp1 * (1.0 - t*t) / 4.0

	if !derivs {
		return
	}

	dSdR[0][0] = -0.125*sm1*tm1*(-r-s-t-2.0) - 0.125*rm1*sm1*tm1
	dSdR[1][0] = 0.125*sm1*tm1*(r-s-t-2.0) + 0.125*rp1*sm1*tm1
	dSdR[2][0] = 0.125*sp1*tm1*(r+s-t-2.0) + 0.125*rp1*sp1*tm1
	dSdR[3][0] = -0.125*sp1*tm1*(-r+s-t-2.0) - 0.125*rm1*sp1*tm1
	dSdR[4][0] = -0.125*sm1*tp1*(-r-s+t-2.0) - 0.125*rm1*sm1*tp1
	dSdR[5][0] = 0.125*sm1*tp1*(r-s+t-2.0) + 0.125*rp1*sm1*tp1
	dSdR[6][0] = 0.125*sp1*tp1*(r+s+t-2.0) + 0.125*rp1*sp1*tp1
	dSdR[7][0] = -0.125*sp1*tp1*(-r+s+t-2.0) - 0.125*rm1*sp1*tp1
	dSdR[8][0] = -0.5 * r * sm1 * tm1
	dSdR[9][0] = 0.25 * (1.0 - s*s) * tm1
	dSdR[10][0] = -0.5 * r * sp1 * tm1
	dSdR[11][0] = -0.25 * (1.0 - s*s) * tm1
	dSdR[12][0] = -0.5 * r * sm1 * tp1
	dSdR[13][0] = 0.25 * (1.0 - s*s) * tp1
	dSdR[14][0] = -0.5 * r * sp1 * tp1
	dSdR[15][0] = -0.25 * (1.0 - s*s) * tp1
	dSdR[16][0] = -0.25 * sm1 * (1.0 - t*t)
	dSdR[17][0] = 0.25 * sm1 * (1.0 - t*t)
	dSdR[18][0] = 0.25 * sp1 * (1.0 - t*t)
	dSdR[19][0] = -0.25 * sp1 * (1.0 - t*t)

	dSdR[0][1] = -0.125*rm1*tm1*(-r-s-t-2.0) - 0.125*rm1*sm1*tm1
	dSdR[1][1] = -0.125*rp1*tm1*(r-s-t-2.0) - 0.125*rp1*sm1*tm1
	dSdR[2][1] = 0.125*rp1*tm1*(r+s-t-2.0) + 0.125*rp1*sp1*tm1
	dSdR[3][1] = 0.125*rm1*tm1*(-r+s-t-2.0) + 0.125*rm1*sp1*tm1
	dSdR[4][1] = -0.125*rm1*tp1*(-r-s+t-2.0) - 0.125*rm1*sm1*tp1
	dSdR[5][1] = -0.125*rp1*tp1*(r-s+t-2.0) - 0.125*rp1*sm1*tp1
	dSdR[6][1] = 0.125*rp1*tp1*(r+s+t-2.0) + 0.125*rp1*sp1*tp1
	dSdR[7][1] = 0.125*rm1*tp1*(-r+s+t-2.0) + 0.125*rm1*sp1*tp1
	dSdR[8][1] = -0.25 * (1.0 - r*r) * tm1
	dSdR[9][1] = -0.5 * s * rp1 * tm1
	dSdR[10][1] = 0.25 * (1.0 - r*r) * tm1
	dSdR[11][1] = -0.5 * s * rm1 * tm1
	dSdR[12][1] = -0.25 * (1.0 - r*r) * tp1
	dSdR[13][1] = -0.5 * s * rp1 * tp1
	dSdR[14][1] = 0.25 * (1.0 - r*r) * tp1
	dSdR[15][1] = -0.5 * s * rm1 * tp1
	dSdR[16][1] = -0.25 * rm1 * (1.0 - t*t)
	dSdR[17][1] = -0.25 * rp1 * (1.0 - t*t)
	dSdR[18][1] = 0.25 * rp1 * (1.0 - t*t)
	dSdR[19][1] = 0.25 * rm1 * (1.0 - t*t)

	dSdR[0][2] = -0.125*rm1*sm1*(-r-s-t-2.0) - 0.125*rm1*sm1*tm1
	dSdR[1][2] = -0.125*rp1*sm1*(r-s-t-2.0) - 0.125*rp1*sm1*tm1
	dSdR[2][2] = -0.125*rp1*sp1*(r+s-t-2.0) - 0.125*rp1*sp1*tm1
	dSdR[3][2] = -0.125*rm1*sp1*(-r+s-t-2.0) - 0.125*rm1*sp1*tm1
	dSdR[4][2] = 0.125*rm1*sm1*(-r-s+t-2.0) + 0.125*rm1*sm1*tp1
	dSdR[5][2] = 0.125*rp1*sm1*(r-s+t-2.0) + 0.125*rp1*sm1*tp1
	dSdR[6][2] = 0.125*rp1*sp1*(r+s+t-2.0) + 0.125*rp1*sp1*tp1
	dSdR[7][2] = 0.125*rm1*sp1*(-r+s+t-2.0) + 0.125*rm1*sp1*tp1
	dSdR[8][2] = -0.25 * (1.0 - r*r) * sm1
	dSdR[9][2] = -0.25 * rp1 * (1.0 - s*s)
	dSdR[10][2] = -0.25 * (1.0 - r*r) * sp1
	dSdR[11][2] = -0.25 * rm1 * (1.0 - s*s)
	dSdR[12][2] = 0.25 * (1.0 - r*r) * sm1
	dSdR[13][2] = 0.25 * rp1 * (1.0 - s*s)
	dSdR[14][2] = 0.25 * (1.0 - r*r) * sp1
	dSdR[15][2] = 0.25 * rm1 * (1.0 - s*s)
	dSdR[16][2] = -0.5 * t * rm1 * sm1
	dSdR[17][2] = -0.5 * t * rp1 * sm1
	dSdR[18][2] = -0.5 * t * rp1 * sp1
	dSdR[19][2] = -0.5 * t * rm1 * sp1
}
