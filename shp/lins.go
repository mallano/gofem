// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package shp

// shapes
var lin2, lin3, lin4, lin5 Shape

// register shapes
func init() {

	// lin2
	lin2.Type = "lin2"
	lin2.Func = Lin2
	lin2.BasicType = "lin2"
	lin2.Gndim = 1
	lin2.Nverts = 2
	lin2.VtkCode = VTK_LINE
	lin2.NatCoords = [][]float64{
		{-1, 1},
	}
	lin2.init_scratchpad()
	factory["lin2"] = &lin2
	ipsfactory["lin2_0"] = ips_lin_2
	ipsfactory["lin2_2"] = ips_lin_2

	// lin3
	lin3.Type = "lin3"
	lin3.Func = Lin3
	lin3.BasicType = "lin3"
	lin3.Gndim = 1
	lin3.Nverts = 3
	lin3.VtkCode = VTK_QUADRATIC_EDGE
	lin3.NatCoords = [][]float64{
		{-1, 1, 0},
	}
	lin3.init_scratchpad()
	factory["lin3"] = &lin3
	ipsfactory["lin3_0"] = ips_lin_2
	ipsfactory["lin3_2"] = ips_lin_2
	ipsfactory["lin3_3"] = ips_lin_3

	// lin4
	lin4.Type = "lin4"
	lin4.Func = Lin4
	lin4.BasicType = "lin4"
	lin4.Gndim = 1
	lin4.Nverts = 4
	lin4.VtkCode = VTK_POLY_LINE
	lin4.NatCoords = [][]float64{
		{-1, 1, -1.0 / 3.0, 1.0 / 3.0},
	}
	lin4.init_scratchpad()
	factory["lin4"] = &lin4
	ipsfactory["lin4_0"] = ips_lin_2
	ipsfactory["lin4_2"] = ips_lin_2
	ipsfactory["lin4_3"] = ips_lin_3
	ipsfactory["lin4_5"] = ips_lin_5

	// lin5
	lin5.Type = "lin5"
	lin5.Func = Lin2
	lin5.BasicType = "lin5"
	lin5.Gndim = 1
	lin5.Nverts = 2
	lin5.VtkCode = VTK_POLY_LINE
	lin5.NatCoords = [][]float64{
		{-1, 1, 0, -0.5, 0.5},
	}
	lin5.init_scratchpad()
	factory["lin5"] = &lin5
	ipsfactory["lin5_0"] = ips_lin_3
	ipsfactory["lin5_3"] = ips_lin_3
	ipsfactory["lin5_5"] = ips_lin_5
}

// Lin2 calculates the shape functions (S) and derivatives of shape functions (dSdR) of lin2
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Lin2(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*
	   -1     0    +1
	    0-----------1-->r
	*/
	S[0] = 0.5 * (1.0 - r)
	S[1] = 0.5 * (1.0 + r)

	if !derivs {
		return
	}

	dSdR[0][0] = -0.5
	dSdR[1][0] = 0.5
}

// Lin3 calculates the shape functions (S) and derivatives of shape functions (dSdR) of lin3
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Lin3(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*
	   -1     0    +1
	    0-----2-----1-->r
	*/
	S[0] = 0.5 * (r*r - r)
	S[1] = 0.5 * (r*r + r)
	S[2] = 1.0 - r*r

	if !derivs {
		return
	}

	dSdR[0][0] = r - 0.5
	dSdR[1][0] = r + 0.5
	dSdR[2][0] = -2.0 * r
}

// Lin4 calculates the shape functions (S) and derivatives of shape functions (dSdR) of lin4
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Lin4(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*
	   -1                  +1
	   @------@-----@------@  --> r
	   0      2     3      1
	*/
	S[0] = (-9.0*r*r*r + 9.0*r*r + r - 1.0) / 16.0
	S[1] = (9.0*r*r*r + 9.0*r*r - r - 1.0) / 16.0
	S[2] = (27.0*r*r*r - 9.0*r*r - 27.0*r + 9.0) / 16.0
	S[3] = (-27.0*r*r*r - 9.0*r*r + 27.0*r + 9.0) / 16.0

	if !derivs {
		return
	}

	dSdR[0][0] = 1.0 / 16.0 * (-27.0*r*r + 18.0*r + 1.0)
	dSdR[1][0] = 1.0 / 16.0 * (27.0*r*r + 18.0*r - 1.0)
	dSdR[2][0] = 1.0 / 16.0 * (81.0*r*r - 18.0*r - 27.0)
	dSdR[3][0] = 1.0 / 16.0 * (-81.0*r*r - 18.0*r + 27.0)
}

// Lin5 calculates the shape functions (S) and derivatives of shape functions (dSdR) of lin5
// elements at {r,s,t} natural coordinates. The derivatives are calculated only if derivs==true.
func Lin5(S []float64, dSdR [][]float64, r, s, t float64, derivs bool) {
	/*
	    @-----@-----@-----@-----@-> r
	    0     3     2     4     1
	    |           |           |
	   r=-1  -1/2   r=0  1/2   r=+1
	*/
	S[0] = (r - 1.0) * (1.0 - 2.0*r) * r * (-1.0 - 2.0*r) / 6.0
	S[1] = (1.0 - 2.0*r) * r * (-1.0 - 2.0*r) * (1.0 + r) / 6.0
	S[2] = (1.0 - r) * (1.0 - 2.0*r) * (-1.0 - 2.0*r) * (-1.0 - r)
	S[3] = 4.0 * (1.0 - r) * (1.0 - 2.0*r) * r * (-1.0 - r) / 3.0
	S[4] = 4.0 * (1.0 - r) * r * (-1.0 - 2.0*r) * (-1.0 - r) / 3.0

	if !derivs {
		return
	}

	dSdR[0][0] = -((1.0-2.0*r)*(r-1.0)*r)/3.0 - ((-2.0*r-1.0)*(r-1.0)*r)/3.0 + ((-2.0*r-1.0)*(1.0-2.0*r)*r)/6.0 + ((-2.0*r-1.0)*(1.0-2.0*r)*(r-1.0))/6.0
	dSdR[1][0] = -((1.0-2.0*r)*r*(r+1.0))/3.0 - ((-2.0*r-1.0)*r*(r+1.0))/3.0 + ((-2.0*r-1.0)*(1.0-2.0*r)*(r+1.0))/6.0 + ((-2.0*r-1.0)*(1.0-2.0*r)*r)/6.0
	dSdR[2][0] = -2.0*(1.0-2.0*r)*(-r-1.0)*(1.0-r) - 2.0*(-2.0*r-1.0)*(-r-1.0)*(1.0-r) - (-2.0*r-1.0)*(1.0-2.0*r)*(1.0-r) - (-2.0*r-1.0)*(1.0-2.0*r)*(-r-1.0)
	dSdR[3][0] = -(8.0*(-r-1.0)*(1.0-r)*r)/3.0 - (4.0*(1.0-2.0*r)*(1.0-r)*r)/3.0 - (4.0*(1.0-2.0*r)*(-r-1.0)*r)/3.0 + (4.0*(1.0-2.0*r)*(-r-1.0)*(1.0-r))/3.0
	dSdR[4][0] = -(8.0*(-r-1.0)*(1.0-r)*r)/3.0 - (4.0*(-2.0*r-1.0)*(1.0-r)*r)/3.0 - (4.0*(-2.0*r-1.0)*(-r-1.0)*r)/3.0 + (4.0*(-2.0*r-1.0)*(-r-1.0)*(1.0-r))/3.0
}
