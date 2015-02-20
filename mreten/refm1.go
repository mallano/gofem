// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package mreten

import (
	"math"

	"github.com/cpmech/gosl/chk"
	"github.com/cpmech/gosl/fun"
)

// RefM1 implements a nonlinear liquid retention model based on the concept of references [1,2,3]
//  References:
//   [1] Pedroso DM, Sheng D and Zhao, J (2009) The concept of reference curves for constitutive
//       modelling in soil mechanics, Computers and Geotechnics, 36(1-2), 149-165,
//       http://dx.doi.org/10.1016/j.compgeo.2008.01.009
//   [2] Pedroso DM and Williams DJ (2010) A novel approach for modelling soil-water
//       characteristic curves with hysteresis, Computers and Geotechnics, 37(3), 374-380,
//       http://dx.doi.org/10.1016/j.compgeo.2009.12.004
//   [3] Pedroso DM and Williams DJ (2011) Automatic Calibration of soil-water characteristic
//       curves using genetic algorithms. Computers and Geotechnics, 38(3), 330-340,
//       http://dx.doi.org/10.1016/j.compgeo.2010.12.004
type RefM1 struct {

	// parameters
	λd, λw, xrd, xrw float64 // parameters
	βd, βw, β1, β2   float64 // parameters
	y0, yr, α        float64 // initial y and parameters
	c1d, c2d, c3d    float64 // constants
	c1w, c2w, c3w    float64 // constants
	nowet            bool    // disregard wetting path
	model            int     // model number:  -1:Cc(pc)  0:default  1:Cc(pc,sl)
	A, pcr           float64 // Amplitude and pc_ref for non-default model

	// temporary variables (calculated by wetting or drying functions)
	Dw, λwb, yw      float64 // [temporary] wetting
	Dd, λdb, β2b, yd float64 // [temporary] drying
	D, λb            float64 // [temporary] wetting/drying
}

// add model to factory
func init() {
	allocators["ref-m1"] = func() Model { return new(RefM1) }
}

// Init initialises model
func (o *RefM1) Init(prms fun.Prms) (err error) {

	o.model, o.A, o.pcr = 0, 70.0, 0.02
	for _, p := range prms {
		switch p.N {
		case "lamd":
			o.λd = p.V
		case "lamw":
			o.λw = p.V
		case "xrd":
			o.xrd = p.V
		case "xrw":
			o.xrw = p.V
		case "yr":
			o.yr = p.V
		case "betd":
			o.βd = p.V
		case "betw":
			o.βw = p.V
		case "bet1":
			o.β1 = p.V
		case "bet2":
			o.β2 = p.V
		case "alp":
			o.α = p.V
		case "nowet":
			if p.V > 0.0 {
				o.nowet = true
			}
		case "A":
			o.A = p.V
		case "pcr":
			o.pcr = p.V
		default:
			return chk.Err("ref-m1: parameter named %q is incorrect\n", p.N)
		}
	}

	if o.nowet {
		o.λw, o.xrw, o.βw, o.β1 = o.λd, o.xrd, o.β2, o.βd
	}

	o.y0 = 1.0
	o.c1d = o.βd * o.λd
	o.c2d = math.Exp(o.βd * o.yr)
	o.c3d = math.Exp(o.βd*(o.y0+o.λd*o.xrd)) - o.c2d*math.Exp(o.c1d*o.xrd)
	o.c1w = -o.βw * o.λw
	o.c2w = math.Exp(-o.βw * o.y0)
	o.c3w = math.Exp(-o.βw*o.λw*o.xrw) - o.c2w*math.Exp(o.c1w*o.xrw)
	return
}

// GetPrms gets (an example) of parameters
func (o RefM1) GetPrms(example bool) fun.Prms {
	return []*fun.Prm{
		&fun.Prm{N: "lamd", V: 3},
		&fun.Prm{N: "lamw", V: 3},
		&fun.Prm{N: "xrd", V: 2.0},
		&fun.Prm{N: "xrw", V: 2.0},
		&fun.Prm{N: "yr", V: 0.005},
		&fun.Prm{N: "betd", V: 2},
		&fun.Prm{N: "betw", V: 2},
		&fun.Prm{N: "bet1", V: 2},
		&fun.Prm{N: "bet2", V: 2},
		&fun.Prm{N: "alp", V: 0.5},
		&fun.Prm{N: "nowet", V: 0, Inact: true},
	}
}

// SlMin returns sl_min
func (o RefM1) SlMin() float64 {
	return o.yr
}

// compute Cc(pc,sl) := dsl/dpc
func (o *RefM1) Cc(pc, sl float64, wet bool) (Ccval float64, err error) {
	if pc <= 0 {
		return 0, nil
	}
	if sl < o.yr {
		sl = o.yr
	}
	x := math.Log(1.0 + pc)
	if wet && !o.nowet {
		o.wetting(x, sl)
	} else {
		o.drying(x, sl)
	}
	Ccval = -o.λb / (1.0 + pc)
	return
}

// L computes L = ∂Cc/∂pc
func (o RefM1) L(pc, sl float64, wet bool) (float64, error) {
	if pc <= 0 {
		return 0, nil
	}
	if sl < o.yr {
		sl = o.yr
	}
	x := math.Log(1.0 + pc)
	var DλbDx float64
	if wet && !o.nowet {
		o.wetting(x, sl)
		DywDx := -o.λw - o.c1w*o.c2w*math.Exp(o.c1w*x)/(o.βw*(o.c3w+o.c2w*math.Exp(o.c1w*x)))
		DλbDx = o.β1 * o.λb * DywDx
	} else {
		o.drying(x, sl)
		DydDx := -o.λd + o.c1d*o.c2d*math.Exp(o.c1d*x)/(o.βd*(o.c3d+o.c2d*math.Exp(o.c1d*x)))
		DλbDx = -o.β2b * o.λb * DydDx
	}
	den := 1.0 + pc
	den2 := den * den
	return (o.λb - DλbDx) / den2, nil
}

// J computes J = ∂Cc/∂sl
func (o RefM1) J(pc, sl float64, wet bool) (float64, error) {
	if pc <= 0 {
		return 0, nil
	}
	if sl < o.yr {
		sl = o.yr
	}
	x := math.Log(1.0 + pc)
	var DλbDy float64
	if wet && !o.nowet {
		o.wetting(x, sl)
		DλbDy = (o.βw*(o.λwb-o.λw) - o.λwb*o.β1) * math.Exp(-o.β1*o.D)
	} else {
		o.drying(x, sl)
		Dβ2bDy := o.α * o.β2 * math.Pow(max(sl, 0.0), o.α-1.0)
		DλbDy = (o.βd*(o.λd-o.λdb) + o.λdb*(o.β2b-Dβ2bDy*o.D)) * math.Exp(-o.β2b*o.D)
	}
	return -DλbDy / (1.0 + pc), nil
}

// derivatives
func (o *RefM1) Derivs(pc, sl float64, wet bool) (L, Lx, J, Jx, Jy float64, err error) {
	if pc <= 0 {
		return
	}
	if sl < o.yr {
		sl = o.yr
	}
	x := math.Log(1.0 + pc)
	var DλbDx, DλbDy, D2λbDx2, D2λbDyDx, D2λbDy2 float64
	if wet && !o.nowet {
		o.wetting(x, sl)
		DywDx := -o.λw - o.c1w*o.c2w*math.Exp(o.c1w*x)/(o.βw*(o.c3w+o.c2w*math.Exp(o.c1w*x)))
		DλbDx = o.β1 * o.λb * DywDx
		DλbDy = (o.βw*(o.λwb-o.λw) - o.λwb*o.β1) * math.Exp(-o.β1*o.D)
		D2ywDx2 := -o.c1w * o.c1w * o.c2w * o.c3w * math.Exp(o.c1w*x) / (o.βw * math.Pow(o.c3w+o.c2w*math.Exp(o.c1w*x), 2.0))
		D2λbDx2 = o.β1 * (DλbDx*DywDx + o.λb*D2ywDx2)
		D2λbDyDx = o.β1 * DλbDy * DywDx
		DλwbDy := -o.λw * math.Exp(-o.βw*o.Dw) * o.βw
		D2λbDy2 = (o.βw*DλwbDy-DλwbDy*o.β1)*math.Exp(-o.β1*o.D) - (o.βw*(o.λwb-o.λw)-o.λwb*o.β1)*math.Exp(-o.β1*o.D)*o.β1
	} else {
		o.drying(x, sl)
		DydDx := -o.λd + o.c1d*o.c2d*math.Exp(o.c1d*x)/(o.βd*(o.c3d+o.c2d*math.Exp(o.c1d*x)))
		DλbDx = -o.β2b * o.λb * DydDx
		Dβ2bDy := o.α * o.β2 * math.Pow(max(sl, 0.0), o.α-1.0)
		DλbDy = (o.βd*(o.λd-o.λdb) + o.λdb*(o.β2b-Dβ2bDy*o.D)) * math.Exp(-o.β2b*o.D)
		D2ydDx2 := o.c1d * o.c1d * o.c2d * o.c3d * math.Exp(o.c1d*x) / (o.βd * math.Pow(o.c3d+o.c2d*math.Exp(o.c1d*x), 2.0))
		D2λbDx2 = -o.β2b * (DλbDx*DydDx + o.λb*D2ydDx2)
		D2λbDyDx = -(o.β2b*DλbDy + o.λb*Dβ2bDy) * DydDx
		DλdbDy := o.λd * math.Exp(-o.βd*o.Dd) * o.βd
		Dβ2bDy2 := o.α * o.β2 * math.Pow(max(sl, 0.0), o.α-2.0) * (o.α - 1.0)
		D2λbDy2 = (-o.βd*DλdbDy+DλdbDy*(o.β2b-Dβ2bDy*o.D)+o.λdb*(2.0*Dβ2bDy-Dβ2bDy2*o.D))*math.Exp(-o.β2b*o.D) + (o.βd*(o.λd-o.λdb)+o.λdb*(o.β2b-Dβ2bDy*o.D))*math.Exp(-o.β2b*o.D)*(-Dβ2bDy*o.D+o.β2b)
	}
	den := 1.0 + pc
	den2 := den * den
	L = (o.λb - DλbDx) / den2
	Lx = ((DλbDx-D2λbDx2)/den2 - 2.0*L) / den
	J = -DλbDy / den
	Jx = -(D2λbDyDx/den + J) / den
	Jy = -D2λbDy2 / den
	return
}

// auxiliary //////////////////////////////////////////////////////////////////////////////////////

func (o *RefM1) drying(x, y float64) {
	o.Dd = max(y-o.yr, 0.0)
	o.λdb = o.λd * (1.0 - math.Exp(-o.βd*o.Dd))
	o.yd = -o.λd*x + math.Log(o.c3d+o.c2d*math.Exp(o.c1d*x))/o.βd
	o.D = max(o.yd-y, 0.0)
	o.β2b = o.β2 * math.Pow(max(y, 0.0), o.α)
	o.λb = o.λdb * math.Exp(-o.β2b*o.D)
	return
}

func (o *RefM1) wetting(x, y float64) {
	o.Dw = max(o.y0-y, 0.0)
	o.λwb = o.λw * (1.0 - math.Exp(-o.βw*o.Dw))
	o.yw = -o.λw*x - math.Log(o.c3w+o.c2w*math.Exp(o.c1w*x))/o.βw
	o.D = max(y-o.yw, 0.0)
	o.λb = o.λwb * math.Exp(-o.β1*o.D)
	return
}
