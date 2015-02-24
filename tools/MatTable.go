// Copyright 2012 Dorival de Moraes Pedroso. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bytes"
	"flag"
	"strings"

	"github.com/cpmech/gofem/inp"
	"github.com/cpmech/gosl/io"
)

func main() {

	defer func() {
		if err := recover(); err != nil {
			io.PfRed("Some error has happened: %v\n", err)
		}
	}()

	// input data
	matfn := "materials.mat"

	skip := "gref nowet α"

	// parse flags
	flag.Parse()
	if len(flag.Args()) > 0 {
		matfn = flag.Arg(0)
	}
	if len(flag.Args()) > 1 {
		skip = flag.Arg(1)
	}

	// skip parameters
	skipp := make(map[string]bool)
	for _, key := range io.SplitKeys(skip) {
		skipp[key] = true
	}

	// file key
	fnk := io.FnKey(matfn)

	// print input data
	io.Pforan("Input data\n")
	io.Pfblue2("  matfn = %v\n", matfn)
	io.Pfblue2("  skip  = %v\n", skip)

	// Read
	mdb := inp.ReadMat("", matfn)

	// Get max number of parameters
	nmaxprms := 0
	for _, mdat := range mdb.Materials {
		n := len(mdat.Prms)
		if n > nmaxprms {
			nmaxprms = n
		}
	}

	// header
	b := new(bytes.Buffer)
	io.Ff(b, "\\documentclass[12pt,a4paper]{article}\n")
	io.Ff(b, "\\usepackage[margin=2.0cm,footskip=0.5cm]{geometry}\n")
	io.Ff(b, "\\usepackage[labelfont=bf,tableposition=top,aboveskip=4pt]{caption}\n")
	io.Ff(b, "\\usepackage{tabularx}\n")
	io.Ff(b, "\\usepackage{booktabs}\n")
	io.Ff(b, "\n")
	io.Ff(b, "\\title{Materials Table}\n")
	io.Ff(b, "\\author{GoFem MatTable tool}\n")
	io.Ff(b, "\n")
	io.Ff(b, "\\begin{document}\n")
	io.Ff(b, "\\maketitle\n")

	// table with parameters
	io.Ff(b, "\n")
	io.Ff(b, "\\begin{table} \\centering\n")
	io.Ff(b, "\\caption{Parameters from %s}\n", matfn)
	io.Ff(b, "\\begin{tabularx}{\\linewidth}[c]{l %s} \\toprule\n", strings.Repeat("c", nmaxprms))

	for i, mdat := range mdb.Materials {
		necols := nmaxprms - len(mdat.Prms) // number of empty cols

		// mat name
		io.Ff(b, "  %-20s", mdat.Name)

		// parameters names
		for _, param := range mdat.Prms {
			io.Ff(b, " &%12s", ToTex(param.N))
		}

		io.Ff(b, " %s", strings.Repeat(" &", necols))
		io.Ff(b, " \\\\\n")

		// values
		io.Ff(b, "  %-20s", "")
		for _, param := range mdat.Prms {
			io.Ff(b, " &%12s", NumFormat(param.V))
		}

		io.Ff(b, " %s", strings.Repeat(" &", necols))
		io.Ff(b, " \\\\\n")

		// units
		io.Ff(b, "  %-20s", "")
		for _, param := range mdat.Prms {
			if param.U == "" {
				io.Ff(b, " &%12v", "-")
			} else {
				io.Ff(b, " &%12v", UnitFormat(param.U))
			}
		}

		io.Ff(b, " %s", strings.Repeat(" &", necols))
		io.Ff(b, " \\\\\n")

		if i < len(mdb.Materials)-1 {
			io.Ff(b, "  \\\\\n")
		}
	}

	// footer
	io.Ff(b, "  \\bottomrule\n\\end{tabularx}\n\\label{tab:prms}\n\\end{table}\n")
	io.Ff(b, "\\end{document}\n")
	io.WriteFileV(fnk+".tex", b)

}

func NumFormat(num float64) (str string) {
	if num == 0 {
		return io.Sf("%g", num)
	}
	if num >= 10000 {
		mant, exp := mantexp(num)
		if mant == 1 {
			return io.Sf("$10^{%g}$", exp)
		} else {
			return io.Sf("$%g\\times10^{%g}$", mant, exp)

		}
	}

	return io.Sf("%g", num)
}

func mantexp(num float64) (mant, exp float64) {
	mant = num
	exp = 0
	for mant >= 10 {
		exp += 1
		mant /= 10
	}
	return
}

func UnitFormat(unit string) (tex string) {
	if strings.Contains(unit, "m3") {
		return io.Sf("$\\mathrm{%s}$", strings.Replace(unit, "m3", "m^3", -1))
	}
	if strings.Contains(unit, "m2") {
		return io.Sf("$\\mathrm{%s}$", strings.Replace(unit, "m2", "m^2", -1))
	}
	return unit
}

func ToTex(prm string) (tex string) {
	switch prm {
	case "ρS":
		tex = "\\rho^s"
	case "ν", "nu":
		tex = "\\nu"
	case "μ", "mu":
		tex = "\\mu"
	case "φ":
		tex = "\\phi"
	case "nf":
		tex = "n_f"
	case "kl":
		tex = "k_\\ell"
	case "kg":
		tex = "k_g"
	case "βl":
		tex = "\\beta_\\ell"
	case "βg":
		tex = "\\beta_g"
	case "ρ", "rho":
		tex = "\\rho"
	case "ρL":
		tex = "\\rho^\\ell"
	case "ρG":
		tex = "\\rho^g"
	case "κ", "kappa", "kap":
		tex = "\\kappa"
	case "ks":
		tex = "k_s"
	case "k1":
		tex = "k_1"
	case "k2":
		tex = "k_2"
	case "kh":
		tex = "k_h"
	case "Kl":
		tex = "K_\\ell"
	case "RΘg":
		tex = "R_{\\Theta g}"
	case "λd":
		tex = "\\lambda_d"
	case "λw":
		tex = "\\lambda_w"
	case "yr":
		tex = "y_r"
	case "xrd":
		tex = "x_{rd}"
	case "xrw":
		tex = "x_{rw}"
	case "βd":
		tex = "\\beta_d"
	case "βw":
		tex = "\\beta_w"
	case "β1":
		tex = "\\beta_1"
	case "β2":
		tex = "\\beta_2"
	case "α", "alp", "alpha":
		tex = "\\alpha"
	case "αl":
		tex = "\\alpha_\\ell"
	case "λ0l":
		tex = "\\lambda_{0\\ell}"
	case "λ1l":
		tex = "\\lambda_{1\\ell}"
	case "λ0g":
		tex = "\\lambda_{0g}"
	case "λ1g":
		tex = "\\lambda_{1g}"
	case "tauy0":
		tex = "\\tau_{y0}"
	default:
		tex = prm
	}
	return "$" + tex + "$"
}
