// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be fndim in the LICENSE file.

package mporous

// scratchpad for TPM calculations
var TPM struct {

	// other porous media data
	Sg float64 // sl: saturation of gas
	Pc float64 // pc: capillary pressure: pg - pl
	P  float64 // p: averaged pressue of fluids in porous

	// n variables
	Ns float64 // Ns: volume fraction of solids
	Nf float64 // Nf: volume fraction of fluids: liquid + gas
	Nl float64 // Nl: volume fraction of liquid
	Ng float64 // Ng: volume fraction of gas

	// ρ (partial) variables
	Rhol float64 // ρl: partial density of liquid
	Rhog float64 // ρg: partial density of gas
	Rhos float64 // ρs: partial density of solids
	Rho  float64 // ρ: partial density of mixture

	// conductivity and retention models variables
	Klr float64 // relative liquid conductivity
	Kgr float64 // relative gas conductivity
	Ccb float64 // liquid retention model = dsl/dpc
	Cl  float64 // liquid compresssibility
	Cg  float64 // gas compressibility

	// liquid balance related coefficients
	Cpl float64 // Cpl coefficient
	Cpg float64 // Cpg coefficient
	Cvs float64 // Cvs coefficient

	// gas balance related coefficients
	Dpl float64 // Dpl coefficient
	Dpg float64 // Dpg coefficient
	Dvs float64 // Dvs coefficient

	// liquid retention model derivative
	Ccd float64 // dCc/dpc
}

// Derivatives
var D struct {

	// liquid conductivity
	klr_pl float64 // ∂klr/∂pl
	klr_pg float64 // ∂klr/∂pg

	// gas conductivity
	kgr_pl float64 // ∂kgr/∂pl
	kgr_pg float64 // ∂kgr/∂pg

	// mixture density
	rho_pl float64 // ∂ρ/∂pl
	rho_pg float64 // ∂ρ/∂pg

	// 1 liquid balance related coefficients
	Cpl_pl  float64 // ∂Cpl/∂pl
	Cpl_pg  float64 // ∂Cpl/∂pg
	Cpl_usM float64 // ∂Cpl/∂us multiplier

	// 2 liquid balance related coefficients
	Cpg_pl  float64 // ∂Cpg/∂pl
	Cpg_pg  float64 // ∂Cpg/∂pg
	Cpg_usM float64 // ∂Cpg/∂us multiplier

	// 3 liquid balance related coefficients
	Cvs_pl float64 // ∂Cvs/∂pl
	Cvs_pg float64 // ∂Cvs/∂pg

	// 1 gas balance related coefficients
	Dpl_pl  float64 // ∂Dpl/∂pl
	Dpl_pg  float64 // ∂Dpl/∂pg
	Dpl_usM float64 // ∂Dpl/∂us multiplier

	// 2 gas balance related coefficients
	Dpg_pl  float64 // ∂Dpg/∂pl
	Dpg_pg  float64 // ∂Dpg/∂pg
	Dpg_usM float64 // ∂Dpg/∂us multiplier

	// 3 gas balance related coefficients
	Dvs_pl float64 // ∂Dvs/∂pl
	Dvs_pg float64 // ∂Dvs/∂pg
}

// CalcLGS calculates TPM variables for models with liquid, gas and solid
func CalcLGS(divus float64, sta *State, mdl *Model, derivs bool) (err error) {

	// other porous media data
	TPM.Sg = 1.0 - sta.Sl
	TPM.Pc = sta.Pg - sta.Pl
	TPM.P = sta.Sl*sta.Pl + TPM.Sg*sta.Pg

	// n variables
	TPM.Ns = (1.0 - divus) * sta.Ns0
	TPM.Nf = 1.0 - TPM.Ns
	TPM.Nl = TPM.Nf * sta.Sl
	TPM.Ng = TPM.Nf * TPM.Sg

	// ρ (partial) variables
	TPM.Rhol = TPM.Nl * sta.RhoL
	TPM.Rhog = TPM.Ng * sta.RhoG
	TPM.Rhos = TPM.Ns * mdl.RhoS0
	TPM.Rho = TPM.Rhol + TPM.Rhog + TPM.Rhos

	// conductivity and retention models variables
	TPM.Klr = mdl.Cnd.Klr(sta.Sl)
	TPM.Kgr = mdl.Cnd.Kgr(TPM.Sg)
	TPM.Ccb, err = mdl.Ccb(sta)
	if err != nil {
		return
	}
	TPM.Cl = mdl.Cl
	TPM.Cg = mdl.Cg

	// liquid balance related coefficients
	TPM.Cpl = TPM.Nf * (sta.Sl*TPM.Cl - sta.RhoL*TPM.Ccb)
	TPM.Cpg = TPM.Nf * sta.RhoL * TPM.Ccb
	TPM.Cvs = sta.Sl * sta.RhoL

	// gas balance related coefficients
	TPM.Dpl = TPM.Nf * sta.RhoG * TPM.Ccb
	TPM.Dpg = TPM.Nf * (TPM.Sg*TPM.Cg - sta.RhoG*TPM.Ccb)
	TPM.Dvs = TPM.Sg * sta.RhoG

	// derivatives
	if derivs {

		// auxiliary variables
		ns0, nf := sta.Ns0, TPM.Nf
		sl, sg := sta.Sl, TPM.Sg
		ρL, ρG := sta.RhoL, sta.RhoG
		Cl, Cg := mdl.Cl, mdl.Cg
		Ccb, Ccd := TPM.Ccb, TPM.Ccd

		// liquid conductivity
		D.klr_pl = -mdl.Cnd.DklrDsl(sl) * Ccb
		D.klr_pg = +mdl.Cnd.DklrDsl(sl) * Ccb

		// gas conductivity
		D.kgr_pl = +mdl.Cnd.DkgrDsg(sg) * Ccb
		D.kgr_pg = -mdl.Cnd.DkgrDsg(sg) * Ccb

		// liquid retention model
		TPM.Ccd, err = mdl.Ccd(sta)
		if err != nil {
			return
		}

		// mixture density
		D.rho_pl = nf * (sl*Cl - ρL*Ccb + ρG*Ccb)
		D.rho_pg = nf * (sg*Cg - ρG*Ccb + ρL*Ccb)

		// 1 liquid balance related coefficients
		D.Cpl_pl = nf * (ρL*Ccd - 2.0*Ccb*Cl)
		D.Cpl_pg = nf * (Ccb*Cl - ρL*Ccd)
		D.Cpl_usM = (sl*Cl - ρL*Ccb) * ns0

		// 2 liquid balance related coefficients
		D.Cpg_pl = nf * (Cl*Ccb - ρL*Ccd)
		D.Cpg_pg = nf * ρL * Ccd
		D.Cpg_usM = ρL * Ccb * ns0

		// 3 liquid balance related coefficients
		D.Cvs_pl = sl*Cl - Ccb*ρL
		D.Cvs_pg = Ccb * ρL

		// 1 gas balance related coefficients
		D.Dpl_pl = -nf * ρG * Ccd
		D.Dpl_pg = nf * (ρG*Ccd + Cg*Ccb)
		D.Dpl_usM = ρG * Ccb * ns0

		// 2 gas balance related coefficients
		D.Dpg_pl = nf * (Ccb*Cg + ρG*Ccd)
		D.Dpg_pg = -nf * (ρG*Ccd + 2.0*Cg*Ccb)
		D.Dpg_usM = (sg*Cg - ρG*Ccb) * ns0

		// 3 gas balance related coefficients
		D.Dvs_pl = Ccb * ρG
		D.Dvs_pg = sg*Cg - Ccb*ρG
	}
	return
}
