// Copyright 2012 Dorival Pedroso & Raul Durand. All rights reserved.
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
	Cc  float64 // liquid retention model = ∂sl/∂pc
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

	// --- derivatives -------------------------------------

	// liquid conductivity
	DklrDpl float64 // ∂klr/∂pl
	DklrDpg float64 // ∂klr/∂pg

	// gas conductivity
	DkgrDpl float64 // ∂kgr/∂pl
	DkgrDpg float64 // ∂kgr/∂pg

	// liquid retention model
	DCcDpc float64 // ∂Cc/∂pc

	// mixture density
	DrhoDpl float64 // ∂ρ/∂pl
	DrhoDpg float64 // ∂ρ/∂pg

	// 1 liquid balance related coefficients
	DCplDpl  float64 // ∂Cpl/∂pl
	DCplDpg  float64 // ∂Cpl/∂pg
	DCplDusM float64 // ∂Cpl/∂us multiplier

	// 2 liquid balance related coefficients
	DCpgDpl  float64 // ∂Cpg/∂pl
	DCpgDpg  float64 // ∂Cpg/∂pg
	DCpgDusM float64 // ∂Cpg/∂us multiplier

	// 3 liquid balance related coefficients
	DCvsDpl float64 // ∂Cvs/∂pl
	DCvsDpg float64 // ∂Cvs/∂pg

	// 1 gas balance related coefficients
	DDplDpl  float64 // ∂Dpl/∂pl
	DDplDpg  float64 // ∂Dpl/∂pg
	DDplDusM float64 // ∂Dpl/∂us multiplier

	// 2 gas balance related coefficients
	DDpgDpl  float64 // ∂Dpg/∂pl
	DDpgDpg  float64 // ∂Dpg/∂pg
	DDpgDusM float64 // ∂Dpg/∂us multiplier

	// 3 gas balance related coefficients
	DDvsDpl float64 // ∂Dvs/∂pl
	DDvsDpg float64 // ∂Dvs/∂pg
}

// CalcLGS calculates TPM variables for models with liquid, gas and solid
func CalcLGS(divus float64, sta *StateLG, mdl *Model, derivs bool) (err error) {

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
	TPM.Rhos = TPM.Ns * mdl.RhoS
	TPM.Rho = TPM.Rhol + TPM.Rhog + TPM.Rhos

	// conductivity and retention models variables
	TPM.Klr = mdl.Cnd.Klr(sta.Sl)
	TPM.Kgr = mdl.Cnd.Kgr(TPM.Sg)
	TPM.Cc = mdl.Cc(TPM.Pc, sta.Sl, sta.Wet, sta.Dpc)
	TPM.Cl = mdl.Cl
	TPM.Cg = mdl.Cg

	// liquid balance related coefficients
	TPM.Cpl = TPM.Nf * (sta.Sl*TPM.Cl - sta.RhoL*TPM.Cc)
	TPM.Cpg = TPM.Nf * sta.RhoL * TPM.Cc
	TPM.Cvs = sta.Sl * sta.RhoL

	// gas balance related coefficients
	TPM.Dpl = TPM.Nf * sta.RhoG * TPM.Cc
	TPM.Dpg = TPM.Nf * (TPM.Sg*TPM.Cg - sta.RhoG*TPM.Cc)
	TPM.Dvs = TPM.Sg * sta.RhoG

	// derivatives
	if derivs {

		// auxiliary variables
		ns0, nf := sta.Ns0, TPM.Nf
		sl, sg := sta.Sl, TPM.Sg
		ρL, ρG := sta.RhoL, sta.RhoG
		Cl, Cg := mdl.Cl, mdl.Cg
		Cc := TPM.Cc

		// liquid conductivity
		TPM.DklrDpl = -mdl.Cnd.DklrDsl(sl) * Cc
		TPM.DklrDpg = +mdl.Cnd.DklrDsl(sl) * Cc

		// gas conductivity
		TPM.DkgrDpl = +mdl.Cnd.DkgrDsg(sg) * Cc
		TPM.DkgrDpg = -mdl.Cnd.DkgrDsg(sg) * Cc

		// liquid retention model
		Ccd := mdl.DCcDpc(TPM.Pc, sta.Sl, sta.Wet, sta.Dpc)
		TPM.DCcDpc = Ccd

		// mixture density
		TPM.DrhoDpl = nf * (sl*Cl - ρL*Cc + ρG*Cc)
		TPM.DrhoDpg = nf * (sg*Cg - ρG*Cc + ρL*Cc)

		// 1 liquid balance related coefficients
		TPM.DCplDpl = nf * (ρL*Ccd - 2.0*Cc*Cl)
		TPM.DCplDpg = nf * (Cc*Cl - ρL*Ccd)
		TPM.DCplDusM = (sl*Cl - ρL*Cc) * ns0

		// 2 liquid balance related coefficients
		TPM.DCpgDpl = nf * (Cl*Cc - ρL*Ccd)
		TPM.DCpgDpg = nf * ρL * Ccd
		TPM.DCpgDusM = ρL * Cc * ns0

		// 3 liquid balance related coefficients
		TPM.DCvsDpl = sl*Cl - Cc*ρL
		TPM.DCvsDpg = Cc * ρL

		// 1 gas balance related coefficients
		TPM.DDplDpl = -nf * ρG * Ccd
		TPM.DDplDpg = nf * (ρG*Ccd + Cg*Cc)
		TPM.DDplDusM = ρG * Cc * ns0

		// 2 gas balance related coefficients
		TPM.DDpgDpl = nf * (Cc*Cg + ρG*Ccd)
		TPM.DDpgDpg = -nf * (ρG*Ccd + 2.0*Cg*Cc)
		TPM.DDpgDusM = (sg*Cg - ρG*Cc) * ns0

		// 3 gas balance related coefficients
		TPM.DDvsDpl = Cc * ρG
		TPM.DDvsDpg = sg*Cg - Cc*ρG
	}
	return
}

// CalcLS calculates TPM variables for models with liquid and solid
//  With: pg = sg = ρG = 0
func CalcLS(divus float64, sta *StateLG, mdl *Model, derivs bool) (err error) {

	// other porous media data
	TPM.Pc = -sta.Pl
	TPM.P = sta.Sl * sta.Pl

	// n variables
	TPM.Ns = (1.0 - divus) * sta.Ns0
	TPM.Nf = 1.0 - TPM.Ns
	TPM.Nl = TPM.Nf * sta.Sl

	// ρ (partial) variables
	TPM.Rhol = TPM.Nl * sta.RhoL
	TPM.Rhos = TPM.Ns * mdl.RhoS
	TPM.Rho = TPM.Rhol + TPM.Rhos

	// conductivity and retention models variables
	TPM.Klr = mdl.Cnd.Klr(sta.Sl)
	TPM.Cc = mdl.Cc(TPM.Pc, sta.Sl, sta.Wet, sta.Dpc)
	TPM.Cl = mdl.Cl

	// liquid balance related coefficients
	TPM.Cpl = TPM.Nf * (sta.Sl*TPM.Cl - sta.RhoL*TPM.Cc)
	TPM.Cvs = sta.Sl * sta.RhoL

	// derivatives
	if derivs {

		// auxiliary variables
		ns0, nf := sta.Ns0, TPM.Nf
		sl := sta.Sl
		ρL := sta.RhoL
		Cl := mdl.Cl
		Cc := TPM.Cc

		// liquid conductivity
		TPM.DklrDpl = -mdl.Cnd.DklrDsl(sl) * Cc
		TPM.DklrDpg = +mdl.Cnd.DklrDsl(sl) * Cc

		// liquid retention model
		Ccd := mdl.DCcDpc(TPM.Pc, sta.Sl, sta.Wet, sta.Dpc)
		TPM.DCcDpc = Ccd

		// mixture density
		TPM.DrhoDpl = nf * (sl*Cl - ρL*Cc)

		// 1 liquid balance related coefficients
		TPM.DCplDpl = nf * (ρL*Ccd - 2.0*Cc*Cl)
		TPM.DCplDusM = (sl*Cl - ρL*Cc) * ns0

		// 2 liquid balance related coefficients
		TPM.DCpgDpl = nf * (Cl*Cc - ρL*Ccd)
		TPM.DCpgDusM = ρL * Cc * ns0

		// 3 liquid balance related coefficients
		TPM.DCvsDpl = sl*Cl - Cc*ρL
	}
	return
}

// CalcL calculates TPM variables for models with liquid only
//  With: pg = sg = ρG = 0 and divus = 0
func CalcL(sta *StateLG, mdl *Model, derivs bool) (err error) {

	// other porous media data
	TPM.Pc = -sta.Pl
	TPM.P = sta.Sl * sta.Pl

	// n variables
	TPM.Ns = sta.Ns0
	TPM.Nf = 1.0 - TPM.Ns
	TPM.Nl = TPM.Nf * sta.Sl

	// ρ (partial) variables
	TPM.Rhol = TPM.Nl * sta.RhoL
	TPM.Rhos = TPM.Ns * mdl.RhoS
	TPM.Rho = TPM.Rhol + TPM.Rhos

	// conductivity and retention models variables
	TPM.Klr = mdl.Cnd.Klr(sta.Sl)
	TPM.Cc = mdl.Cc(TPM.Pc, sta.Sl, sta.Wet, sta.Dpc)
	TPM.Cl = mdl.Cl

	// liquid balance related coefficients
	TPM.Cpl = TPM.Nf * (sta.Sl*TPM.Cl - sta.RhoL*TPM.Cc)

	// derivatives
	if derivs {

		// auxiliary variables
		nf := TPM.Nf
		sl := sta.Sl
		ρL := sta.RhoL
		Cl := mdl.Cl
		Cc := TPM.Cc

		// liquid conductivity
		TPM.DklrDpl = -mdl.Cnd.DklrDsl(sl) * Cc

		// liquid retention model
		Ccd := mdl.DCcDpc(TPM.Pc, sta.Sl, sta.Wet, sta.Dpc)
		TPM.DCcDpc = Ccd

		// mixture density
		TPM.DrhoDpl = nf * (sl*Cl - ρL*Cc)

		// 1 liquid balance related coefficients
		TPM.DCplDpl = nf * (ρL*Ccd - 2.0*Cc*Cl)
	}
	return
}
