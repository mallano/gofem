// Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be fndim in the LICENSE file.

package fem

import (
	"math"

	"github.com/cpmech/gosl/la"
	"github.com/cpmech/gosl/tsr"
)

const SQ2 = math.Sqrt2

func IpAddToKt(Kt [][]float64, nne, ndim int, coef float64, G, D [][]float64) {
	if ndim == 3 {
		for m := 0; m < nne; m++ {
			for n := 0; n < nne; n++ {
				Kt[0+m*3][0+n*3] += coef * (G[m][2]*G[n][2]*D[5][5] + G[m][2]*G[n][1]*D[5][3] + SQ2*G[m][2]*G[n][0]*D[5][0] + G[m][1]*G[n][2]*D[3][5] + G[m][1]*G[n][1]*D[3][3] + SQ2*G[m][1]*G[n][0]*D[3][0] + SQ2*G[m][0]*G[n][2]*D[0][5] + SQ2*G[m][0]*G[n][1]*D[0][3] + 2.0*G[m][0]*G[n][0]*D[0][0]) / 2.0
				Kt[0+m*3][1+n*3] += coef * (G[m][2]*G[n][2]*D[5][4] + G[m][2]*G[n][0]*D[5][3] + SQ2*G[m][2]*G[n][1]*D[5][1] + G[m][1]*G[n][2]*D[3][4] + G[m][1]*G[n][0]*D[3][3] + SQ2*G[m][1]*G[n][1]*D[3][1] + SQ2*G[m][0]*G[n][2]*D[0][4] + SQ2*G[m][0]*G[n][0]*D[0][3] + 2.0*G[m][0]*G[n][1]*D[0][1]) / 2.0
				Kt[0+m*3][2+n*3] += coef * (G[m][2]*G[n][0]*D[5][5] + G[m][2]*G[n][1]*D[5][4] + SQ2*G[m][2]*G[n][2]*D[5][2] + G[m][1]*G[n][0]*D[3][5] + G[m][1]*G[n][1]*D[3][4] + SQ2*G[m][1]*G[n][2]*D[3][2] + SQ2*G[m][0]*G[n][0]*D[0][5] + SQ2*G[m][0]*G[n][1]*D[0][4] + 2.0*G[m][0]*G[n][2]*D[0][2]) / 2.0
				Kt[1+m*3][0+n*3] += coef * (G[m][2]*G[n][2]*D[4][5] + G[m][2]*G[n][1]*D[4][3] + SQ2*G[m][2]*G[n][0]*D[4][0] + G[m][0]*G[n][2]*D[3][5] + G[m][0]*G[n][1]*D[3][3] + SQ2*G[m][0]*G[n][0]*D[3][0] + SQ2*G[m][1]*G[n][2]*D[1][5] + SQ2*G[m][1]*G[n][1]*D[1][3] + 2.0*G[m][1]*G[n][0]*D[1][0]) / 2.0
				Kt[1+m*3][1+n*3] += coef * (G[m][2]*G[n][2]*D[4][4] + G[m][2]*G[n][0]*D[4][3] + SQ2*G[m][2]*G[n][1]*D[4][1] + G[m][0]*G[n][2]*D[3][4] + G[m][0]*G[n][0]*D[3][3] + SQ2*G[m][0]*G[n][1]*D[3][1] + SQ2*G[m][1]*G[n][2]*D[1][4] + SQ2*G[m][1]*G[n][0]*D[1][3] + 2.0*G[m][1]*G[n][1]*D[1][1]) / 2.0
				Kt[1+m*3][2+n*3] += coef * (G[m][2]*G[n][0]*D[4][5] + G[m][2]*G[n][1]*D[4][4] + SQ2*G[m][2]*G[n][2]*D[4][2] + G[m][0]*G[n][0]*D[3][5] + G[m][0]*G[n][1]*D[3][4] + SQ2*G[m][0]*G[n][2]*D[3][2] + SQ2*G[m][1]*G[n][0]*D[1][5] + SQ2*G[m][1]*G[n][1]*D[1][4] + 2.0*G[m][1]*G[n][2]*D[1][2]) / 2.0
				Kt[2+m*3][0+n*3] += coef * (G[m][0]*G[n][2]*D[5][5] + G[m][0]*G[n][1]*D[5][3] + SQ2*G[m][0]*G[n][0]*D[5][0] + G[m][1]*G[n][2]*D[4][5] + G[m][1]*G[n][1]*D[4][3] + SQ2*G[m][1]*G[n][0]*D[4][0] + SQ2*G[m][2]*G[n][2]*D[2][5] + SQ2*G[m][2]*G[n][1]*D[2][3] + 2.0*G[m][2]*G[n][0]*D[2][0]) / 2.0
				Kt[2+m*3][1+n*3] += coef * (G[m][0]*G[n][2]*D[5][4] + G[m][0]*G[n][0]*D[5][3] + SQ2*G[m][0]*G[n][1]*D[5][1] + G[m][1]*G[n][2]*D[4][4] + G[m][1]*G[n][0]*D[4][3] + SQ2*G[m][1]*G[n][1]*D[4][1] + SQ2*G[m][2]*G[n][2]*D[2][4] + SQ2*G[m][2]*G[n][0]*D[2][3] + 2.0*G[m][2]*G[n][1]*D[2][1]) / 2.0
				Kt[2+m*3][2+n*3] += coef * (G[m][0]*G[n][0]*D[5][5] + G[m][0]*G[n][1]*D[5][4] + SQ2*G[m][0]*G[n][2]*D[5][2] + G[m][1]*G[n][0]*D[4][5] + G[m][1]*G[n][1]*D[4][4] + SQ2*G[m][1]*G[n][2]*D[4][2] + SQ2*G[m][2]*G[n][0]*D[2][5] + SQ2*G[m][2]*G[n][1]*D[2][4] + 2.0*G[m][2]*G[n][2]*D[2][2]) / 2.0
			}
		}
	} else {
		for m := 0; m < nne; m++ {
			for n := 0; n < nne; n++ {
				Kt[0+m*2][0+n*2] += coef * (G[m][1]*G[n][1]*D[3][3] + SQ2*G[m][1]*G[n][0]*D[3][0] + SQ2*G[m][0]*G[n][1]*D[0][3] + 2.0*G[m][0]*G[n][0]*D[0][0]) / 2.0
				Kt[0+m*2][1+n*2] += coef * (G[m][1]*G[n][0]*D[3][3] + SQ2*G[m][1]*G[n][1]*D[3][1] + SQ2*G[m][0]*G[n][0]*D[0][3] + 2.0*G[m][0]*G[n][1]*D[0][1]) / 2.0
				Kt[1+m*2][0+n*2] += coef * (G[m][0]*G[n][1]*D[3][3] + SQ2*G[m][0]*G[n][0]*D[3][0] + SQ2*G[m][1]*G[n][1]*D[1][3] + 2.0*G[m][1]*G[n][0]*D[1][0]) / 2.0
				Kt[1+m*2][1+n*2] += coef * (G[m][0]*G[n][0]*D[3][3] + SQ2*G[m][0]*G[n][1]*D[3][1] + SQ2*G[m][1]*G[n][0]*D[1][3] + 2.0*G[m][1]*G[n][1]*D[1][1]) / 2.0
			}
		}
	}
}

func IpStrains(εs []float64, nne, ndim int, u []float64, Umap []int, G [][]float64) {
	var r, c int
	var εsij float64
	for i := 0; i < ndim; i++ {
		for j := i; j < ndim; j++ { // note: j := i => only diagonal and above
			εsij = 0
			for m := 0; m < nne; m++ {
				r, c = i+m*ndim, j+m*ndim
				εsij += (u[Umap[r]]*G[m][j] + u[Umap[c]]*G[m][i]) / 2.0
			}
			if i != j {
				εsij *= SQ2
			}
			εs[tsr.T2MI[i][j]] = εsij
		}
	}
}

func IpStrainsAndInc(εs, Δεs []float64, nne, ndim int, u, Δu []float64, Umap []int, G [][]float64) {
	var r, c int
	var εsij, Δεsij float64
	for i := 0; i < ndim; i++ {
		for j := i; j < ndim; j++ { // note: j := i => only diagonal and above
			εsij, Δεsij = 0, 0
			for m := 0; m < nne; m++ {
				r, c = i+m*ndim, j+m*ndim
				εsij += (u[Umap[r]]*G[m][j] + u[Umap[c]]*G[m][i]) / 2.0
				Δεsij += (Δu[Umap[r]]*G[m][j] + Δu[Umap[c]]*G[m][i]) / 2.0
			}
			if i != j {
				εsij *= SQ2
				Δεsij *= SQ2
			}
			εs[tsr.T2MI[i][j]] = εsij
			Δεs[tsr.T2MI[i][j]] = Δεsij
		}
	}
}

// DerivSig returns the derivative of σ (Mandel) with respect to displacement at nodes
//  Note: DσDun = ∂σ/∂un  [nσ][ndim]
func DerivSig(DσDun [][]float64, n, ndim int, G, D [][]float64) {
	if ndim == 3 {
		DσDun[0][0] = (G[n][2]*D[0][5]*SQ2 + G[n][1]*D[0][3]*SQ2 + 2.0*G[n][0]*D[0][0]) / 2.0
		DσDun[0][1] = (G[n][2]*D[0][4]*SQ2 + G[n][0]*D[0][3]*SQ2 + 2.0*G[n][1]*D[0][1]) / 2.0
		DσDun[0][2] = (G[n][0]*D[0][5]*SQ2 + G[n][1]*D[0][4]*SQ2 + 2.0*G[n][2]*D[0][2]) / 2.0
		DσDun[1][0] = (G[n][2]*D[1][5]*SQ2 + G[n][1]*D[1][3]*SQ2 + 2.0*G[n][0]*D[1][0]) / 2.0
		DσDun[1][1] = (G[n][2]*D[1][4]*SQ2 + G[n][0]*D[1][3]*SQ2 + 2.0*G[n][1]*D[1][1]) / 2.0
		DσDun[1][2] = (G[n][0]*D[1][5]*SQ2 + G[n][1]*D[1][4]*SQ2 + 2.0*G[n][2]*D[1][2]) / 2.0
		DσDun[2][0] = (G[n][2]*D[2][5]*SQ2 + G[n][1]*D[2][3]*SQ2 + 2.0*G[n][0]*D[2][0]) / 2.0
		DσDun[2][1] = (G[n][2]*D[2][4]*SQ2 + G[n][0]*D[2][3]*SQ2 + 2.0*G[n][1]*D[2][1]) / 2.0
		DσDun[2][2] = (G[n][0]*D[2][5]*SQ2 + G[n][1]*D[2][4]*SQ2 + 2.0*G[n][2]*D[2][2]) / 2.0
		DσDun[3][0] = (G[n][0]*D[3][0]*SQ2 + G[n][2]*D[3][5] + G[n][1]*D[3][3]) / SQ2
		DσDun[3][1] = (G[n][1]*D[3][1]*SQ2 + G[n][2]*D[3][4] + G[n][0]*D[3][3]) / SQ2
		DσDun[3][2] = (G[n][2]*D[3][2]*SQ2 + G[n][0]*D[3][5] + G[n][1]*D[3][4]) / SQ2
		DσDun[4][0] = (G[n][0]*D[4][0]*SQ2 + G[n][2]*D[4][5] + G[n][1]*D[4][3]) / SQ2
		DσDun[4][1] = (G[n][1]*D[4][1]*SQ2 + G[n][2]*D[4][4] + G[n][0]*D[4][3]) / SQ2
		DσDun[4][2] = (G[n][2]*D[4][2]*SQ2 + G[n][0]*D[4][5] + G[n][1]*D[4][4]) / SQ2
		DσDun[5][0] = (G[n][0]*D[5][0]*SQ2 + G[n][2]*D[5][5] + G[n][1]*D[5][3]) / SQ2
		DσDun[5][1] = (G[n][1]*D[5][1]*SQ2 + G[n][2]*D[5][4] + G[n][0]*D[5][3]) / SQ2
		DσDun[5][2] = (G[n][2]*D[5][2]*SQ2 + G[n][0]*D[5][5] + G[n][1]*D[5][4]) / SQ2
	} else {
		DσDun[0][0] = (G[n][0]*D[0][0]*SQ2 + G[n][1]*D[0][3]) / SQ2
		DσDun[0][1] = (G[n][1]*D[0][1]*SQ2 + G[n][0]*D[0][3]) / SQ2
		DσDun[1][0] = (G[n][0]*D[1][0]*SQ2 + G[n][1]*D[1][3]) / SQ2
		DσDun[1][1] = (G[n][1]*D[1][1]*SQ2 + G[n][0]*D[1][3]) / SQ2
		DσDun[2][0] = (G[n][0]*D[2][0]*SQ2 + G[n][1]*D[2][3]) / SQ2
		DσDun[2][1] = (G[n][1]*D[2][1]*SQ2 + G[n][0]*D[2][3]) / SQ2
		DσDun[3][0] = (G[n][0]*D[3][0]*SQ2 + G[n][1]*D[3][3]) / SQ2
		DσDun[3][1] = (G[n][1]*D[3][1]*SQ2 + G[n][0]*D[3][3]) / SQ2
	}
}

func IpBmatrix(B [][]float64, ndim, nne int, G [][]float64, axisym bool, radius float64, S []float64) {
	if ndim == 3 {
		for i := 0; i < nne; i++ {
			B[0][0+i*3] = G[i][0]
			B[1][1+i*3] = G[i][1]
			B[2][2+i*3] = G[i][2]
			B[3][0+i*3] = G[i][1] / SQ2
			B[4][1+i*3] = G[i][2] / SQ2
			B[5][2+i*3] = G[i][0] / SQ2
			B[3][1+i*3] = G[i][0] / SQ2
			B[4][2+i*3] = G[i][1] / SQ2
			B[5][0+i*3] = G[i][2] / SQ2
		}
		return
	}
	if axisym {
		for i := 0; i < nne; i++ {
			B[0][0+i*2] = G[i][0]
			B[1][1+i*2] = G[i][1]
			B[2][0+i*2] = S[i] / radius
			B[3][0+i*2] = G[i][1] / SQ2
			B[3][1+i*2] = G[i][0] / SQ2
		}
		return
	}
	for i := 0; i < nne; i++ {
		B[0][0+i*2] = G[i][0]
		B[1][1+i*2] = G[i][1]
		B[3][0+i*2] = G[i][1] / SQ2
		B[3][1+i*2] = G[i][0] / SQ2
	}
}

func IpStrainsAndIncB(εs, Δεs []float64, nσ, nu int, B [][]float64, u, Δu []float64, Umap []int) {
	for i := 0; i < nσ; i++ {
		εs[i], Δεs[i] = 0, 0
		for j := 0; j < nu; j++ {
			εs[i] += B[i][j] * u[Umap[j]]
			Δεs[i] += B[i][j] * Δu[Umap[j]]
		}
	}
}

func IpBmatrix_sparse(B *la.Triplet, ndim, nne int, G [][]float64, axisym bool, radius float64, S []float64) {
	B.Start()
	if ndim == 3 {
		for i := 0; i < nne; i++ {
			B.Put(0, 0+i*3, G[i][0])
			B.Put(1, 1+i*3, G[i][1])
			B.Put(2, 2+i*3, G[i][2])
			B.Put(3, 0+i*3, G[i][1]/SQ2)
			B.Put(4, 1+i*3, G[i][2]/SQ2)
			B.Put(5, 2+i*3, G[i][0]/SQ2)
			B.Put(3, 1+i*3, G[i][0]/SQ2)
			B.Put(4, 2+i*3, G[i][1]/SQ2)
			B.Put(5, 0+i*3, G[i][2]/SQ2)
		}
		return
	}
	if axisym {
		for i := 0; i < nne; i++ {
			B.Put(0, 0+i*2, G[i][0])
			B.Put(1, 1+i*2, G[i][1])
			B.Put(2, 0+i*2, S[i]/radius)
			B.Put(3, 0+i*2, G[i][1]/SQ2)
			B.Put(3, 1+i*2, G[i][0]/SQ2)
		}
		return
	}
	for i := 0; i < nne; i++ {
		B.Put(0, 0+i*2, G[i][0])
		B.Put(1, 1+i*2, G[i][1])
		B.Put(3, 0+i*2, G[i][1]/SQ2)
		B.Put(3, 1+i*2, G[i][0]/SQ2)
	}
}

// Ivs2sigmas converts ivs map to a matrix with σ values [nip][nsig]
func Ivs2sigmas(nip, ndim int, ivs map[string][]float64) (σ [][]float64) {
	σ = la.MatAlloc(nip, 2*ndim)
	for key, vals := range ivs {
		for i := 0; i < nip; i++ {
			switch key {
			case "sx":
				σ[i][0] = vals[i]
			case "sy":
				σ[i][1] = vals[i]
			case "sz":
				σ[i][2] = vals[i]
			case "sxy":
				σ[i][3] = vals[i]
			case "syz":
				if ndim == 3 {
					σ[i][4] = vals[i]
				}
			case "szx":
				if ndim == 3 {
					σ[i][5] = vals[i]
				}
			}
		}
	}
	return
}

func StressKeys() []string {
	if Global.Ndim == 2 {
		return []string{"sx", "sy", "sz", "sxy"}
	}
	return []string{"sx", "sy", "sz", "sxy", "syz", "szx"}
}
