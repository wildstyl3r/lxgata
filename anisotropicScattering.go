package lxgata

import (
	"math"
	"math/rand"
)

// following Hagelaar's MCIG paper

func CoulombScatteringAngleSample(energy, uParameter, transitionEnergy float64) (cosChi float64) {
	beta := math.Sqrt(1 - transitionEnergy/energy)
	eta := uParameter / (8. * beta * beta * energy)
	r := rand.Float64()
	return 1. - 2.*eta*r/(eta+1-r)
}

func BornScatteringAngleSample(energy, transitionEnergy float64) (cosChi float64) {
	beta := math.Sqrt(1 - transitionEnergy/energy)
	r := rand.Float64()
	return 1 - (1-beta)*(1-beta)/(2.*beta)*
		(math.Pow(math.Abs((1.+beta)/(1-beta)), 2.*r)-1.)
}

func CoulombNormalization(energy, transitionEnergy, uParameter float64) (denominator float64) {
	beta := math.Sqrt(1 - transitionEnergy/energy)
	eta := uParameter / (8 * beta * beta * energy)
	return 2 * eta * ((1+eta)*math.Log(1.+1./eta) - 1)
}

func BornDifferentialCrossSection(beta, chi float64) float64 {
	return beta / (2 * math.Pi * math.Log(math.Abs((1+beta)/(1-beta))) * (1 + beta*beta - 2*beta*math.Cos(chi)))
}

func BornNormalization(energy, transitionEnergy float64, anglePieces int) (denominator float64) {
	beta := math.Sqrt(1 - transitionEnergy/energy)
	//using trapezoid method
	firstBase := 0.
	deltaChi := math.Pi / float64(anglePieces)
	sumList := make([]float64, anglePieces)
	for i := 1; i <= anglePieces; i++ {
		angle := float64(i) * deltaChi
		secondBase := (1. - math.Cos(angle)) * BornDifferentialCrossSection(beta, angle) * math.Sin(angle)
		sumList[i-1] = 0.5 * (firstBase + secondBase)
		firstBase = secondBase
	}
	return SumFloat64Slice(sumList)
}
