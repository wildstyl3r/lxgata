package lxgata

import (
	"math"
	"math/rand"
)

// following Hagelaar's MCIG paper

type AtomicNumber int

func CoulombScatteringAngleSample(energy, uParameter, transitionEnergy float64, z AtomicNumber) (cosChi float64) {
	r := rand.Float64()
	eta := CoulombEta(energy, transitionEnergy, uParameter, z)
	return 1. - 2.*eta*r/(eta+1-r)
}

func CoulombEta(energy, transitionEnergy, uParameter float64, z AtomicNumber) float64 {
	if z == 0 {
		beta := math.Sqrt(1 - transitionEnergy/energy)
		return uParameter / (8. * beta * beta * energy)
	} else {
		return 1.89 * math.Pow(float64(z), 2./3.) / energy
	}
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
