// Package provides functions to load electron cross sections from files in BOLSIG / LXCat format and to interpolate cross section functions.
// Energy is measured in electronvolts, cross sections in m^2.
package lxgata

import (
	"fmt"
	"strings"
)

type CollisionType string

const ELASTIC, EFFECTIVE, EXCITATION, ATTACHMENT, IONIZATION, ROTATION, DEEXCITATION, ION_ISOTROPIC, ION_BACKSCATTERING CollisionType = "ELASTIC", "EFFECTIVE", "EXCITATION", "ATTACHMENT", "IONIZATION", "ROTATION", "DEEXCITATION", "ION_ISOTROPIC", "ION_BACKSCATTERING"

// Cross section point holds cross section value in [m^2] at energy [eV]
type CrossSectionPoint struct {
	Energy, Value             float64
	_NextValDiffPerEnergyDiff float64
}

type ExcitationType int

const (
	Ordinary ExcitationType = iota
	Rotational
	Vibrational
)

type Collision struct {
	Type            CollisionType
	Excitation      ExcitationType
	MassRatio       float64 // ratio of electron mass to target particle, if applicable
	Threshold       float64 // value of energy [eV], below which collision can not occur
	Data            []CrossSectionPoint
	StatWeightRatio float64           // statistical weight ratio of the upper state to the lower state (for excitations)
	LowerEnergy     float64           // energy of lower state of rotational process (for rotations)
	LowerStatWeight float64           // statistical weight of lower state of rotational process (for rotations)
	UpperEnergy     float64           // energy of upper state of rotational process (for rotations)
	UpperStatWeight float64           // statistical weight of upper state of rotational process (for rotations)
	Info            map[string]string // any additional fields found in collision description
	Species         string            // target particle species
	SpeciesMass     float64
	Outcome         string
}

// CrossSectionAt calculates cross section at given energy as linear interpolation of piecewise linear cross section function.
// If the energy is below first or beyond last data point, it assumes cross section to be constant at corresponding values.
func (p *Collision) CrossSectionAt(energy float64) (value float64) {
	if p.Type == DEEXCITATION {
		//threshold was set as negative to that of forward process for simplicity
		if energy < 1e-8 {
			return 0
		}
		factor := p.StatWeightRatio * (energy - p.Threshold) / energy
		energy -= p.Threshold
		var l, r uint = 0, uint(len(p.Data))

		for c := (l + r) / 2; l+1 < r; c = (l + r) / 2 {
			if energy < p.Data[c].Energy {
				r = c
			} else {
				l = c
			}
		}
		if l == 0 || r == uint(len(p.Data)) {
			value = p.Data[l].Value
		} else {
			value = p.Data[l].Value + (energy-p.Data[l].Energy)*p.Data[l]._NextValDiffPerEnergyDiff
		}
		return factor * value
	} else {
		var l, r uint = 0, uint(len(p.Data))

		for c := (l + r) / 2; l+1 < r; c = (l + r) / 2 {
			if energy < p.Data[c].Energy {
				r = c
			} else {
				l = c
			}
		}
		if l == 0 || r == uint(len(p.Data)) {
			return p.Data[l].Value
		} else {
			return p.Data[l].Value + (energy-p.Data[l].Energy)*p.Data[l]._NextValDiffPerEnergyDiff
		}
	}
}

func (p Collision) String() string {
	return fmt.Sprintf("Cross section of %v %v. Threshold: %v", p.Species, strings.ToLower(string(p.Type)), p.Threshold)
}
