// Package provides functions to load electron cross sections from files in BOLSIG / LXCat format and to interpolate cross section functions.
// Energy is measured in electronvolts, cross sections in m^2.
package lxgata

import (
	"fmt"
	"strings"
)

type CollisionType string

const ELASTIC, EFFECTIVE, EXCITATION, ATTACHMENT, IONIZATION, ROTATION CollisionType = "ELASTIC", "EFFECTIVE", "EXCITATION", "ATTACHMENT", "IONIZATION", "ROTATION"

// Cross section point holds cross section value in [m^2] at energy [eV]
type CrossSectionPoint struct {
	Energy, Value float64
}

type Collision struct {
	Type            CollisionType
	MassRatio       float64 // ratio of electron mass to target particle, if applicable
	Species         string  // target particle species
	Data            []CrossSectionPoint
	Threshold       float64           // value of energy [eV], below which collision can not occur
	StatWeightRatio float64           // statistical weight ratio of the upper state to the lower state (for excitations)
	LowerEnergy     float64           // energy of lower state of rotational process (for rotations)
	LowerStatWeight float64           // statistical weight of lower state of rotational process (for rotations)
	UpperEnergy     float64           // energy of upper state of rotational process (for rotations)
	UpperStatWeight float64           // statistical weight of upper state of rotational process (for rotations)
	Info            map[string]string // any additional fields found in collision description
}

// CrossSectionAt calculates cross section at given energy as linear interpolation of piecewise linear cross section function.
// If the energy is below first or beyond last data point, it assumes cross section to be constant at corresponding values.
func (p *Collision) CrossSectionAt(energy float64) float64 {
	l, r := 0, len(p.Data)

	for c := (l + r) / 2; l+1 < r; c = (l + r) / 2 {
		if energy < p.Data[c].Energy {
			r = c
		} else {
			l = c
		}
	}
	if l == 0 || r == len(p.Data) {
		return p.Data[l].Value
	} else {
		w := (energy - p.Data[l].Energy) / (p.Data[r].Energy - p.Data[l].Energy)
		return p.Data[l].Value + (p.Data[r].Value-p.Data[l].Value)*w
	}
}

func (p Collision) String() string {
	return fmt.Sprintf("Cross section of %v %v. Threshold: %v", p.Species, strings.ToLower(string(p.Type)), p.Threshold)
}
