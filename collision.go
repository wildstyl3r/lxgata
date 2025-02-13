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
	Energy, Value             float64
	_NextValDiffPerEnergyDiff float64
}

type Collision struct {
	Type            CollisionType
	MassRatio       float64 // ratio of electron mass to target particle, if applicable
	Species         string  // target particle species
	ExpandedData    []float64
	ExpandedDiff    float64
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
	if p.ExpandedData == nil {
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
	} else {
		energy -= p.Threshold
		if energy < 0 {
			return 0.
		}
		index := int(energy / float64(p.ExpandedDiff))
		if index >= len(p.ExpandedData) {
			return float64(p.ExpandedData[len(p.ExpandedData)-1])
		}
		return p.ExpandedData[index] + (energy-p.ExpandedDiff*float64(index))/p.ExpandedDiff*(p.ExpandedData[index+1]-p.ExpandedData[index])
	}
}

func (p Collision) String() string {
	return fmt.Sprintf("Cross section of %v %v. Threshold: %v", p.Species, strings.ToLower(string(p.Type)), p.Threshold)
}

func (p *Collision) Expand(diff float64) {
	p.ExpandedDiff = diff
	p.ExpandedData = make([]float64, int((p.Data[len(p.Data)-1].Energy-p.Threshold)/diff)+1)
	dataIndex := 0
	var share, curEnergy float64
	for expIndex := 0; expIndex < len(p.ExpandedData); expIndex++ {
		curEnergy = p.Threshold + float64(expIndex)*diff
		//?? dataIndex < len(p.Data) &&
		for dataIndex+1 < len(p.Data) && curEnergy > p.Data[dataIndex+1].Energy {
			dataIndex++
		}
		if dataIndex+1 == len(p.Data) {
			p.ExpandedData[expIndex] = p.Data[dataIndex].Value
		} else {
			share = (curEnergy - p.Data[dataIndex].Energy) / (p.Data[dataIndex+1].Energy - p.Data[dataIndex].Energy)
			p.ExpandedData[expIndex] = p.Data[dataIndex].Value + share*(p.Data[dataIndex+1].Value-p.Data[dataIndex].Value)
		}
	}
	// for dataIndex+1 < len(p.Data) {
	// 	for p.Threshold+expIndex*diff < p.Data[dataIndex+1].Energy && expIndex < len(p.ExpandedData) {
	// 		cs := float64(1.?)
	// 		p.ExpandedData[expIndex] = cs
	// 		expIndex++
	// 	}
	// 	dataIndex++
	// }
}
