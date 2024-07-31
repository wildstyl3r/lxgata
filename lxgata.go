// Package provides functions to load electron cross sections from files in BOLSIG / LXCat format and to interpolate cross section functions.
// Energy is measured in electronvolts, cross sections in m^2.
package lxgata

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
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

type Collisions []Collision

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

func LoadCrossSections(fileName string) (Collisions, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	setProcessTypes := map[string]struct{}{string(ELASTIC): {}, string(EFFECTIVE): {}, string(EXCITATION): {}, string(ATTACHMENT): {}, string(IONIZATION): {}, string(ROTATION): {}}

	var collisions []Collision

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		if len(tokens) == 0 {
			continue
		}

		if _, ok := setProcessTypes[tokens[0]]; ok {
			collisionType := CollisionType(tokens[0])

			scanner.Scan()
			species, _, _ := strings.Cut(scanner.Text(), " ")

			scanner.Scan()
			parameters := strings.Fields(strings.Trim(scanner.Text(), " "))

			var massRatio, threshold, lowerEnergy, upperEnergy, lowerStatWeight, upperStatWeight float64
			statWeightRatio := 1.

			switch collisionType {
			case ELASTIC:
				massRatio, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return nil, err
				}
			case EFFECTIVE:
				massRatio, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return nil, err
				}
			case EXCITATION:
				threshold, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return nil, err
				}
				if len(parameters) > 1 {
					statWeightRatio, err = strconv.ParseFloat(parameters[1], 64)
					if err != nil {
						return nil, err
					}
				}
			case IONIZATION:
				threshold, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return nil, err
				}
			case ROTATION:
				lowerEnergy, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return nil, err
				}
				lowerStatWeight, err = strconv.ParseFloat(parameters[1], 64)
				if err != nil {
					return nil, err
				}
				scanner.Scan()
				parameters = strings.Fields(scanner.Text())
				upperEnergy, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return nil, err
				}
				upperStatWeight, err = strconv.ParseFloat(parameters[1], 64)
				if err != nil {
					return nil, err
				}
			}

			info := make(map[string]string)
			for !strings.HasPrefix(scanner.Text(), "-----") {
				key, val, found := strings.Cut(scanner.Text(), ":")
				if found {
					info[strings.Trim(key, " ")] = strings.Trim(val, " ")
				}
				scanner.Scan()
			}
			scanner.Scan()

			var data []CrossSectionPoint
			for !strings.HasPrefix(scanner.Text(), "-----") {
				crossSectionPoint := strings.Fields(scanner.Text())
				energy, err := strconv.ParseFloat(crossSectionPoint[0], 64)
				if err != nil {
					return nil, err
				}

				crossSection, err := strconv.ParseFloat(crossSectionPoint[1], 64)
				if err != nil {
					return nil, err
				}
				if !(collisionType == IONIZATION || collisionType == EXCITATION || collisionType == ROTATION) || threshold < energy {
					data = append(data, CrossSectionPoint{energy, crossSection})
				}
				scanner.Scan()
			}

			if (collisionType == IONIZATION || collisionType == EXCITATION || collisionType == ROTATION) && data[0].Value > 0. {
				data = append([]CrossSectionPoint{{threshold, 0.}}, data...)
			}

			collisions = append(collisions, Collision{
				Type:            collisionType,
				MassRatio:       massRatio,
				Species:         species,
				Data:            data,
				Threshold:       threshold,
				StatWeightRatio: statWeightRatio,
				LowerEnergy:     lowerEnergy,
				LowerStatWeight: lowerStatWeight,
				UpperEnergy:     upperEnergy,
				UpperStatWeight: upperStatWeight,
				Info:            info,
			})
		}
	}
	return collisions, nil
}

func (colls Collisions) TotalCrossSectionAt(energy float64) float64 {
	var result float64
	for _, collision := range colls {
		result += collision.CrossSectionAt(energy)
	}
	return result
}
