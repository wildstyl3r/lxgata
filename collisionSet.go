// Package provides functions to load electron cross sections from files in BOLSIG / LXCat format and to interpolate cross section functions.
// Energy is measured in electronvolts, cross sections in m^2.
package lxgata

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
)

type Collisions []Collision

// LoadCrossSections loads cross section data from file in LXCat/BOLSIG format
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
					data = append(data, CrossSectionPoint{energy, crossSection, 0.})
				}
				scanner.Scan()
			}

			if (collisionType == IONIZATION || collisionType == EXCITATION || collisionType == ROTATION) && data[0].Value > 0. {
				data = append([]CrossSectionPoint{{threshold, 0., 0.}}, data...)
			}

			for i := 0; i+1 < len(data); i++ {
				data[i]._NextValDiffPerEnergyDiff = (data[i+1].Value - data[i].Value) / (data[i+1].Energy - data[i].Energy)
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

// TotalCrossSectionAt returns total cross section at given energy for all species and processes in collisions set
func (colls Collisions) TotalCrossSectionAt(energy float64) float64 {
	var result float64
	for i := range colls {
		result += colls[i].CrossSectionAt(energy)
	}
	return result
}

func (colls Collisions) CrossSectionsAt(energy float64) (result []float64) {
	result = make([]float64, len(colls))
	for i := range colls {
		result[i] = colls[i].CrossSectionAt(energy)
	}
	return
}

// TotalCrossSectionOfKindAt returns summed cross section of given type at given energy for all species and processes in collision set
func (colls Collisions) TotalCrossSectionOfKindAt(t CollisionType, energy float64) float64 {
	var result float64
	for i := range colls {
		if colls[i].Type == t {
			result += colls[i].CrossSectionAt(energy)
		}
	}
	return result
}

func (colls Collisions) MinThreshold() float64 {
	var result = math.MaxFloat64
	for i := range colls {
		if colls[i].Threshold != 0 && colls[i].Threshold < result {
			result = colls[i].Threshold
		}
	}
	return result
}

func (colls Collisions) MinThresholdOfKind(t CollisionType) float64 {
	var result = math.MaxFloat64
	for i := range colls {
		if colls[i].Threshold != 0 && colls[i].Threshold < result && colls[i].Type == t {
			result = colls[i].Threshold
		}
	}
	return result
}

// SurplusCrossSection returns sum of maximum values of cross sections over all processes in collision set
// Can be used to estimate lower bound on mean free path
func (colls Collisions) SurplusCrossSection() float64 {
	var result float64
	for i := range colls {
		var max float64
		for d := range colls[i].Data {
			if max < colls[i].Data[d].Value {
				max = colls[i].Data[d].Value
			}
		}
		result += max
	}
	return result
}

func (colls Collisions) MakeEnergyGrid(minStep, maxEnergy float64) []float64 {
	for i := range colls {
		maxEnergy = max(maxEnergy, colls[i].Data[len(colls[i].Data)-1].Energy)
	}
	nSteps := int(maxEnergy / minStep)
	finestGrid := make([]float64, nSteps)
	for i := range finestGrid {
		finestGrid[i] = colls.TotalCrossSectionAt(float64(i) * minStep)
	}
	maxCSChange := 0.
	for i := range nSteps - 1 {
		maxCSChange = max(maxCSChange, math.Abs(finestGrid[i+1]-finestGrid[i]))
	}
	peaks := map[int]struct{}{}
	for i := range nSteps - 2 {
		if (finestGrid[i+1]-finestGrid[i])*(finestGrid[i+2]-finestGrid[i+1]) < 0 {
			peaks[i+1] = struct{}{}
		}
	}

	gridIndicies := []int{0}
	for i := range nSteps - 1 {
		if _, peak := peaks[i+1]; peak {
			gridIndicies = append(gridIndicies, i+1)
			continue
		}
		lastIndex := gridIndicies[len(gridIndicies)-1]
		if i != lastIndex && math.Abs(finestGrid[i+1]-finestGrid[lastIndex]) > maxCSChange {
			gridIndicies = append(gridIndicies, i)
		}
	}
	gridIndicies = append(gridIndicies, nSteps-1)

	grid := make([]float64, len(gridIndicies))
	for i := range gridIndicies {
		grid[i] = float64(gridIndicies[i]) * minStep
	}
	return grid
}
