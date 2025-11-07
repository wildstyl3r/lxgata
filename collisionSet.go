// Package provides functions to load electron cross sections from files in BOLSIG / LXCat format and to interpolate cross section functions.
// Energy is measured in electronvolts, cross sections in m^2.
package lxgata

import (
	"bufio"
	"fmt"
	"math"
	"math/rand"
	"os"
	"slices"
	"sort"
	"strconv"
	"strings"
)

type ScatteringMode int

const Hartree float64 = 27.0211386 // [eV]

const (
	Incorrect ScatteringMode = iota
	Isotropic
	Coulomb
	Born
)

type Collisions struct {
	elasticScatteringMode              ScatteringMode
	inelasticScatteringMode            ScatteringMode
	UParameter                         float64 // parameter regulating the shape of differential cross section in the Coulomb model, u_eta in Hagelaar's MCIG paper, equation (32) at page 10
	Processes                          []Collision
	TotalCrossSectionAtCache           []float64
	TotalCrossSectionEnergyStep        float64
	TotalCrossSectionEnergyStepInverse float64
	TotalCrossSectionUpTo              float64
}

// LoadCrossSections loads cross section data from file in LXCat/BOLSIG format
func LoadCrossSections(fileName string, forMonteCarlo bool, totalCrossSectionEnergyStep, totalCrossSectionUpTo float64, elasticScatteringMode, inelasticScatteringMode ScatteringMode, uParameter float64, z AtomicNumber) (Collisions, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return Collisions{}, err
	}
	defer file.Close()

	setProcessTypes := map[string]struct{}{string(ELASTIC): {}, string(EFFECTIVE): {}, string(EXCITATION): {}, string(ATTACHMENT): {}, string(IONIZATION): {}, string(ROTATION): {}}

	var collisions = Collisions{
		elasticScatteringMode:   elasticScatteringMode,
		inelasticScatteringMode: inelasticScatteringMode,
		UParameter:              uParameter,
	}

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
			species, outcome, _ := strings.Cut(scanner.Text(), " ")

			scanner.Scan()
			parameters := strings.Fields(strings.Trim(scanner.Text(), " "))

			var massRatio, threshold, lowerEnergy, upperEnergy, lowerStatWeight, upperStatWeight float64
			statWeightRatio := 1.

			switch collisionType {
			case ELASTIC:
				massRatio, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return Collisions{}, err
				}
			case EFFECTIVE:
				massRatio, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return Collisions{}, err
				}
			case EXCITATION:
				threshold, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return Collisions{}, err
				}
				if len(parameters) > 1 {
					statWeightRatio, err = strconv.ParseFloat(parameters[1], 64)
					if err != nil {
						return Collisions{}, err
					}
				}
			case IONIZATION:
				threshold, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return Collisions{}, err
				}
			case ROTATION:
				lowerEnergy, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return Collisions{}, err
				}
				lowerStatWeight, err = strconv.ParseFloat(parameters[1], 64)
				if err != nil {
					return Collisions{}, err
				}
				scanner.Scan()
				parameters = strings.Fields(scanner.Text())
				upperEnergy, err = strconv.ParseFloat(parameters[0], 64)
				if err != nil {
					return Collisions{}, err
				}
				upperStatWeight, err = strconv.ParseFloat(parameters[1], 64)
				if err != nil {
					return Collisions{}, err
				}
			}

			excitationType := Ordinary
			info := make(map[string]string)
			for !strings.HasPrefix(scanner.Text(), "-----") {
				key, val, found := strings.Cut(scanner.Text(), ":")
				if found {
					info[strings.Trim(key, " ")] = strings.Trim(val, " ")
					if strings.Contains(val, "VIB") {
						if excitationType == Ordinary || excitationType == Vibrational {
							excitationType = Vibrational
						} else {
							println("Ambiguity in excitation type: [", key, ": ", val, "], already set as Rotational")
						}
					}
					if strings.Contains(val, "ROT") {
						if excitationType == Ordinary || excitationType == Rotational {
							excitationType = Rotational
						} else {
							println("Ambiguity in excitation type: [", key, ": ", val, "], already set as Vibrational")
						}
					}
				}
				scanner.Scan()
			}
			scanner.Scan()

			var data []CrossSectionPoint
			for !strings.HasPrefix(scanner.Text(), "-----") {
				crossSectionPoint := strings.Fields(scanner.Text())
				energy, err := strconv.ParseFloat(crossSectionPoint[0], 64)
				if err != nil {
					return Collisions{}, err
				}

				crossSection, err := strconv.ParseFloat(crossSectionPoint[1], 64)
				if err != nil {
					return Collisions{}, err
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

			collisions.Processes = append(collisions.Processes, Collision{
				Type:            collisionType,
				Excitation:      excitationType,
				MassRatio:       massRatio,
				Species:         species,
				Outcome:         outcome,
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
	if forMonteCarlo {
		for i := range len(collisions.Processes) {
			if collisions.Processes[i].Type == EFFECTIVE {
				collisions.Processes = append(append(collisions.Processes[:i], Collision{
					Type:      ELASTIC,
					MassRatio: collisions.Processes[i].MassRatio,
					Species:   collisions.Processes[i].Species,
					Data:      collisions.CalculateElasticFromEffective(),
					Info:      collisions.Processes[i].Info,
				}), collisions.Processes[i+1:]...)
				break
			}
		}
	}
	switch collisions.elasticScatteringMode {
	case Incorrect:
		panic("elastic scattering mode not set")
	case Coulomb:
		for i := range len(collisions.Processes) {
			if collisions.Processes[i].Type == ELASTIC {
				for j := range collisions.Processes[i].Data {
					if j == 0 {
						break
					}
					collisions.Processes[i].Data[j].Value /= CoulombNormalization(collisions.Processes[i].Data[j].Energy, 0., Hartree, z)
				}
				break
			}
		}

	case Born:
	case Isotropic:
	}

	if totalCrossSectionEnergyStep != 0 {
		tcsCache := make([]float64, int(totalCrossSectionUpTo/totalCrossSectionEnergyStep))
		for i := range tcsCache {
			tcsCache[i] = collisions.TotalCrossSectionAt((float64(i) + 0.5) * totalCrossSectionEnergyStep)
		}
		collisions.TotalCrossSectionUpTo = totalCrossSectionUpTo
		collisions.TotalCrossSectionAtCache = tcsCache
		collisions.TotalCrossSectionEnergyStep = totalCrossSectionEnergyStep
		collisions.TotalCrossSectionEnergyStepInverse = 1. / totalCrossSectionEnergyStep

	}

	return collisions, nil
}

func Argsort(s []float64, abs bool) (indices []int) {
	indices = make([]int, len(s))
	for i := range indices {
		indices[i] = i
	}

	if abs {
		sort.Slice(indices, func(i, j int) bool {
			return math.Abs(s[indices[i]]) < math.Abs(s[indices[j]])
		})
	} else {
		sort.Slice(indices, func(i, j int) bool {
			return s[indices[i]] < s[indices[j]]
		})
	}

	return indices
}

func SumFloat64Slice(arr []float64) (sum float64) { // Kahan's algorithm
	compensation := 0.
	summationOrder := Argsort(arr, true)
	for i := range summationOrder {
		y := arr[summationOrder[i]] - compensation
		temp := sum + y
		compensation = (temp - sum) - y
		sum = temp
	}
	return sum
}

func (colls Collisions) CalculateElasticFromEffective() []CrossSectionPoint {
	effectiveIndex := -1
	for i := range colls.Processes {
		if colls.Processes[i].Type == EFFECTIVE {
			effectiveIndex = i
		}
	}
	if effectiveIndex == -1 {
		return nil
	}

	elasticData := make([]CrossSectionPoint, len(colls.Processes[effectiveIndex].Data))
	for i := range colls.Processes[effectiveIndex].Data {
		inelasticSumTerms := make([]float64, 0, len(colls.Processes)-1)
		energy := colls.Processes[effectiveIndex].Data[i].Energy
		for process := range colls.Processes {
			if colls.Processes[process].Type != EFFECTIVE && colls.Processes[process].Type != ELASTIC {
				inelasticSumTerms = append(inelasticSumTerms, colls.Processes[process].CrossSectionAt(energy))
			}
		}
		elasticData[i] = CrossSectionPoint{Energy: energy, Value: colls.Processes[effectiveIndex].Data[i].Value - SumFloat64Slice(inelasticSumTerms)}
	}
	for i := 0; i+1 < len(elasticData); i++ {
		elasticData[i]._NextValDiffPerEnergyDiff = (elasticData[i+1].Value - elasticData[i].Value) / (elasticData[i+1].Energy - elasticData[i].Energy)
	}
	return elasticData
}

// TotalCrossSectionAt returns total cross section at given energy for all species and processes in collisions set
func (colls Collisions) TotalCrossSectionAt(energy float64) float64 {
	if energy < colls.TotalCrossSectionUpTo {
		return colls.TotalCrossSectionAtCache[int(energy*colls.TotalCrossSectionEnergyStepInverse)]
	}

	var result float64
	for i := range colls.Processes {
		result += colls.Processes[i].CrossSectionAt(energy)
	}
	return result
}

func (colls Collisions) CrossSectionsAt(energy float64) (result []float64) {
	result = make([]float64, len(colls.Processes))
	for i := range colls.Processes {
		result[i] = colls.Processes[i].CrossSectionAt(energy)
	}
	return
}

// TotalCrossSectionOfKindAt returns summed cross section of given type at given energy for all species and processes in collision set
func (colls Collisions) TotalCrossSectionOfKindAt(t CollisionType, energy float64) float64 {
	var result float64
	for i := range colls.Processes {
		if colls.Processes[i].Type == t {
			result += colls.Processes[i].CrossSectionAt(energy)
		}
	}
	return result
}

func (colls Collisions) MinThreshold() float64 {
	var result = math.MaxFloat64
	for i := range colls.Processes {
		if colls.Processes[i].Threshold != 0 && colls.Processes[i].Threshold < result {
			result = colls.Processes[i].Threshold
		}
	}
	return result
}

func (colls Collisions) MinThresholdOfKind(t CollisionType) float64 {
	var result = math.MaxFloat64
	for i := range colls.Processes {
		if colls.Processes[i].Threshold != 0 && colls.Processes[i].Threshold < result && colls.Processes[i].Type == t {
			result = colls.Processes[i].Threshold
		}
	}
	return result
}

// SurplusCrossSection returns sum of maximum values of cross sections over all processes in collision set
// Can be used to estimate lower bound on mean free path
func (colls Collisions) SurplusCrossSection() float64 {
	var result float64
	for i := range colls.Processes {
		var max float64
		for d := range colls.Processes[i].Data {
			if max < colls.Processes[i].Data[d].Value {
				max = colls.Processes[i].Data[d].Value
			}
		}
		result += max
	}
	return result
}

func (colls Collisions) MakeEnergyGrid(minStep, maxEnergy float64) []float64 {
	for i := range colls.Processes {
		maxEnergy = max(maxEnergy, colls.Processes[i].Data[len(colls.Processes[i].Data)-1].Energy)
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

func (colls Collisions) SampleScatteringAngleCos(energy, transitionEnergy float64, collisionType CollisionType, z AtomicNumber) (cosChi float64) {
	if collisionType == ELASTIC {
		switch colls.elasticScatteringMode {
		case Born:
			panic("lxgata does not support Born elastic scattering\n")
		case Coulomb:
			return CoulombScatteringAngleSample(energy, colls.UParameter, transitionEnergy, z)
		case Isotropic:
			return 1. - 2.*rand.Float64()
		default:
			panic(fmt.Sprintf("unexpected lxgata.ScatteringMode: %#v", colls.inelasticScatteringMode))
		}
	} else {
		switch colls.inelasticScatteringMode {
		case Born:
			return BornScatteringAngleSample(energy, transitionEnergy)
		case Coulomb:
			return CoulombScatteringAngleSample(energy, colls.UParameter, transitionEnergy, z)
		case Isotropic:
			return 1. - 2.*rand.Float64()
		default:
			panic(fmt.Sprintf("unexpected lxgata.ScatteringMode: %#v", colls.inelasticScatteringMode))
		}
	}
}

func (colls Collisions) GetTypes() (types []CollisionType) {
	set := map[CollisionType]struct{}{}
	for i := range colls.Processes {
		set[colls.Processes[i].Type] = struct{}{}
	}
	for collType := range set {
		types = append(types, collType)
	}
	slices.Sort(types)
	return types
}
