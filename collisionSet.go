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

	"github.com/pokeyaro/algo-ds-go/ds/tree/segment"
)

type ScatteringMode int

const Hartree float64 = 27.0211386 // [eV]

const (
	Incorrect ScatteringMode = iota
	Isotropic
	Coulomb
	Born
)

type Species struct {
	ShareOfUnity float64
	UParameter   float64
}

type Collisions struct {
	elasticScatteringMode   ScatteringMode
	inelasticScatteringMode ScatteringMode
	// UParameter               float64 // parameter regulating the shape of differential cross section in the Coulomb model, u_eta in Hagelaar's MCIG paper, equation (32) at page 10
	Processes                []Collision
	Species                  map[string]Species
	TotalCrossSectionAtCache []float64
	MaxTCSRangeTree          *segment.SegmentTree[float64]
	FixedStepTable           [][]float64
	EnergyStep               float64
	EnergyStepInverse        float64
	TotalCrossSectionUpTo    float64
	MaxTCS                   float64
}

type maxOverFloat64 struct{}

func (maxOverFloat64) Merge(a, b float64) float64 { return max(a, b) }

// LoadCrossSections loads cross section data from file in LXCat/BOLSIG format
func LoadCrossSections(fileName string, forMonteCarlo bool, totalCrossSectionEnergyStep, totalCrossSectionUpTo float64, elasticScatteringMode, inelasticScatteringMode ScatteringMode, uParameter float64, z AtomicNumber, species map[string]Species) (Collisions, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return Collisions{}, err
	}
	defer file.Close()

	setProcessTypes := map[string]struct{}{string(ELASTIC): {}, string(EFFECTIVE): {}, string(EXCITATION): {}, string(ATTACHMENT): {}, string(IONIZATION): {}, string(ROTATION): {}}

	var collisions = Collisions{
		elasticScatteringMode:   elasticScatteringMode,
		inelasticScatteringMode: inelasticScatteringMode,
		Species:                 species,
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
			species, outcome, _ := strings.Cut(scanner.Text(), "->")
			species, makeInverse := strings.CutSuffix(species, "<")
			species = strings.Trim(species, " ")
			outcome = strings.Trim(outcome, " ")

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

			if makeInverse && collisions.Processes[len(collisions.Processes)-1].Type == EXCITATION {
				collisions.Processes = append(collisions.Processes, Collision{
					Type:            DEEXCITATION,
					Excitation:      excitationType,
					MassRatio:       massRatio,
					Species:         outcome,
					Outcome:         species,
					Data:            data,
					Threshold:       -threshold,
					StatWeightRatio: 1. / statWeightRatio,
				})
			}
		}
	}
	if forMonteCarlo {
		for i := range len(collisions.Processes) {
			if collisions.Processes[i].Type == EFFECTIVE {
				collisions.Processes = append(append(collisions.Processes[:i], Collision{
					Type:      ELASTIC,
					MassRatio: collisions.Processes[i].MassRatio,
					Species:   collisions.Processes[i].Species,
					Data:      collisions.CalculateElasticFromEffective(collisions.Processes[i].Species),
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
		for p := range len(collisions.Processes) {
			if collisions.Processes[p].Type == ELASTIC {
				for j := range collisions.Processes[p].Data {
					if j == 0 {
						break
					}
					collisions.Processes[p].Data[j].Value /= CoulombNormalization(collisions.Processes[p].Data[j].Energy, 0., Hartree, z)
				}
				for i := 0; i+1 < len(collisions.Processes[p].Data); i++ {
					collisions.Processes[p].Data[i]._NextValDiffPerEnergyDiff = (collisions.Processes[p].Data[i+1].Value - collisions.Processes[p].Data[i].Value) / (collisions.Processes[p].Data[i+1].Energy - collisions.Processes[p].Data[i].Energy)
				}
				break
			}
		}

	case Born:
	case Isotropic:
	}

	if totalCrossSectionEnergyStep != 0 && totalCrossSectionUpTo != 0 {
		thresholdUpperCells := make([]int, len(collisions.Processes))
		for p := range collisions.Processes {
			thresholdUpperCells[p] = int(math.Ceil(collisions.Processes[p].Threshold / totalCrossSectionEnergyStep))
		}

		numberOfSteps := int(totalCrossSectionUpTo / totalCrossSectionEnergyStep)
		fixedStepTable := make([][]float64, numberOfSteps)
		tcsCache := make([]float64, numberOfSteps)
		for step := range fixedStepTable {
			fixedStepTable[step] = make([]float64, len(collisions.Processes))
			for process := range fixedStepTable[step] {
				share := collisions.Species[collisions.Processes[process].Species].ShareOfUnity
				if step != thresholdUpperCells[process] {
					fixedStepTable[step][process] = share * collisions.Processes[process].CrossSectionAt(float64(step)*totalCrossSectionEnergyStep)
				} else {
					fixedStepTable[step][process] = 0
				}
			}
			tcsCache[step] = SumFloat64Slice(fixedStepTable[step])
		}
		collisions.FixedStepTable = fixedStepTable
		collisions.TotalCrossSectionUpTo = totalCrossSectionUpTo
		collisions.TotalCrossSectionAtCache = tcsCache
		collisions.EnergyStep = totalCrossSectionEnergyStep
		collisions.EnergyStepInverse = 1. / totalCrossSectionEnergyStep
		if len(collisions.TotalCrossSectionAtCache) > 0 {
			collisions.MaxTCS = slices.Max(collisions.TotalCrossSectionAtCache)
			collisions.MaxTCSRangeTree = segment.NewSegmentTree(collisions.TotalCrossSectionAtCache, maxOverFloat64{})
		}
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

func (colls *Collisions) CalculateElasticFromEffective(species string) []CrossSectionPoint {
	effectiveIndex := -1
	for i := range colls.Processes {
		if colls.Processes[i].Type == EFFECTIVE && colls.Processes[i].Species == species {
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
			if colls.Processes[process].Type != EFFECTIVE && colls.Processes[process].Type != ELASTIC && colls.Processes[process].Species == species {
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
func (colls *Collisions) TotalCrossSectionAt(energy float64) float64 {
	step := math.Floor(energy * colls.EnergyStepInverse)
	if i := int(step); energy < colls.TotalCrossSectionUpTo && i+1 < len(colls.TotalCrossSectionAtCache) {
		if step < 1e-3 {
			return 0
		}
		delta := energy - step*colls.EnergyStep

		return linearInterpolation(colls.TotalCrossSectionAtCache[i], colls.TotalCrossSectionAtCache[i+1], delta/step)
	}

	var result float64
	for i := range colls.Processes {
		share := colls.Species[colls.Processes[i].Species].ShareOfUnity
		result += share * colls.Processes[i].CrossSectionAt(energy)
	}
	return result
}

func linearInterpolation(a, b, t float64) float64 {
	if t < 0.5 {
		return math.FMA(b-a, t, a)
	} else {
		return math.FMA(b-a, 1-t, b)
	}
}

func (colls *Collisions) SampleWithNullCollision(energy, totalCSPrimed float64) *Collision {
	choice := rand.Float64() * totalCSPrimed
	accum := 0.
	step := math.Floor(energy * colls.EnergyStepInverse)
	if stepI := int(step); energy < colls.TotalCrossSectionUpTo && stepI+1 < len(colls.FixedStepTable) {
		t := (energy - step*colls.EnergyStep) / colls.EnergyStep
		for p := range colls.Processes {
			accum += linearInterpolation(colls.FixedStepTable[stepI][p], colls.FixedStepTable[stepI+1][p], t)
			if choice < accum {
				return &colls.Processes[p]
			}
		}
	} else {
		for p := range colls.Processes {
			share := colls.Species[colls.Processes[p].Species].ShareOfUnity
			accum += share * colls.Processes[p].CrossSectionAt(energy)
			if choice < accum {
				return &colls.Processes[p]
			}
		}
	}
	return nil
}

func (colls *Collisions) CrossSectionsAt(energy float64) (result []float64) {
	step := math.Floor(energy * colls.EnergyStepInverse)
	if i := int(step); energy < colls.TotalCrossSectionUpTo && i+1 < len(colls.FixedStepTable) {
		result = make([]float64, len(colls.FixedStepTable[i]))
		if step < 1e-3 {
			return result
		}
		t := (energy - step*colls.EnergyStep) / step
		for j := range result {
			result[j] = linearInterpolation(colls.FixedStepTable[i][j], colls.FixedStepTable[i+1][j], t)
		}
		return
	}
	result = make([]float64, len(colls.Processes))
	for p := range colls.Processes {
		share := colls.Species[colls.Processes[p].Species].ShareOfUnity
		result[p] = share * colls.Processes[p].CrossSectionAt(energy)
	}
	return
}

func (colls *Collisions) MinThreshold() float64 {
	var result = math.MaxFloat64
	for i := range colls.Processes {
		if colls.Processes[i].Threshold != 0 && colls.Processes[i].Threshold < result {
			result = colls.Processes[i].Threshold
		}
	}
	return result
}

func (colls *Collisions) MinThresholdOfKind(t CollisionType) float64 {
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
func (colls *Collisions) SurplusCrossSection() float64 {
	if colls.MaxTCS != 0 {
		return colls.MaxTCS
	}
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

func (colls *Collisions) MaxTotalCrossSectionOverRange(e1, e2 float64) float64 {
	if len(colls.TotalCrossSectionAtCache) == 0 {
		return colls.SurplusCrossSection()
	}
	eIndex1, eIndex2 := int(math.Floor(e1/colls.EnergyStep)), min(int(math.Ceil(e2/colls.EnergyStep)), len(colls.TotalCrossSectionAtCache)-1)
	return colls.MaxTCSRangeTree.Query(eIndex1, eIndex2)
}

func (colls *Collisions) MakeEnergyGrid(minStep, maxEnergy float64) []float64 {
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

func (colls *Collisions) SampleScatteringAngleCos(energy, transitionEnergy float64, collisionType CollisionType, z AtomicNumber, species string) (cosChi float64) {
	if collisionType == ELASTIC {
		switch colls.elasticScatteringMode {
		case Born:
			panic("lxgata does not support Born elastic scattering\n")
		case Coulomb:
			return CoulombScatteringAngleSample(energy, colls.Species[species].UParameter, transitionEnergy, z)
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
			return CoulombScatteringAngleSample(energy, colls.Species[species].UParameter, transitionEnergy, z)
		case Isotropic:
			return 1. - 2.*rand.Float64()
		default:
			panic(fmt.Sprintf("unexpected lxgata.ScatteringMode: %#v", colls.inelasticScatteringMode))
		}
	}
}

func (colls *Collisions) GetTypes() (types []CollisionType) {
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
