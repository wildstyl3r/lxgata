package lxgata

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type ProcessType string

const ELASTIC, EFFECTIVE, EXCITATION, ATTACHMENT, IONIZATION, ROTATION ProcessType = "ELASTIC", "EFFECTIVE", "EXCITATION", "ATTACHMENT", "IONIZATION", "ROTATION"

// Cross section in [m^2] at energy [eV]
type CrossSectionPoint struct {
	Energy, Value float64
}

type Process struct {
	Type            ProcessType
	MassRatio       float64
	Species         string
	Data            []CrossSectionPoint
	Threshold       float64
	LowerEnergy     float64
	LowerStatWeight float64
	UpperEnergy     float64
	UpperStatWeight float64
	Info            map[string]string
}

func (p *Process) CrossSectionAt(energy float64) float64 {
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

func (p Process) String() string {
	return fmt.Sprintf("Cross section of %v %v. Threshold: %v", p.Species, strings.ToLower(string(p.Type)), p.Threshold)
}

func LoadCrossSections(fileName string) ([]Process, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	setProcessTypes := map[string]struct{}{string(ELASTIC): {}, string(EFFECTIVE): {}, string(EXCITATION): {}, string(ATTACHMENT): {}, string(IONIZATION): {}, string(ROTATION): {}}

	var processes []Process

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		if len(tokens) == 0 {
			continue
		}

		if _, ok := setProcessTypes[tokens[0]]; ok {
			processType := ProcessType(tokens[0])
			scanner.Scan()
			species, _, _ := strings.Cut(scanner.Text(), " ")
			scanner.Scan()
			parameters := strings.Fields(strings.Trim(scanner.Text(), " "))
			massRatio, threshold, lowerEnergy, upperEnergy, lowerStatWeight, upperStatWeight := 0., 0., 0., 0., 0., 0.
			switch processType {
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
				if !(processType == IONIZATION || processType == EXCITATION || processType == ROTATION) || threshold < energy {
					data = append(data, CrossSectionPoint{energy, crossSection})
				}
				scanner.Scan()
			}

			if (processType == IONIZATION || processType == EXCITATION || processType == ROTATION) && data[0].Value > 0. {
				data = append([]CrossSectionPoint{{threshold, 0.}}, data...)
			}

			processes = append(processes, Process{
				Type:            processType,
				MassRatio:       massRatio,
				Species:         species,
				Data:            data,
				Threshold:       threshold,
				LowerEnergy:     lowerEnergy,
				LowerStatWeight: lowerStatWeight,
				UpperEnergy:     upperEnergy,
				UpperStatWeight: upperStatWeight,
				Info:            info,
			})
		}
	}
	return processes, nil
}
