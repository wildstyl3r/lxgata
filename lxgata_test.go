package lxgata

import (
	"math"
	"testing"
)

func TestLoadCrossSections(t *testing.T) {
	cs, err := LoadCrossSections("LXCat_format_test.txt")
	if err != nil {
		t.Fatalf("Unable to load cross sections %v", err.Error())
	}
	if len(cs) != 4 {
		t.Fatalf("Not all processes were loaded")
	}
	for _, p := range cs {
		if (p.Type == IONIZATION || p.Type == EXCITATION || p.Type == ROTATION) && p.Data[0].Value != 0. {
			t.Fatalf("Zero beyond threshold energy is not added for process %v", p)
		}
	}

	//todo: add more tests
}

func TestCrossSectionAt(t *testing.T) {
	collision := Collision{Data: []CrossSectionPoint{
		{1.0e+1, 0},
		{1.6e+1, 1},
		{3.2e+1, 2},
		{5.0e+2, 7},
	}}
	eps := 1e-5
	tests := map[string]struct {
		input  float64
		result float64
	}{
		"energy < any in Collision.Data": {
			input:  1.0e+0,
			result: 0,
		},
		"energy exactly at lowest point in Collision.Data": {
			input:  1.0e+1,
			result: 0,
		},
		"energy between two Collision.Data points": {
			input:  266.,
			result: 4.5,
		},
		"energy exactly at maximum in Collision.Data": {
			input:  5.0e+2,
			result: 7,
		},
		"energy > any in Collision.Data": {
			input:  6.0e+2,
			result: 7,
		},
	}

	for name, test := range tests {
		test := test // NOTE: uncomment for Go < 1.22, see /doc/faq#closures_and_goroutines
		t.Run(name, func(t *testing.T) {
			t.Parallel()
			if got, expected := collision.CrossSectionAt(test.input), test.result; math.Abs(got-expected) > eps {
				t.Fatalf("%v (c).CrossSectionAt(%v) returned %v; expected %v", name, test.input, got, expected)
			}
		})
	}
}
