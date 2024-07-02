package lxgata

import (
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
