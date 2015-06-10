package remez

import (
	"fmt"
	"math"
	"testing"
)

func almostEqual(v1, v2, e float64) bool {
	return math.Abs(v1-v2) < e
}

func TestRemezBandPass(t *testing.T) {
	bands := []float64{
		0.0, 0.1,
		0.2, 0.4,
		0.45, 0.5,
	}
	des := []float64{
		0.0, 1.0, 0.0,
	}
	weight := []float64{
		1.0, 1.0, 1.0,
	}
	gridDensity := 16
	maxIterations := 24
	h, err := Remez(13, bands, des, weight, BandPass, maxIterations, gridDensity)
	if err != nil {
		t.Fatal(err)
	}
	// if err != nil {
	// 	t.Fatal(err)
	// }
	expected := []float64{
		0.06124956814984416,
		0.1257619102752248,
		-0.06533641166341043,
		0.04658698524593622,
		-0.28615809726680885,
		-0.09231464024316056,
		0.5804898815607504,
		-0.09231464024316056,
		-0.28615809726680885,
		0.04658698524593622,
		-0.06533641166341043,
		0.1257619102752248,
		0.06124956814984416,
	}
	for i, exp := range expected {
		if !almostEqual(h[i], exp, 0.0000000001) {
			t.Errorf("Remez h[%d]=%f, want %f", i, h[i], exp)
		}
	}
}

func TestRemezDifferentiator(t *testing.T) {
	bands := []float64{
		0.0, 0.1,
		0.2, 0.4,
		0.45, 0.5,
	}
	des := []float64{
		0.0, 1.0, 0.0,
	}
	weight := []float64{
		1.0, 1.0, 1.0,
	}
	gridDensity := 16
	maxIterations := 24
	h, err := Remez(13, bands, des, weight, Differentiator, maxIterations, gridDensity)
	if err != nil {
		t.Fatal(err)
	}
	// if err != nil {
	// 	t.Fatal(err)
	// }
	expected := []float64{
		0.04114286, 0.01801784, 0.01148911, 0.00607054, -0.08632247,
		0.11984129, 0.0, -0.11984129, 0.08632247, -0.00607054,
		-0.01148911, -0.01801784, -0.04114286,
	}
	for i, exp := range expected {
		if !almostEqual(h[i], exp, 0.000001) {
			t.Errorf("Remez h[%d]=%f, want %f", i, h[i], exp)
		}
	}
}

func TestRemezHilbert(t *testing.T) {
	bands := []float64{
		0.0, 0.1,
		0.2, 0.4,
		0.45, 0.5,
	}
	des := []float64{
		0.0, 1.0, 0.0,
	}
	weight := []float64{
		1.0, 1.0, 1.0,
	}
	gridDensity := 16
	maxIterations := 24
	h, err := Remez(13, bands, des, weight, Hilbert, maxIterations, gridDensity)
	if err != nil {
		t.Fatal(err)
	}
	// if err != nil {
	// 	t.Fatal(err)
	// }
	expected := []float64{
		0.12346521, -0.00405045, -0.0227233, -0.08463579, -0.17251609,
		0.44315734, 0.0, -0.44315734, 0.17251609, 0.08463579,
		0.0227233, 0.00405045, -0.12346521,
	}
	for i, exp := range expected {
		if !almostEqual(h[i], exp, 0.000001) {
			t.Errorf("Remez h[%d]=%f, want %f", i, h[i], exp)
		}
	}
}

func TestRemez2(t *testing.T) {
	sr := 2000.0
	bands := []float64{
		0.0, 400.0 / sr,
		500.0 / sr, 0.5,
	}
	des := []float64{
		1.0, 0.0,
	}
	weight := []float64{
		1.0, 1.0,
	}
	filterType := BandPass
	gridDensity := 16
	maxIterations := 50
	h, err := Remez(21, bands, des, weight, filterType, gridDensity, maxIterations)
	if err != nil {
		t.Fatal(err)
	}
	for i, c := range h {
		fmt.Printf("%3d %f\n", i, c)
	}
}
