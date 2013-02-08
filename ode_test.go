package bvp

import (
	"github.com/skelterjohn/go.matrix"
	"testing"
)

func TestCalculateQOnCreateODE(t *testing.T) {
	if MattheijODE.Q != 2 {
		t.Errorf("Detection of q was incorrect, expected %d, got %d", 2, MattheijODE.Q)
	}
}

func TestChecksXDim(t *testing.T) {
	ODEf, err := MattheijODE.F(matrix.Zeros(3, 1), 0, matrix.Zeros(2, 1))

	if err != nil || ODEf == nil {
		t.Errorf("Incorrect x dim checking. Correct dims used.")
	}

	ODEf, err = MattheijODE.F(matrix.Zeros(4, 1), 0, matrix.Zeros(2, 1))
	if err == nil {
		t.Errorf("Incorrect x dim checking. Incorrect rows not detected.")
	}

	ODEf, err = MattheijODE.F(matrix.Zeros(3, 7), 0, matrix.Zeros(2, 1))
	if err == nil {
		t.Errorf("Incorrect x dim checking. Incorrect cols not detected.")
	}
}

func TestChecksBetaDim(t *testing.T) {
	ODEf, err := MattheijODE.F(matrix.Zeros(3, 1), 0, matrix.Zeros(2, 1))

	if err != nil || ODEf == nil {
		t.Errorf("Incorrect beta dim checking. Correct dims used.")
	}

	ODEf, err = MattheijODE.F(matrix.Zeros(3, 1), 0, matrix.Zeros(4, 1))
	if err == nil {
		t.Errorf("Incorrect beta dim checking. Incorrect rows not detected.")
	}

	ODEf, err = MattheijODE.F(matrix.Zeros(3, 1), 0, matrix.Zeros(2, 4))
	if err == nil {
		t.Errorf("Incorrect beta dim checking. Incorrect cols not detected.")
	}
}
