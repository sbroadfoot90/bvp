package bvp

import (
	"github.com/skelterjohn/go.matrix"
	"testing"
)

var (
	MattheijDfdbeta = make([]func(matrix.Matrix, float64, matrix.Matrix) matrix.Matrix, 0, 0)

	MattheijODE = NewODE(MattheijF, MattheijDfdx, MattheijDfdbeta, 3)
)

func TestDimensionError(t *testing.T) {

	testError := NewDimensionError("x", 4, 6, 5, 7)

	if testError.Error() != "Variable x, received dimensions (5, 7), expected dimensions (4, 6)" {
		t.Errorf("DimensionError not returning expected string representation")
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
