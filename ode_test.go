package bvp

import (
	"github.com/skelterjohn/go.matrix"
	"math"
	"testing"
)

var (
	MattheijDfdbeta = []func(matrix.Matrix, float64, matrix.Matrix) matrix.Matrix{MattheijDfdbeta1, MattheijDfdbeta2}

	MattheijODE = NewODE(MattheijF, MattheijDfdx, MattheijDfdbeta, 3)
)

func MattheijDfdbeta1(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix {
	dadbeta1 := matrix.Zeros(3, 3)

	dadbeta1.Set(0, 0, -math.Cos(beta.Get(1, 0)*t))
	dadbeta1.Set(2, 0, -math.Sin(beta.Get(1, 0)*t))
	dadbeta1.Set(1, 1, 1)
	dadbeta1.Set(0, 2, math.Sin(beta.Get(1, 0)*t))
	dadbeta1.Set(2, 2, math.Cos(beta.Get(1, 0)*t))

	return matrix.Product(dadbeta1, x)
}

func MattheijDfdbeta2(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix {
	dadbeta2 := matrix.Zeros(3, 3)

	dadbeta2.Set(0, 0, t*beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))
	dadbeta2.Set(2, 0, -t*beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	dadbeta2.Set(0, 2, t*beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	dadbeta2.Set(2, 2, -t*beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))

	return matrix.Product(dadbeta2, x)
}

func TestDimensionError(t *testing.T) {

	testError := NewDimensionError("x", 4, 6, 5, 7)

	if testError.Error() != "Variable x, received dimensions (5, 7), expected dimensions (4, 6)" {
		t.Errorf("DimensionError not returning expected string representation")
	}
}

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
