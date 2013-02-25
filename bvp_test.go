package bvp

import (
	"github.com/skelterjohn/go.matrix"
	"math"
	"testing"
)

func TestChecksInitialGuessLength(t *testing.T) {
	n := 1000

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(3, 1)
	}

	B0 := matrix.Zeros(3, 3)
	B1 := matrix.Zeros(3, 3)
	beta := matrix.Zeros(2, 1)
	b := matrix.Zeros(3, 1)

	_, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Incorrect initial guess length check. Correct dims used.")
	}

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess[1:100], timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect initial guess length check. Incorrect dims used.")
	}
}

func TestChecksInitialGuessDims(t *testing.T) {
	n := 1000

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(3, 1)
	}

	B0 := matrix.Zeros(3, 3)
	B1 := matrix.Zeros(3, 3)
	beta := matrix.Zeros(2, 1)
	b := matrix.Zeros(3, 1)

	_, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Incorrect initial guess length check. Correct dims used.")
	}

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(4, 1)
	}

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect initial guess length check. Incorrect rows used.")
	}

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(3, 2)
	}

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect initial guess length check. Incorrect rows used.")
	}
}

func TestChecksB0B1Dims(t *testing.T) {
	n := 1000

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(3, 1)
	}

	B0 := matrix.Zeros(3, 3)
	B1 := matrix.Zeros(3, 3)
	beta := matrix.Zeros(2, 1)
	b := matrix.Zeros(3, 1)

	_, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Incorrect boundary matrices check. Correct dims used.")
	}

	B0 = matrix.Zeros(2, 3)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect boundary matrices check. Incorrect B0 rows used.")
	}

	B0 = matrix.Zeros(3, 2)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect boundary matrices check. Incorrect B0 cols used.")
	}

	B1 = matrix.Zeros(2, 3)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect boundary matrices check. Incorrect B1 rows used.")
	}

	B1 = matrix.Zeros(3, 2)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect boundary matrices check. Incorrect B1 cols used.")
	}
}

func TestChecksBetaDims(t *testing.T) {
	n := 1000

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(3, 1)
	}

	B0 := matrix.Zeros(3, 3)
	B1 := matrix.Zeros(3, 3)
	beta := matrix.Zeros(2, 1)
	b := matrix.Zeros(3, 1)

	_, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Incorrect b check. Correct dims used.")
	}

	beta = matrix.Zeros(3, 1)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect b check. Incorrect b rows used.")
	}

	beta = matrix.Zeros(3, 1)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect b check. Incorrect b cols used.")
	}
}

func TestChecksBDims(t *testing.T) {
	n := 1000

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(3, 1)
	}

	B0 := matrix.Zeros(3, 3)
	B1 := matrix.Zeros(3, 3)
	beta := matrix.Zeros(2, 1)
	b := matrix.Zeros(3, 1)

	_, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Incorrect b check. Correct dims used.")
	}

	b = matrix.Zeros(2, 1)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect b check. Incorrect b rows used.")
	}

	b = matrix.Zeros(3, 2)

	_, err = NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err == nil {
		t.Errorf("Incorrect b check. Incorrect b cols used.")
	}
}

func TestConstraintVector(t *testing.T) {
	n := 10001

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		timeMesh[i] = float64(i) / (float64(n) - 1)
		initialGuess[i] = matrix.Scaled(matrix.Ones(3, 1), math.Exp(timeMesh[i]))
	}

	B0 := matrix.MakeDenseMatrix([]float64{1, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3)
	B1 := matrix.MakeDenseMatrix([]float64{0, 0, 0, 0, 1, 0, 0.8415, 0, 0.5403}, 3, 3)
	beta := matrix.MakeDenseMatrix([]float64{19, 2}, 2, 1)
	b := matrix.Sum(matrix.Product(B0, initialGuess[0]), matrix.Product(B1, initialGuess[n-1]))

	MattheijBVP, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Error creating Mattheij BVP")
	}

	sol, err := ConstraintVector(&MattheijBVP)

	if err != nil {
		t.Errorf("Error calculating constraint vector")
	}

	if !matrix.ApproxEquals(sol, matrix.Zeros(n*3, 1), 1e-2) {
		t.Errorf("Constraint vector not close to zero for actual solution")
	}

	for i := 0; i < n; i++ {
		MattheijBVP.X[i] = matrix.Numbers(3, 1, 340)
	}

	sol, err = ConstraintVector(&MattheijBVP)

	if err != nil {
		t.Errorf("Error calculating constraint vector")
	}

	if matrix.ApproxEquals(sol, matrix.Zeros(n*3, 1), 1e-2) {
		t.Errorf("Constraint vector too close to zero for non-solution")
	}
}


func TestConstraintMatrix(t *testing.T) {
	n := 10001

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		timeMesh[i] = float64(i) / (float64(n) - 1)
		initialGuess[i] = matrix.Scaled(matrix.Ones(3, 1), math.Exp(timeMesh[i]))
	}

	B0 := matrix.MakeDenseMatrix([]float64{1, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3)
	B1 := matrix.MakeDenseMatrix([]float64{0, 0, 0, 0, 1, 0, 0.8415, 0, 0.5403}, 3, 3)
	beta := matrix.MakeDenseMatrix([]float64{19, 2}, 2, 1)
	b := matrix.Sum(matrix.Product(B0, initialGuess[0]), matrix.Product(B1, initialGuess[n-1]))

	MattheijBVP, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Error creating Mattheij BVP")
	}

	_, err = ConstraintMatrix(&MattheijBVP)

	if err != nil {
		t.Errorf("Error calculating constraint vector")
	}
}