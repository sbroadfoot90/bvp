package bvp

import (
	"github.com/skelterjohn/go.matrix"
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
