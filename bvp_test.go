package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
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

func TestConstraintVectorBlocks(t *testing.T) {
	n := 101

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

	_, err = ConstraintVectorBlocks(&MattheijBVP)
}

func TestConstraintMatrixBlocks(t *testing.T) {
	n := 11

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

	A, B, err := ConstraintMatrixBlocks(&MattheijBVP)

	if err != nil {
		t.Errorf("Error calculating constraint matrix")
	}
	for i := 0; i < (n - 1); i++ {
		if A[i].Rows() != 3 {
			t.Errorf("Incorrect number of rows for A")
		}
		if A[i].Cols() != 3 {
			t.Errorf("Incorrect number of cols for A")
		}
		if B[i].Rows() != 3 {
			t.Errorf("Incorrect number of rows for B")
		}
		if B[i].Cols() != 3 {
			t.Errorf("Incorrect number of cols for B")
		}
	}
}

func TestGetDelta(t *testing.T) {
	n := 11

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

	_, err = getDelta(&MattheijBVP)

	if err != nil {
		t.Errorf("Error calculating delta")
	}
	// WriteMatrices(delta, "delta.csv")
}

func TestSolveMattheij(t *testing.T) {
	n := 1001

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

	err = (&MattheijBVP).Solve()

	if err != nil {
		t.Errorf("Error solving")
	}
	// WriteMatrices(MattheijBVP.X, "x.csv")
}

func TestMattheijBCs(t *testing.T) {
	n := 1001

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)

	for i := 0; i < n; i++ {
		timeMesh[i] = float64(i) / (float64(n) - 1)
		initialGuess[i] = matrix.Scaled(matrix.Ones(3, 1), math.Exp(timeMesh[i]))
	}

	B0 := matrix.MakeDenseMatrix([]float64{1, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3)
	B1 := matrix.MakeDenseMatrix([]float64{0, 0, 0, 0, 1, 0, 0.8415, 0, 0.5403}, 3, 3)
	beta := matrix.MakeDenseMatrix([]float64{2, 2}, 2, 1)
	b := matrix.Sum(matrix.Product(B0, initialGuess[0]), matrix.Product(B1, initialGuess[n-1]))

	MattheijBVP, err := NewBVPWithInitialGuess(MattheijODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Error creating Mattheij BVP")
	}

	SetOptimalBoundaryMatrices(&MattheijBVP)

	// WriteMatrix(MattheijBVP.B0, "B0.csv")
	// WriteMatrix(MattheijBVP.B1, "B1.csv")
	// WriteMatrix(MattheijBVP.B, "B.csv")
}

func TestSolveLorenz(t *testing.T) {
	n := 3001

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)
	tf := 3.0

	for i := 0; i < n; i++ {
		timeMesh[i] = float64(i) * tf / (float64(n) - 1)
		initialGuess[i] = matrix.Scaled(matrix.Ones(3, 1), math.Exp(timeMesh[i]))
	}

	B0 := matrix.MakeDenseMatrix([]float64{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, 3)
	B1 := matrix.MakeDenseMatrix([]float64{0, 0, 0, 0, 0, 0, 0, 0, 0}, 3, 3)
	beta := matrix.MakeDenseMatrix([]float64{10, 28, 8. / 3.}, 3, 1)
	b := matrix.MakeDenseMatrix([]float64{1, 1, 30}, 3, 1)
	//b := matrix.Sum(matrix.Product(B0, initialGuess[0]), matrix.Product(B1, initialGuess[n-1]))

	LorenzBVP, err := NewBVPWithInitialGuess(LorenzODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Error creating Mattheij BVP")
	}

	err = (&LorenzBVP).SolveIVP(b)

	if err != nil {
		t.Errorf("Error solving")
	}
	// WriteMatrices(LorenzBVP.X, "xlorenz.csv")
	err = (&LorenzBVP).Solve()

	if err != nil {
		t.Errorf("Error solving")
	}

	// WriteMatrices(LorenzBVP.X, "xlorenz.csv")
}

func TestSolveLorenzBVP(t *testing.T) {
	n := 3001

	timeMesh := make([]float64, n, n)
	initialGuess := make([]matrix.Matrix, n, n)
	tf := 3.0

	for i := 0; i < n; i++ {
		timeMesh[i] = float64(i) * tf / (float64(n) - 1)
		initialGuess[i] = matrix.Ones(3, 1)
	}

	B0 := matrix.MakeDenseMatrix([]float64{-0.0155, 0.0084, -0.2942, 0.0483, 0.0061, -0.9790, 0.1958, 0.1958, 0.0574}, 3, 3)
	B1 := matrix.MakeDenseMatrix([]float64{-0.4504, -0.7696, 0.3434, 0, -0.4694, -0.3611, 0, 0, 0}, 3, 3)
	beta := matrix.MakeDenseMatrix([]float64{10, 28, 8. / 3.}, 3, 1)
	x0 := matrix.MakeDenseMatrix([]float64{1, 1, 30}, 3, 1)
	b := matrix.MakeDenseMatrix([]float64{16.824831499260988, -40.13868387826423, 2.1136}, 3, 1)

	LorenzBVP, err := NewBVPWithInitialGuess(LorenzODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Error creating Mattheij BVP")
	}

	err = (&LorenzBVP).SolveIVP(x0)

	if err != nil {
		t.Errorf("Error solving")
	}
	// WriteMatrices(LorenzBVP.X, "xinitlorenz.csv")

	err = (&LorenzBVP).Solve()

	if err != nil {
		t.Errorf("Error solving")
	}

	// WriteMatrices(LorenzBVP.X, "xlorenz.csv")
}
