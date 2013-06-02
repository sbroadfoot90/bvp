package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
	"testing"
)

func TestLorenzSurfPlotB(t *testing.T) {
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
	// b := matrix.MakeDenseMatrix([]float64{-0.7924167886319804, -34.40243412416355, 2.1136}, 3, 1)
	//b := matrix.Sum(matrix.Product(B0, initialGuess[0]), matrix.Product(B1, initialGuess[n-1]))

	LorenzBVP, err := NewBVPWithInitialGuess(LorenzODE, initialGuess, timeMesh, B0, B1, beta, b)

	if err != nil {
		t.Errorf("Error creating Mattheij BVP")
	}

	err = (&LorenzBVP).SolveIVP(x0)

	if err != nil {
		t.Errorf("Error solving")
	}

	err = (&LorenzBVP).Solve()

	if err != nil {
		t.Errorf("Error solving")
	}

	// O := matrix.MakeDenseMatrix([]float64{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3, 3)
	// b1, b2, costMatrix := Surfplotb(LorenzBVP, O, 2, 0, 16.724831499260988, 16.924831499260988, 1, -40.23868387826423, -40.03868387826423)

	// WriteMatrix(b1, "b1hr.csv")
	// WriteMatrix(b2, "b2hr.csv")
	// WriteMatrix(costMatrix, "costMatrixhr.csv")
}
