package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
	"os"
	"testing"
)

func TestWriteMatrix(t *testing.T) {

	a, err := matrix.ParseMatlab("[1 2 3; 4 2 1]")

	if err != nil {
		t.Errorf(err.Error())
		return
	}
	err = WriteMatrix(a, "temp.csv")
	defer os.Remove("temp.csv")

	if err != nil {
		t.Errorf(err.Error())
	}

	tempCsv, err := os.Open("temp.csv")

	defer tempCsv.Close()

	if err != nil {
		t.Errorf(err.Error())
	}
	expectedContents := []byte("1,2,3\n4,2,1\n")

	tempCsvContents := make([]byte, 12, 12)

	n, err := tempCsv.Read(tempCsvContents)
	if n != len(expectedContents) {
		t.Error("Length of temp.csv is incorrect, expected: %d, got %d", len(expectedContents), n)
	}

	for i := range tempCsvContents {
		if expectedContents[i] != tempCsvContents[i] {
			t.Error("Character %d is incorrect, expected: %c, got %c", i, expectedContents[i], tempCsvContents[i])
		}
	}

	return
}

func TestLinearFDfdxDfdbeta(t *testing.T) {
	x := matrix.Zeros(3, 1)
	x.Set(0, 0, 1)
	x.Set(1, 0, 0.1)
	x.Set(2, 0, 0.14)
	time := 0.5
	beta := matrix.Zeros(2, 1)
	beta.Set(0, 0, 19)
	beta.Set(1, 0, 2)

	if !matrix.ApproxEquals(MattheijF(x, time, beta), matrix.Sum(matrix.Product(MattheijA(time, beta), x), MattheijQ(time, beta)), 10e-8) {
		t.Error("f function in linearFDfdxDfdbeta is broken")
	}

	if !matrix.ApproxEquals(MattheijDfdx(x, time, beta), MattheijA(time, beta), 10e-8) {
		t.Error("dfdx function in linearFDfdxDfdbeta is broken")
	}
}

// 
// func TestSparseQR(t *testing.T) {
// 	n := 100
// 	A := matrix.NormalsSparse(n, n, 10*n)
// 
// 	Q, R := SparseQR(A)
// 
// 	QR, err := Q.TimesSparse(R)
// 
// 	if err != nil {
// 		t.Error("Dimension mismatch in Q and R")
// 	}
// 
// 	if !matrix.ApproxEquals(A, QR, 10e-8) {
// 		t.Error("Q*R does not equal A")
// 	}
// 
// 	if !matrix.ApproxEquals(R, R.U(), 10e-8) {
// 		t.Error("R is not upper triangular")
// 	}
// 
// 	QTQ, err := Q.TimesSparse(Q.Transpose())
// 
// 	if err != nil {
// 		t.Error("Dimension mismatch in Q and R")
// 	}
// 
// 	if !matrix.ApproxEquals(QTQ, sparseEye(n), 10e-8) {
// 		t.Error("Q is not orthogonal")
// 	}
// }
// 
// func TestSolveUpper(t *testing.T) {
// 	n := 25
// 	A := matrix.NormalsSparse(n, n, n*n*n).U()
// 	b := matrix.NormalsSparse(n, 1, n)
// 	x := solveUpper(A, b)
// 
// 	Ax, err := A.Times(x)
// 
// 	if err != nil {
// 		t.Error("Dimension mismatch in A and x")
// 	}
// 
// 	if !matrix.ApproxEquals(Ax, b, 10e-4) {
// 		t.Error("A*x does not equal b")
// 	}
// 
// }
// 
// 
// 
// func TestSolve(t *testing.T) {
// 	n := 205
// 	A := matrix.NormalsSparse(n, n, n*n*n)
// 	b := matrix.NormalsSparse(n, 1, n)
// 	x := solve(A, b)
// 
// 	Ax, err := A.Times(x)
// 
// 	if err != nil {
// 		t.Error("Dimension mismatch in A and x")
// 	}
// 
// 	if !matrix.ApproxEquals(Ax, b, 10e-4) {
// 		t.Error("A*x does not equal b")
// 	}
// 
// }
