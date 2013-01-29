package bvp

import (
	"github.com/skelterjohn/go.matrix"
	"math"
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

func MattheijA(t float64, beta matrix.Matrix) matrix.Matrix {
	a := matrix.Zeros(3, 3)

	a.Set(0, 0, 1-beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	a.Set(2, 0, 1-beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))
	a.Set(1, 1, beta.Get(0, 0))
	a.Set(0, 2, -1+beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))
	a.Set(2, 2, 1+beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	return a
}

func MattheijQ(t float64, beta matrix.Matrix) matrix.Matrix {
	q := matrix.Zeros(3, 1)
	q.Set(0, 0, math.Exp(t)*(-1+19*(math.Cos(2*t)-math.Sin(2*t))))
	q.Set(2, 0, math.Exp(t)*(1-19*(math.Cos(2*t)+math.Sin(2*t))))

	return q
}

var MattheijF, MattheijDfdx = LinearFDfdx(MattheijA, MattheijQ)

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
