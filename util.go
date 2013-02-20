package bvp

import (
	"encoding/csv"
	"github.com/skelterjohn/go.matrix"
	"os"
	"strconv"
)

func WriteMatrix(mat matrix.MatrixRO, filename string) (err error) {
	file, err := os.Create(filename)
	defer file.Close()

	if err != nil {
		return
	}

	csvWriter := csv.NewWriter(file)

	strMatrix := make([][]string, mat.Rows())

	for rowIndex := range strMatrix {
		strMatrix[rowIndex] = make([]string, mat.Cols())
		for colIndex := range strMatrix[rowIndex] {
			strMatrix[rowIndex][colIndex] = strconv.FormatFloat(mat.Get(rowIndex, colIndex), 'f', -1, 64)
		}
	}
	err = csvWriter.WriteAll(strMatrix)
	csvWriter.Flush()
	return
}

func LinearFDfdx(
	a, q func(float64, matrix.MatrixRO) matrix.Matrix,
) (
	f, dfdx func(matrix.MatrixRO, float64, matrix.MatrixRO) matrix.Matrix,
) {

	f = func(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
		return matrix.Sum(matrix.Product(a(t, beta), x), q(t, beta))
	}

	dfdx = func(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
		return a(t, beta)
	}

	return
}
