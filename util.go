package bvp

import (
	"encoding/csv"
	"github.com/skelterjohn/go.matrix"
	"os"
	"strconv"
)

func writeMatrix(mat matrix.MatrixRO, filename string) (err error) {
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
