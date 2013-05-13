package bvp

import (
	"encoding/csv"
	"github.com/sbroadfoot90/go.matrix"
	"math"
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

// func SparseQR(A *matrix.SparseMatrix) (Q, R *matrix.SparseMatrix, err error) {
// 	m := A.Rows()
// 	// n := A.Cols()
// 	
// 	Q = sparseEye(m)
// 	R = A.Copy()
// 	
// 	for i := 0; i < m - 1; i++ {
// 		x := R.GetColVector(i).Copy()
// 		
// 		if i != 0 {
// 			for j := 0; j < i; j++ {
// 				x.Set(j, 0, 0)
// 			}
// 		}
// 		x.Set(i, 0, x.Get(i, 0)+sign(x.Get(i, 0))*x.TwoNorm())
// 		
// 		x.Scale(1/x.TwoNorm())
// 
// 		var xtx *matrix.SparseMatrix
// 		xtx, err = x.TimesSparse(x.Transpose())
// 		
// 		if err != nil {
// 			err = MatrixError(fmt.Sprintf("Error in outerproduct: %s", err.Error()))
// 			return
// 		}
// 		xtx.Scale(2)
// 		
// 		var Qi *matrix.SparseMatrix
// 		Qi, err = sparseEye(m).MinusSparse(xtx)
// 		
// 		if err != nil {
// 			err = MatrixError(fmt.Sprintf("Error in calculating Qi: %s", err.Error()))
// 			return
// 		}
// 
// 		Q, err = Q.TimesSparse(Qi)
// 		if err != nil {
// 			err = MatrixError(fmt.Sprintf("Error in calculating Q: %s", err.Error()))
// 			return
// 		}
// 		
// 		R, err = Qi.TimesSparse(R)
// 		if err != nil {
// 			err = MatrixError(fmt.Sprintf("Error in calculating R: %s", err.Error()))
// 			return
// 		}
// 	}
// 	
// 	return
// }

func sparseEye(span int) (I *matrix.SparseMatrix) {
	I = matrix.ZerosSparse(span, span)

	for i := 0; i < span; i++ {
		I.Set(i, i, 1)
	}

	return
}

func sign(x float64) float64 {
	if x < 0 {
		return -1
	}
	return 1
}

func SparseQR(A *matrix.SparseMatrix) (Q, R *matrix.SparseMatrix) {
	m := A.Rows()
	n := A.Cols()
	QR := A.Copy()
	Q = matrix.ZerosSparse(m, n)
	R = matrix.ZerosSparse(m, n)
	i, j, k := 0, 0, 0
	norm := float64(0.0)
	s := float64(0.0)

	for k = 0; k < n; k++ {
		norm = 0
		for i = k; i < m; i++ {
			norm = math.Hypot(norm, QR.Get(i, k))
		}

		if norm != 0.0 {
			if QR.Get(k, k) < 0 {
				norm = -norm
			}

			for i = k; i < m; i++ {
				QR.Set(i, k, QR.Get(i, k)/norm)
			}
			QR.Set(k, k, QR.Get(k, k)+1.0)

			for j = k + 1; j < n; j++ {
				s = 0.0
				for i = k; i < m; i++ {
					s += QR.Get(i, k) * QR.Get(i, j)
				}
				s = -s / QR.Get(k, k)
				for i = k; i < m; i++ {
					QR.Set(i, j, QR.Get(i, j)+s*QR.Get(i, k))

					if i < j {
						R.Set(i, j, QR.Get(i, j))
					}
				}

			}
		}

		R.Set(k, k, -norm)

	}

	//Q Matrix:
	i, j, k = 0, 0, 0

	for k = n - 1; k >= 0; k-- {
		Q.Set(k, k, 1.0)
		for j = k; j < n; j++ {
			if QR.Get(k, k) != 0 {
				s = 0.0
				for i = k; i < m; i++ {
					s += QR.Get(i, k) * Q.Get(i, j)
				}
				s = -s / QR.Get(k, k)
				for i = k; i < m; i++ {
					Q.Set(i, j, Q.Get(i, j)+s*QR.Get(i, k))
				}
			}
		}
	}

	return
}

func solveUpper(A matrix.MatrixRO, b matrix.MatrixRO) *matrix.DenseMatrix {
	x := make([]float64, A.Cols())
	for i := A.Rows() - 1; i >= 0; i-- {
		x[i] = b.Get(i, 0)
		for j := i + 1; j < A.Cols(); j++ {
			x[i] -= x[j] * A.Get(i, j)
		}
		x[i] /= A.Get(i, i)
	}
	return matrix.MakeDenseMatrix(x, A.Cols(), 1)
}

func solve(A *matrix.SparseMatrix, b matrix.MatrixRO) matrix.Matrix {
	Q, R := SparseQR(A)
	Qtb, _ := Q.Transpose().Times(b)
	
	return solveUpper(R, Qtb)
}