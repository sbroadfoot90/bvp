package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
)

var (
	LorenzDfdbeta = []func(matrix.MatrixRO, float64, matrix.MatrixRO) matrix.Matrix{LorenzDfdbeta1, LorenzDfdbeta2, LorenzDfdbeta3}

	LorenzODE = NewODE(LorenzF, LorenzDfdx, LorenzDfdbeta, 3)
)

func LorenzF(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
	return matrix.MakeDenseMatrix([]float64{
		beta.Get(0, 0) * (x.Get(1, 0) - x.Get(0, 0)),
		x.Get(0, 0)*(beta.Get(1, 0)-x.Get(2, 0)) - x.Get(1, 0),
		x.Get(0, 0)*x.Get(1, 0) - beta.Get(2, 0)*x.Get(2, 0),
	}, 3, 1)
}

func LorenzDfdx(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
	return matrix.MakeDenseMatrix([]float64{
		-beta.Get(0, 0), beta.Get(0, 0), 0,
		beta.Get(1, 0) - x.Get(2, 0), -1, -x.Get(0, 0),
		x.Get(0, 0), x.Get(0, 0), -beta.Get(2, 0),
	}, 3, 3)
}

func LorenzDfdbeta1(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
	return matrix.MakeDenseMatrix([]float64{
		x.Get(1, 0) - x.Get(0, 0),
		(beta.Get(1, 0) - x.Get(2, 0)),
		0,
	}, 3, 1)
}

func LorenzDfdbeta2(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
	return matrix.MakeDenseMatrix([]float64{
		0,
		x.Get(0, 0),
		0,
	}, 3, 1)
}

func LorenzDfdbeta3(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
	return matrix.MakeDenseMatrix([]float64{
		0,
		0,
		-x.Get(2, 0),
	}, 3, 1)
}
