package bvp

import (
	"github.com/skelterjohn/go.matrix"
	"math"
)

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

var (
	MattheijDfdbeta = []func(matrix.Matrix, float64, matrix.Matrix) matrix.Matrix{MattheijDfdbeta1, MattheijDfdbeta2}

	MattheijODE = NewODE(MattheijF, MattheijDfdx, MattheijDfdbeta, 3)
)

func MattheijDfdbeta1(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix {
	dadbeta1 := matrix.Zeros(3, 3)

	dadbeta1.Set(0, 0, -math.Cos(beta.Get(1, 0)*t))
	dadbeta1.Set(2, 0, -math.Sin(beta.Get(1, 0)*t))
	dadbeta1.Set(1, 1, 1)
	dadbeta1.Set(0, 2, math.Sin(beta.Get(1, 0)*t))
	dadbeta1.Set(2, 2, math.Cos(beta.Get(1, 0)*t))

	return matrix.Product(dadbeta1, x)
}

func MattheijDfdbeta2(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix {
	dadbeta2 := matrix.Zeros(3, 3)

	dadbeta2.Set(0, 0, t*beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))
	dadbeta2.Set(2, 0, -t*beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	dadbeta2.Set(0, 2, t*beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	dadbeta2.Set(2, 2, -t*beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))

	return matrix.Product(dadbeta2, x)
}
