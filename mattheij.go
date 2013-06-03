package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
	"math"
)

func MattheijA(t float64, beta matrix.MatrixRO) matrix.Matrix {
	a := matrix.Zeros(3, 3)

	a.Set(0, 0, 1-beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	a.Set(0, 2, 1+beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))
	a.Set(1, 1, beta.Get(0, 0))
	a.Set(2, 0, -1+beta.Get(0, 0)*math.Sin(beta.Get(1, 0)*t))
	a.Set(2, 2, 1+beta.Get(0, 0)*math.Cos(beta.Get(1, 0)*t))
	return a
}

func MattheijQ(t float64, beta matrix.MatrixRO) matrix.Matrix {
	q := matrix.Zeros(3, 1)
	q.Set(0, 0, math.Exp(t)*(-1+19*(math.Cos(2*t)-math.Sin(2*t))))
	q.Set(1, 0, -18*math.Exp(t))
	q.Set(2, 0, math.Exp(t)*(1-19*(math.Cos(2*t)+math.Sin(2*t))))

	return q
}

var MattheijF, MattheijDfdx = LinearFDfdx(MattheijA, MattheijQ)

var MattheijODE = NewODE(MattheijF, MattheijDfdx, 3, 2)