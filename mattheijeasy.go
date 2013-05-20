package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
	"math"
)

func MattheijEasyQ(t float64, beta matrix.MatrixRO) matrix.Matrix {
	q := matrix.Zeros(3, 1)
	q.Set(0, 0, math.Exp(t)*(-1+2*(math.Cos(2*t)-math.Sin(2*t))))
	q.Set(1, 0, -1*math.Exp(t))
	q.Set(2, 0, math.Exp(t)*(1-2*(math.Cos(2*t)+math.Sin(2*t))))

	return q
}

var MattheijEasyF, MattheijEasyDfdx = LinearFDfdx(MattheijA, MattheijEasyQ)

var (
	MattheijEasyODE = NewODE(MattheijEasyF, MattheijEasyDfdx, MattheijDfdbeta, 3)
)
