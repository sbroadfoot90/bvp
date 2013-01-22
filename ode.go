package bvp

import (
	"github.com/skelterjohn/go.matrix"
)

type ODE struct {
	// dx/dt = f(x, t, beta)
	F, DfDx func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix
	DfDbeta []func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix
	
	P,    // number of varibales (length of x)
	Q int // number of parameters (length of beta)
}
