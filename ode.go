/*
    Package bvp implements a simple library to solve boundary value problems.
	The package makes use of the github.com/skelterjohn/go.matrix library.
*/
package bvp

import (
	"github.com/skelterjohn/go.matrix"
)

type ODE struct {
	// dx/dt = F(x, t, beta)
	f, dfdx func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix
	dfdbeta []func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix

	P, // number of varibales (length of x)
	Q int // number of parameters (length of beta)
}

// Evaluates the function f with checking of matrix dimensions
func (o *ODE) F(x matrix.Matrix, t float64, beta matrix.Matrix) (matrix.Matrix, error) {
	return o.f(x, t, beta), nil
}
