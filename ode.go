/*
    Package bvp implements a simple library to solve boundary value problems.
	The package makes use of the github.com/skelterjohn/go.matrix library.
*/
package bvp

import (
	"fmt"
	"github.com/skelterjohn/go.matrix"
)

type ODE struct {
	// dx/dt = F(x, t, beta)
	f, dfdx func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix
	dfdbeta []func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix

	P, // number of varibales (length of x)
	Q int // number of parameters (length of beta)
}

func NewODE(f, dfdx func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix,
	dfdbeta []func(x matrix.Matrix, t float64, beta matrix.Matrix) matrix.Matrix,
	p int) ODE {
	return ODE{f, dfdx, dfdbeta, p, len(dfdbeta)}
}

// Evaluates the function f with checking of matrix dimensions
func (o *ODE) F(x matrix.Matrix, t float64, beta matrix.Matrix) (matrix.Matrix, error) {
	//checking dimensions of input
	if x.Rows() != o.P || x.Cols() != 1 {
		return nil, NewDimensionError("x", o.P, 1, x.Rows(), x.Cols())
	}

	if o.Q != 0 && (beta.Rows() != o.Q || beta.Cols() != 1) {
		return nil, NewDimensionError("beta", o.Q, 1, beta.Rows(), beta.Cols())
	}

	return o.f(x, t, beta), nil
}

type DimensionError struct {
	VariableName                                       string
	ExpectedRow, ExpectedCol, ReceivedRow, ReceivedCol int
}

func (de DimensionError) Error() string {
	return fmt.Sprintf("Variable %s, received dimensions (%d, %d), expected dimensions (%d, %d)", de.VariableName, de.ReceivedRow, de.ReceivedCol, de.ExpectedRow, de.ExpectedCol)
}

func NewDimensionError(variableName string, er, ec, rr, rc int) DimensionError {
	return DimensionError{variableName, er, ec, rr, rc}
}
