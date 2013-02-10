package bvp

import (
	"github.com/skelterjohn/go.matrix"
)

// Boundary conditions on the ODE are B0 x(t_1) + B1 x(t_n) = b
type BVP struct {
	ODE     ODE
	X       []matrix.Matrix
	T       []float64
	B0, B1  matrix.Matrix
	Beta, B matrix.Matrix
}

func NewBVPWithInitialGuess(ode ODE, initialGuess []matrix.Matrix, timeMesh []float64, B0, B1, beta, b matrix.Matrix) (BVP, error) {
	var bvp BVP

	n := len(timeMesh)

	if len(initialGuess) != n {
		return bvp, NewDimensionError("initialGuess", ode.P, n, ode.P, len(initialGuess))
	}

	for i := 0; i < n; i++ {
		if initialGuess[i].Rows() != ode.P || initialGuess[i].Cols() != 1 {
			return bvp, NewDimensionError("initialGuess", ode.P, 1, initialGuess[i].Rows(), initialGuess[i].Cols())
		}
	}

	if B0.Rows() != ode.P || B0.Cols() != ode.P {
		return bvp, NewDimensionError("B0", ode.P, ode.P, B0.Rows(), B0.Cols())
	}

	if B1.Rows() != ode.P || B1.Cols() != ode.P {
		return bvp, NewDimensionError("B1", ode.P, ode.P, B1.Rows(), B1.Cols())
	}

	if beta.Rows() != ode.Q || beta.Cols() != 1 {
		return bvp, NewDimensionError("beta", ode.Q, 1, beta.Rows(), beta.Cols())
	}
	if b.Rows() != ode.P || b.Cols() != 1 {
		return bvp, NewDimensionError("b", ode.P, 1, b.Rows(), b.Cols())
	}

	bvp = BVP{ode, initialGuess, timeMesh, B0, B1, beta, b}
	return bvp, nil
}

func NewBVPWithoutInitialGuess(ode ODE, timeMesh []float64, B0, B1, beta, b matrix.Matrix) (BVP, error) {
	n := len(timeMesh)
	initialGuess := make([]matrix.Matrix, n, n)
	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(ode.P, 1)
	}

	return NewBVPWithInitialGuess(ode, initialGuess, timeMesh, B0, B1, beta, b)
}

func (bvp *BVP) Solve() {

}
