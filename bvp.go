package bvp

import (
	"github.com/skelterjohn/go.matrix"
)

// Boundary conditions on the ODE are B0 x(t_1) + B1 x(t_n) = b
type BVP struct {
	ODE    ODE
	X      []matrix.Matrix
	T      []float64
	B0, B1 matrix.Matrix
	b matrix.Matrix
}

func NewBVPWithInitialGuess(ode ODE, initialGuess []matrix.Matrix, timeMesh []float64, B0, B1, b matrix.Matrix) (BVP, error) {
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
	
	if b.Rows() != ode.P || b.Cols() != 1 {
		return bvp, NewDimensionError("B1", ode.P, 1, b.Rows(), b.Cols())
	}
	
	bvp = BVP{ode, initialGuess, timeMesh, B0, B1, b}
	return bvp, nil
}

func NewBVPWithoutInitialGuess(ode ODE, timeMesh []float64, B0, B1, b matrix.Matrix) (BVP, error) {
	n := len(timeMesh)
	initialGuess := make([]matrix.Matrix, n, n)
	for i := 0; i < n; i++ {
		initialGuess[i] = matrix.Zeros(ode.P, 1)
	}

	return NewBVPWithInitialGuess(ode, initialGuess, timeMesh, B0, B1, b)
}

func (bvp *BVP) solve(bta, b matrix.Matrix) {

}

// 
// properties(GetAccess = 'private', SetAccess = 'private')
//     Hf; %dx/dt = f(x, t, beta)
//     Hdf; %df/dx
//     Hdfdbta; %df/dbeta
//     f2;
//     
// end
// 
// properties(GetAccess = 'public', SetAccess = 'public')
//     p; %number of variables
//     
//     x; %solution
//     t; %time mesh
//     n; %points in mesh
//     dt; %t(i+1) - t(i)
//     
//     bta;
//     
//     q; %number of parameters
//     
//     %boundary conditions
//     B0; %B0x(ti) + B1x(tf) = b
//     B1;
//     b;
//     
//     f2defined;
// end
