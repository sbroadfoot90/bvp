package bvp

import (
	"github.com/sbroadfoot90/go.matrix"
)

// Boundary conditions on the ODE are B0 x(t_1) + B1 x(t_n) = b
type BVP struct {
	ODE     ODE
	X       []matrix.Matrix
	T       []float64
	B0, B1  matrix.Matrix
	Beta, B matrix.Matrix
	N       int
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

	bvp = BVP{ode, initialGuess, timeMesh, B0, B1, beta, b, n}
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
	// tolerance := 1e-4
	// maxiter := 500
	// converged := false
}

func ConstraintVector(bvp *BVP) (constraint *matrix.DenseMatrix, err error) {

	constraint = matrix.Zeros(bvp.ODE.P*bvp.N, 1)

	var constraint_i *matrix.DenseMatrix
	fNow, err := bvp.ODE.F(bvp.X[0], bvp.T[0], bvp.Beta)

	if err != nil {
		return
	}

	for i := 1; i < bvp.N; i++ {

		fBefore := matrix.MakeDenseCopy(fNow)

		fNow, err = bvp.ODE.F(bvp.X[i], bvp.T[i], bvp.Beta)

		if err != nil {
			return
		}

		constraint_i = matrix.Difference(
			matrix.Difference(bvp.X[i], bvp.X[i-1]),
			matrix.Scaled(matrix.Sum(fBefore, fNow), (bvp.T[i]-bvp.T[i-1])/2),
		)

		constraint.SetMatrix((i-1)*bvp.ODE.P, 0, constraint_i)
	}

	constraint.SetMatrix((bvp.N-1)*bvp.ODE.P, 0,
		matrix.Difference(matrix.Sum(matrix.Product(bvp.B0, bvp.X[0]), matrix.Product(bvp.B1, bvp.X[bvp.N-1])), bvp.B))

	return
}

func ConstraintMatrixBlocks(bvp *BVP) (A, B []*matrix.DenseMatrix, err error) {
	A = make([]*matrix.DenseMatrix, bvp.N - 1, bvp.N - 1)
	B = make([]*matrix.DenseMatrix, bvp.N - 1, bvp.N - 1)
	
	dfdxNow, err := bvp.ODE.Dfdx(bvp.X[0], bvp.T[0], bvp.Beta)
	
	if err != nil {
		return
	}

	for i := 1; i < bvp.N; i++ {

		dfdxBefore := matrix.MakeDenseCopy(dfdxNow)

		dfdxNow, err = bvp.ODE.Dfdx(bvp.X[i], bvp.T[i], bvp.Beta)

		if err != nil {
			return
		}

		A[i-1] = matrix.Difference(
			matrix.Scaled(matrix.Eye(bvp.ODE.P), -1),
			matrix.Scaled(dfdxBefore, (bvp.T[i]-bvp.T[i-1])/2),
		)
		
		B[i-1] = matrix.Difference(
			matrix.Eye(bvp.ODE.P),
			matrix.Scaled(dfdxBefore, (bvp.T[i]-bvp.T[i-1])/2),
		)

		// for r := 0; r < bvp.ODE.P; r++ {
		// 			for c := 0; c < bvp.ODE.P; c++ {
		// 				constraint.Set(r+(i-1)*bvp.ODE.P, c+(i-1)*bvp.ODE.P, constraint_i1.Get(r, c))
		// 				constraint.Set(r+(i-1)*bvp.ODE.P, c+i*bvp.ODE.P, constraint_i2.Get(r, c))
		// 			}
		// 		}
	}

	// for r := 0; r < bvp.ODE.P; r++ {
	// 		for c := 0; c < bvp.ODE.P; c++ {
	// 			constraint.Set(r+(bvp.N-1)*bvp.ODE.P, c, bvp.B0.Get(r, c))
	// 			constraint.Set(r+(bvp.N-1)*bvp.ODE.P, c+(bvp.N-1)*bvp.ODE.P, bvp.B1.Get(r, c))
	// 		}
	// 	}

	return
}










func ConstraintMatrix(bvp *BVP) (constraint *matrix.SparseMatrix, err error) {

	constraint = matrix.ZerosSparse(bvp.ODE.P*bvp.N, bvp.ODE.P*bvp.N)

	var constraint_i1 *matrix.DenseMatrix
	var constraint_i2 *matrix.DenseMatrix

	dfdxNow, err := bvp.ODE.Dfdx(bvp.X[0], bvp.T[0], bvp.Beta)

	if err != nil {
		return
	}

	for i := 1; i < bvp.N; i++ {

		dfdxBefore := matrix.MakeDenseCopy(dfdxNow)

		dfdxNow, err = bvp.ODE.Dfdx(bvp.X[i], bvp.T[i], bvp.Beta)

		if err != nil {
			return
		}

		constraint_i1 = matrix.Difference(
			matrix.Scaled(matrix.Eye(bvp.ODE.P), -1),
			matrix.Scaled(dfdxBefore, (bvp.T[i]-bvp.T[i-1])/2),
		)
		constraint_i2 = matrix.Difference(
			matrix.Eye(bvp.ODE.P),
			matrix.Scaled(dfdxBefore, (bvp.T[i]-bvp.T[i-1])/2),
		)

		for r := 0; r < bvp.ODE.P; r++ {
			for c := 0; c < bvp.ODE.P; c++ {
				constraint.Set(r+(i-1)*bvp.ODE.P, c+(i-1)*bvp.ODE.P, constraint_i1.Get(r, c))
				constraint.Set(r+(i-1)*bvp.ODE.P, c+i*bvp.ODE.P, constraint_i2.Get(r, c))
			}
		}
	}

	for r := 0; r < bvp.ODE.P; r++ {
		for c := 0; c < bvp.ODE.P; c++ {
			constraint.Set(r+(bvp.N-1)*bvp.ODE.P, c, bvp.B0.Get(r, c))
			constraint.Set(r+(bvp.N-1)*bvp.ODE.P, c+(bvp.N-1)*bvp.ODE.P, bvp.B1.Get(r, c))
		}
	}

	return
}

func ConstraintMatrixPermuted(bvp *BVP) (constraint *matrix.SparseMatrix, err error) {

	constraint = matrix.ZerosSparse(bvp.ODE.P*bvp.N, bvp.ODE.P*bvp.N)

	var constraint_i1 *matrix.DenseMatrix
	var constraint_i2 *matrix.DenseMatrix

	dfdxNow, err := bvp.ODE.Dfdx(bvp.X[0], bvp.T[0], bvp.Beta)

	if err != nil {
		return
	}

	for i := 1; i < bvp.N; i++ {

		dfdxBefore := matrix.MakeDenseCopy(dfdxNow)

		dfdxNow, err = bvp.ODE.Dfdx(bvp.X[i], bvp.T[i], bvp.Beta)

		if err != nil {
			return
		}

		constraint_i1 = matrix.Difference(
			matrix.Scaled(matrix.Eye(bvp.ODE.P), -1),
			matrix.Scaled(dfdxBefore, (bvp.T[i]-bvp.T[i-1])/2),
		)
		constraint_i2 = matrix.Difference(
			matrix.Eye(bvp.ODE.P),
			matrix.Scaled(dfdxBefore, (bvp.T[i]-bvp.T[i-1])/2),
		)

		for r := 0; r < bvp.ODE.P; r++ {
			for c := 0; c < bvp.ODE.P; c++ {
				if i == 1 {
					constraint.Set(r+(i-1)*bvp.ODE.P, c+(bvp.N-1)*bvp.ODE.P, constraint_i1.Get(r, c))
				} else {
					constraint.Set(r+(i-1)*bvp.ODE.P, c+(i-2)*bvp.ODE.P, constraint_i1.Get(r, c))
				}
				constraint.Set(r+(i-1)*bvp.ODE.P, c+(i-1)*bvp.ODE.P, constraint_i2.Get(r, c))
			}
		}
	}

	for r := 0; r < bvp.ODE.P; r++ {
		for c := 0; c < bvp.ODE.P; c++ {
			constraint.Set(r+(bvp.N-1)*bvp.ODE.P, c+(bvp.N-1)*bvp.ODE.P, bvp.B0.Get(r, c))
			constraint.Set(r+(bvp.N-1)*bvp.ODE.P, c+(bvp.N-2)*bvp.ODE.P, bvp.B1.Get(r, c))
		}
	}

	return
}