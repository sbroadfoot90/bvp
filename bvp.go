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

func (bvp *BVP) Solve() error {
	tolerance := 1e-6
	maxiter := 500

	constraintBlocks, err := ConstraintVectorBlocks(bvp)
	if err != nil {
		return err
	}

	cost := sumOfSquares(constraintBlocks) // compute cost
	for i := 0; i < maxiter; i++ {
		delta, err := getDelta(bvp)

		if !exceedsTolerance(delta, tolerance) {
			return nil
		}

		// save old values
		costold := cost

		xold := make([]matrix.Matrix, bvp.N, bvp.N)
		for i := 0; i < bvp.N; i++ {
			xold[i] = matrix.MakeDenseCopy(bvp.X[i])
		}
		var alpha float64 = 1

		for i := 0; i < bvp.N; i++ {
			bvp.X[i].Subtract(delta[i]) // update x
		}

		constraintBlocks, err = ConstraintVectorBlocks(bvp)
		if err != nil {
			return err
		}
		cost := sumOfSquares(constraintBlocks) // compute cost

		if cost < costold {
			continue
		}

		// revert to old values
		for i := 0; i < bvp.N; i++ {
			bvp.X[i] = matrix.MakeDenseCopy(xold[i])
		}

		var dcost float64 = 0
		constraintBlocks, err = ConstraintVectorBlocks(bvp)
		A, B, err := ConstraintMatrixBlocks(bvp)

		for i := 0; i < bvp.N; i++ {
			var tempMatrix *matrix.DenseMatrix
			if i < bvp.N-1 {
				tempMatrix = matrix.Sum(matrix.Product(A[i], delta[i]), matrix.Product(B[i], delta[i+1]))
			} else {
				tempMatrix = matrix.Sum(matrix.Product(bvp.B0, delta[0]), matrix.Product(bvp.B1, delta[bvp.N-1]))
			}
			for j := 0; j < bvp.ODE.P; j++ {
				dcost -= constraintBlocks[i].Get(j, 0) * tempMatrix.Get(j, 0)
			}
		}

		for {
			if dcost > 0 {
				// display('warning, dcost positive')
				// display(dcost)
				alpha = -alpha
				dcost = -dcost
			} else {
				// scale delta if cost increased
				alpha = alpha / 2
				dcost = dcost / 2
			}

			for i := 0; i < bvp.N; i++ {
				delta[i].Scale(alpha)
			}

			if !exceedsTolerance(delta, tolerance) {
				return nil // converged
			}

			for i := 0; i < bvp.N; i++ {
				bvp.X[i].Subtract(delta[i]) // update x
			}

			constraintBlocks, err = ConstraintVectorBlocks(bvp)
			if err != nil {
				return err
			}
			cost := sumOfSquares(constraintBlocks) // compute cost

			if cost < costold {
				break
			}

			// revert to old values
			for i := 0; i < bvp.N; i++ {
				bvp.X[i] = matrix.MakeDenseCopy(xold[i])
			}
		}
	}

	return ConvergeError("Warning: BVP solver did not terminate")
}

func (bvp *BVP) SolveIVP(initialGuess matrix.Matrix) error {
	bvp.X[0] = initialGuess
	niter := 10

	for i := 1; i < len(bvp.X); i++ {
		fBefore, err := bvp.ODE.F(bvp.X[i-1], bvp.T[i-1], bvp.Beta)
		if err != nil {
			return err
		}
		dt := bvp.T[i] - bvp.T[i-1]
		bvp.X[i] = matrix.Sum(bvp.X[i-1], matrix.Scaled(fBefore, dt))

		for j := 0; j < niter; j++ {
			fNow, err := bvp.ODE.F(bvp.X[i], bvp.T[i], bvp.Beta)
			if err != nil {
				return err
			}
			dfdxNow, err := bvp.ODE.Dfdx(bvp.X[i], bvp.T[i], bvp.Beta)
			if err != nil {
				return err
			}
			C := matrix.Difference(matrix.Eye(bvp.ODE.P), matrix.Scaled(dfdxNow, dt/2))
			c := matrix.Difference(matrix.Difference(bvp.X[i], bvp.X[i-1]), matrix.Scaled(matrix.Sum(fNow, fBefore), dt/2))
			delta, err := C.Solve(c)
			if err != nil {
				return err
			}
			bvp.X[i].Subtract(delta)
		}
	}

	return nil
}

func ConstraintMatrixBlocks(bvp *BVP) (A, B []*matrix.DenseMatrix, err error) {
	A = make([]*matrix.DenseMatrix, bvp.N-1, bvp.N-1)
	B = make([]*matrix.DenseMatrix, bvp.N-1, bvp.N-1)

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
	}
	return
}

func ConstraintVectorBlocks(bvp *BVP) (constraint []*matrix.DenseMatrix, err error) {

	constraint = make([]*matrix.DenseMatrix, bvp.N, bvp.N)

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

		constraint[i-1] = matrix.Difference(
			matrix.Difference(bvp.X[i], bvp.X[i-1]),
			matrix.Scaled(matrix.Sum(fBefore, fNow), (bvp.T[i]-bvp.T[i-1])/2),
		)
	}

	constraint[bvp.N-1] = matrix.Difference(
		matrix.Sum(matrix.Product(bvp.B0, bvp.X[0]), matrix.Product(bvp.B1, bvp.X[bvp.N-1])),
		bvp.B,
	)

	return
}

func getDelta(bvp *BVP) (delta []*matrix.DenseMatrix, err error) {
	// delta = solve(C(x), c(x))
	// C(x) %*% delta = c(x)

	A, B, err := ConstraintMatrixBlocks(bvp)
	if err != nil {
		return
	}

	c, err := ConstraintVectorBlocks(bvp)
	if err != nil {
		return
	}

	C, D, U := rightOrthogonalFactorisation(A, B, bvp.N, bvp.ODE.P)
	rqTransformation(A, B, U, c, bvp.N, bvp.ODE.P)

	hgb1b0 := matrix.Zeros(bvp.ODE.P*2, bvp.ODE.P*2)
	hgb1b0.SetMatrix(0, 0, B[bvp.N-2])
	hgb1b0.SetMatrix(0, bvp.ODE.P, D[bvp.N-2])
	hgb1b0.SetMatrix(bvp.ODE.P, 0, matrix.MakeDenseCopy(bvp.B1))
	hgb1b0.SetMatrix(bvp.ODE.P, bvp.ODE.P, matrix.MakeDenseCopy(bvp.B0))
	smallc := matrix.Zeros(bvp.ODE.P*2, 1)
	smallc.SetMatrix(0, 0, c[bvp.N-2])
	smallc.SetMatrix(bvp.ODE.P, 0, c[bvp.N-1])

	deltaends, err := hgb1b0.SolveDense(smallc)

	delta = make([]*matrix.DenseMatrix, bvp.N, bvp.N)

	for i := 0; i < bvp.N; i++ {
		delta[i] = matrix.Zeros(bvp.ODE.P, 1)
	}

	delta[0].SetMatrix(0, 0, deltaends.GetMatrix(bvp.ODE.P, 0, bvp.ODE.P, 1))
	delta[bvp.N-1].SetMatrix(0, 0, deltaends.GetMatrix(0, 0, bvp.ODE.P, 1))

	rightBackSubstitute(B, C, D, U, c, delta, bvp.N, bvp.ODE.P)

	return
}

func SetOptimalBoundaryMatrices(bvp *BVP) {

	A, B, err := ConstraintMatrixBlocks(bvp)
	if err != nil {
		return
	}

	_, D, _ := rightOrthogonalFactorisation(A, B, bvp.N, bvp.ODE.P)

	bvp.B0, bvp.B1 = calculateBoundaryMatrices(B, D, bvp.N, bvp.ODE.P)

	bvp.B = matrix.Sum(matrix.Product(bvp.B0, bvp.X[0]), matrix.Product(bvp.B1, bvp.X[bvp.N-1]))
}
