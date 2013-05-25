package bvp

import (
	"fmt"
	"github.com/sbroadfoot90/go.matrix"
	"math"
)

func Surfplotb(bvp BVP, O matrix.Matrix, resolution int, b1index int, b1i, b1f float64, b2index int, b2i, b2f float64) (b1, b2, costMatrix *matrix.DenseMatrix) {
	bvp.Solve()

	y := make([]*matrix.DenseMatrix, bvp.N)
	xtrue := make([]*matrix.DenseMatrix, bvp.N)
	for i := 0; i < bvp.N; i++ {
		xtrue[i] = matrix.MakeDenseCopy(bvp.X[i])
		y[i] = matrix.Product(O, bvp.X[i])
	}

	costMatrix = matrix.Zeros(2*resolution+1, 2*resolution+1)

	b1 = matrix.Zeros(2*resolution+1, 1)
	b2 = matrix.Zeros(2*resolution+1, 1)
	for i := 0; i < 2*resolution+1; i++ {
		b1.Set(i, 0, (b1f-b1i)*float64(i)/float64(2*resolution)+b1i)
		b2.Set(i, 0, (b2f-b2i)*float64(i)/float64(2*resolution)+b2i)
	}

	d1 := 0
	d2 := 0
	shell := 0
	side := 4

	for i := 0; i < (2*resolution+1)*(2*resolution+1); i++ {
		bvp.B.Set(b1index, 0, b1.Get(d1+resolution, 0))
		bvp.B.Set(b2index, 0, b2.Get(d2+resolution, 0))
		
		for j := 0; j < bvp.N; j++ {
			bvp.X[j] = matrix.MakeDenseCopy(xtrue[j])
		}
		
		bvp.Solve()

		var cost float64 = 0

		for t := 0; t < bvp.N; t++ {
			residual := matrix.Difference(y[t], matrix.Product(O, bvp.X[t]))
			for j := 0; j < residual.Rows(); j++ {
				for k := 0; k < residual.Cols(); k++ {
					cost += math.Pow(residual.Get(j, k), 2)
				}
			}
		}

		costMatrix.Set(d2+resolution, d1+resolution, cost/float64(2)/float64(bvp.N))
		d1, d2, shell, side = spiral(d1, d2, shell, side)
		fmt.Printf("%d, %d, %d, \n", i, d1, d2)
	}

	return
}
