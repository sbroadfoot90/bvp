package bvp

import (
	"encoding/csv"
	"github.com/sbroadfoot90/go.matrix"
	"math"
	"os"
	"strconv"
)

func WriteMatrix(mat matrix.MatrixRO, filename string) (err error) {
	file, err := os.Create(filename)
	defer file.Close()

	if err != nil {
		return
	}

	csvWriter := csv.NewWriter(file)

	strMatrix := make([][]string, mat.Rows())

	for rowIndex := range strMatrix {
		strMatrix[rowIndex] = make([]string, mat.Cols())
		for colIndex := range strMatrix[rowIndex] {
			strMatrix[rowIndex][colIndex] = strconv.FormatFloat(mat.Get(rowIndex, colIndex), 'f', -1, 64)
		}
	}
	err = csvWriter.WriteAll(strMatrix)
	csvWriter.Flush()
	return
}

func WriteMatrices(mats []matrix.Matrix, filename string) (err error) {
	file, err := os.Create(filename)
	defer file.Close()

	if err != nil {
		return
	}

	csvWriter := csv.NewWriter(file)

	for i := range mats {
		strMatrix := make([][]string, mats[i].Rows())
		for rowIndex := range strMatrix {
			strMatrix[rowIndex] = make([]string, mats[i].Cols())
			for colIndex := range strMatrix[rowIndex] {
				strMatrix[rowIndex][colIndex] = strconv.FormatFloat(mats[i].Get(rowIndex, colIndex), 'f', -1, 64)
			}
		}
		err = csvWriter.WriteAll(strMatrix)
		csvWriter.Flush()
	}
	return
}

func LinearFDfdx(
	a, q func(float64, matrix.MatrixRO) matrix.Matrix,
) (
	f, dfdx func(matrix.MatrixRO, float64, matrix.MatrixRO) matrix.Matrix,
) {

	f = func(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
		return matrix.Sum(matrix.Product(a(t, beta), x), q(t, beta))
	}

	dfdx = func(x matrix.MatrixRO, t float64, beta matrix.MatrixRO) matrix.Matrix {
		return a(t, beta)
	}

	return
}

func RORFAC(A, B []*matrix.DenseMatrix, n, m int) (C, D, U []*matrix.DenseMatrix) {
	//............................................................................
	//Orthogonal factorization of block bidiagonal matrix to upper triangular form
	//  B(1)    ........... A(1)        B(1) C(1) ............D(1)
	//  A(2) B(2)  ........                  B(2) C(2) .......D(2)
	//       A(3) B(3)  ...        ->             B(3) .......D(3)
	//  ........................        ..........................
	//           A(n-1) B(n-1)                        B(n-1) D(n-1)
	//Here the C(i) allows for fill as do entries D(2),... , D(n-1), while
	//U(i) contains diagonal of upper triangular matrix. Transformation is held
	//in lower triangle of B(i) and A(i+1), i=1,2,...n-2. B(n-1) is not
	//swept out at this stage.
	//application is to difference relation A(i)x(i)+B(i)x(i+1)=q(i), D(1)=A(1).
	//............................................................................
	// #real(kind=8),dimension(:,:,:),intent(inout)::A,B,C,D
	// #real(kind=8),dimension(:,:),intent(inout)::U
	// #integer,intent(in)::n,m
	// #!allocate(C(n - 1, m, m))

	C = make([]*matrix.DenseMatrix, n-1, n-1)
	D = make([]*matrix.DenseMatrix, n-1, n-1)
	U = make([]*matrix.DenseMatrix, n-2, n-2)

	for i := 0; i < n-1; i++ {
		C[i] = matrix.Zeros(m, m)
		D[i] = matrix.Zeros(m, m)
	}

	for i := 0; i < n-2; i++ {
		U[i] = matrix.Zeros(m, 1)
	}

	D[0] = A[0].Copy()

	for i := 0; i < n-2; i++ {
		for j := 0; j < m; j++ {
			var ss1 float64 = 0
			for k := j; k < m; k++ {
				ss1 = ss1 + math.Pow(B[i].Get(k, j), 2)
			}

			for k := 0; k < m; k++ {
				ss1 = ss1 + math.Pow(A[i+1].Get(k, j), 2)
			}

			U[i].Set(j, 0, -math.Sqrt(ss1)) //#!assumes B(i,j,j)>0
			if B[i].Get(j, j) < 0 {
				U[i].Set(j, 0, -U[i].Get(j, 0))
			}

			B[i].Set(j, j, B[i].Get(j, j)-U[i].Get(j, 0))
			var s float64 = -U[i].Get(j, 0) * B[i].Get(j, j)

			if j < m-1 {
				for k := j + 1; k < m; k++ {
					ss1 = 0.
					for l := j; l < m; l++ {
						ss1 = ss1 + B[i].Get(l, j)*B[i].Get(l, k)
					}

					for l := 0; l < m; l++ {
						ss1 = ss1 + A[i+1].Get(l, j)*A[i+1].Get(l, k)
					}

					ss1 = ss1 / s

					for l := j; l < m; l++ {
						B[i].Set(l, k, B[i].Get(l, k)-ss1*B[i].Get(l, j))
					}
					for l := 0; l < m; l++ {
						A[i+1].Set(l, k, A[i+1].Get(l, k)-ss1*A[i+1].Get(l, j))
					}
				}
			}

			for k := 0; k < m; k++ {
				var ss1, ss2 float64 = 0, 0
				for l := j; l < m; l++ {
					ss1 = ss1 + B[i].Get(l, j)*C[i].Get(l, k)
					ss2 = ss2 + B[i].Get(l, j)*D[i].Get(l, k)
				}

				for l := 0; l < m; l++ {
					ss1 = ss1 + A[i+1].Get(l, j)*B[i+1].Get(l, k)
					ss2 = ss2 + A[i+1].Get(l, j)*D[i+1].Get(l, k)
				}
				ss1 = ss1 / s
				ss2 = ss2 / s
				for l := j; l < m; l++ {
					C[i].Set(l, k, C[i].Get(l, k)-ss1*B[i].Get(l, j))
					D[i].Set(l, k, D[i].Get(l, k)-ss2*B[i].Get(l, j))
				}
				for l := 0; l < m; l++ {
					B[i+1].Set(l, k, B[i+1].Get(l, k)-ss1*A[i+1].Get(l, j))
					D[i+1].Set(l, k, D[i+1].Get(l, k)-ss2*A[i+1].Get(l, j))
				}
			}
		}
	}

	return
}

func RQTRANS(A, B, U, q []*matrix.DenseMatrix, n, m int) {
	// #!...............................................................................
	// #!application of orthogonal transformation RORFAC to rhs
	// #!...............................................................................
	// #real(kind=8),dimension(:,:,:),intent(inout)::A,B
	// #real(kind=8),dimension(:,:),intent(inout)::U,q
	// #integer,intent(in)::n,m
	for i := 0; i < n-2; i++ {
		for j := 0; j < m; j++ {
			var ss1 float64 = 0.
			for k := j; k < m; k++ {
				ss1 = ss1 + B[i].Get(k, j)*q[i].Get(k, 0)
			}
			for k := 0; k < m; k++ {
				ss1 = ss1 + A[i+1].Get(k, j)*q[i+1].Get(k, 0)
			}
			ss1 = ss1 / (-U[i].Get(j, 0) * B[i].Get(j, j))

			for k := j; k < m; k++ {
				q[i].Set(k, 0, q[i].Get(k, 0)-ss1*B[i].Get(k, j))
			}
			for k := 0; k < m; k++ {
				q[i+1].Set(k, 0, q[i+1].Get(k, 0)-ss1*A[i+1].Get(k, j))
			}
		}
	}
}

func RBKSUB(B, C, D, U, q, xc []*matrix.DenseMatrix, n, m int) {
	// #!..............................................................
	// #!given starting values xc(1),xc(n) perform back substitution
	// #!..............................................................
	// #real(kind=8),dimension(:,:,:),intent(inout)::B,C,D
	// #real(kind=8),dimension(:,:),intent(inout)::U,q,xc
	// #integer,intent(in)::n,m
	for i := n - 3; i >= 0; i-- {
		for j := m - 1; j >= 0; j-- {
			var ss1 float64 = 0.
			for k := 0; k < m; k++ {
				ss1 = ss1 + C[i].Get(j, k)*xc[i+2].Get(k, 0) + D[i].Get(j, k)*xc[0].Get(k, 0)
			}
			ss1 = q[i].Get(j, 0) - ss1
			if j < m-1 {
				for k := j + 1; k < m; k++ {
					ss1 = ss1 - B[i].Get(j, k)*xc[i+1].Get(k, 0)
				}
			}
			xc[i+1].Set(j, 0, ss1/U[i].Get(j, 0))
		}
	}
}

func BMCALC(B, D []*matrix.DenseMatrix, n, m int) (B1, Bn *matrix.DenseMatrix) {
	// !....................................................................
	// !use orthogonal factorization of constraint matrix [B,D] to set up boundary
	// !conditions
	// !....................................................................
	// real(kind=8),dimension(:,:,:),intent(inout)::B,D
	// real(kind=8),dimension(:,:),intent(inout)::UL,B1,Bn
	// integer,intent(in)::n,m
	// 
	// real(kind=8),dimension(:,:),allocatable::ABC
	// real(kind=8),dimension(:),allocatable::UD

	// allocate(ABC(2 * m, 3 * m), UD(m))

	ABC := matrix.Zeros(2*m, 3*m)
	UD := make([]float64, m)

	// !cast input data
	// ![B^T I ]=[A_{1,1} A_{1,2} A_{1,3}]
	// ![D^T  I] [A_{2,1} A_{2,2} A_{2,3}]

	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			ABC.Set(i, j, B[n-2].Get(j, i))   //!A_{1,1}=B^T
			ABC.Set(i+m, j, D[n-2].Get(j, i)) //!A_{2,1}=D^T
		}
		ABC.Set(i, m+i, 1.)     //!A_{1,2}=I
		ABC.Set(m+i, 2*m+i, 1.) //!A_{2,3}=I
	}

	//!make orthogonal factorization
	//![B^T I ]->[U Q_1^T]
	//![D^T  I]  [0 Q_2^T]

	for i := 0; i < m; i++ {
		var ss1 float64 = 0.
		for j := i; j < 2*m; j++ {
			ss1 = ss1 + math.Pow(ABC.Get(j, i), 2)
		}
		UD[i] = -math.Sqrt(ss1) //!assumes ABC(i,i)>0
		if ABC.Get(i, i) < 0. {
			UD[i] = -UD[i]
		}
		ABC.Set(i, i, ABC.Get(i, i)-UD[i])
		for j := i + 1; j < 3*m; j++ {
			ss1 = 0.
			for k := i; k < 2*m; k++ {
				ss1 = ss1 + ABC.Get(k, i)*ABC.Get(k, j)
			}
			ss1 = ss1 / (-UD[i] * ABC.Get(i, i))
			for k := i; k < 2*m; k++ {
				ABC.Set(k, j, ABC.Get(k, j)-ss1*ABC.Get(k, i))
			}
		}
	}
	//!save U^T to apply to rhs q(n-1,*)
	//![ B  D][xc(n)]=[q]-> [Q_1^T][xc(n)]=[U^-Tq]
	//![Q_2^T][xc(1)] [b]   [Q_2^T][xc(1)] [  b  ]

	// UL.Set(0, 0, UD[0])
	// for i := 1; i < m; i++ {
	// 	UL.Set(i, i, UD[i])
	// 	for j := 0; j < i-1; j++ {
	// 		UL.Set(i, j, ABC.Get(j, i))
	// 	}
	// }

	//!copy orthogonal transformation to boundary matrices
	//![Q_1^T]->[ B D ]
	//![Q_2^T]  [Bn B1]
	B1 = matrix.Zeros(3, 3)
	Bn = matrix.Zeros(3, 3)
	for i := 0; i < m; i++ {

		for j := 0; j < m; j++ {
			B[n-2].Set(i, j, ABC.Get(i, m+j))
			D[n-2].Set(i, j, ABC.Get(i, 2*m+j))
			Bn.Set(i, j, ABC.Get(m+i, m+j))
			B1.Set(i, j, ABC.Get(m+i, 2*m+j))
		}
	}

	return B1, Bn
}

func exceedsTolerance(matrices []*matrix.DenseMatrix, tolerance float64) bool {
	for t := range matrices {
		for i := 0; i < matrices[t].Rows(); i++ {
			for j := 0; j < matrices[t].Cols(); j++ {
				if math.Abs(matrices[t].Get(i, j)) > tolerance {
					return true
				}
			}
		}
	}
	return false
}

func sumOfSquares(matrices []*matrix.DenseMatrix) (ss float64) {
	for t := range matrices {
		for i := 0; i < matrices[t].Rows(); i++ {
			for j := 0; j < matrices[t].Cols(); j++ {
				ss += math.Pow(matrices[t].Get(i, j), 2)
			}
		}
	}
	return
}

func sign(x float64) float64 {
	if x < 0 {
		return -1
	}
	return 1
}

func spiral(x, y, shell, side int) (int, int, int, int) {

	if side == 1 {
		x++
		if x == shell {
			side = 2
		}
	} else if side == 2 {
		y--
		if y == -shell {
			side = 3
		}
	} else if side == 3 {
		x--
		if x == -shell {
			side = 4
		}
	} else {
		y++
		if y == shell+1 {
			shell++
			side = 1
		}
	}

	return x, y, shell, side
}