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

    C = make([]*matrix.DenseMatrix, n - 1, n - 1)
	D = make([]*matrix.DenseMatrix, n - 1, n - 1)
    U = make([]*matrix.DenseMatrix, n - 2, n - 2)

    for i := 0; i < n - 1; i++ {
        C[i] = matrix.Zeros(m, m)
        D[i] = matrix.Zeros(m, m)
    }

	for i := 0; i < n - 2; i++ {
        U[i] = matrix.Zeros(m, 1)
    }
    
   D[0] = A[0].Copy()

    for i := 0; i < n - 2; i++ {
        for j := 0; j < m; j++ {
            var ss1 float64 = 0
            for k := j; k < m; k++ {
                ss1 = ss1 + math.Pow(B[i].Get(k, j), 2)
            }
            
            for k := 0; k < m; k++ {
                ss1 = ss1 + math.Pow(A[i + 1].Get(k, j), 2)
            }
			
            U[i].Set(j, 0, -math.Sqrt(ss1))  //#!assumes B(i,j,j)>0
            if ( B[i].Get(j, j) < 0 ) {
                U[i].Set(j, 0, -U[i].Get(j, 0))
            }

            B[i].Set(j, j, B[i].Get(j, j) - U[i].Get(j,0))
            var s float64 = -U[i].Get(j, 0) * B[i].Get(j, j)
            
            if ( j < m - 1 ) {
                for k := j + 1; k < m; k++ {
                    ss1 = 0.
                    for l := j; l < m; l++ {
                        ss1 = ss1 + B[i].Get(l, j) * B[i].Get(l, k)
                    }
                    
                    for l := 0; l < m; l++ {
                        ss1 = ss1 + A[i + 1].Get(l, j) * A[i + 1].Get(l, k)
                    }
                    
                    ss1 = ss1 / s
                
                    for l := j; l < m; l++ {
                        B[i].Set(l, k, B[i].Get(l, k) - ss1 * B[i].Get(l, j))
                    }
                    for l := 0; l < m; l++ {
                        A[i + 1].Set(l, k, A[i + 1].Get(l, k) - ss1 * A[i + 1].Get(l, j))
                    }
                }
            }
            
            for k := 0; k < m; k++ {
                var ss1, ss2 float64 = 0, 0
                for l := j; l < m; l ++ {
                    ss1 = ss1 + B[i].Get(l, j) * C[i].Get(l, k)
                    ss2 = ss2 + B[i].Get(l, j) * D[i].Get(l, k)
                }
                
                for l := 0; l < m; l++ {
                    ss1 = ss1 + A[i + 1].Get(l, j) * B[i + 1].Get(l, k)
                    ss2 = ss2 + A[i + 1].Get(l, j) * D[i + 1].Get(l, k)
                }
                ss1 = ss1 / s
                ss2 = ss2 / s
                for l := j; l < m; l++ {
                    C[i].Set(l, k, C[i].Get(l, k) - ss1 * B[i].Get(l, j))
                    D[i].Set(l, k, D[i].Get(l, k) - ss2 * B[i].Get(l, j))
                }
                for l := 0; l < m; l++ {
                    B[i + 1].Set(l, k, B[i + 1].Get(l, k) - ss1 * A[i + 1].Get(l, j))
                    D[i + 1].Set(l, k, D[i + 1].Get(l, k) - ss2 * A[i + 1].Get(l, j))
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
    for i := 0; i < n - 2; i++ {
        for j := 0; j < m; j++ {
            var ss1 float64 = 0.
            for k := j; k < m; k++ {
                ss1 = ss1 + B[i].Get(k, j) * q[i].Get(k, 0)
            }
            for k := 0; k < m; k++ {
                ss1 = ss1 + A[i + 1].Get(k, j) * q[i + 1].Get(k, 0)
            }
            ss1 = ss1 / (-U[i].Get(j, 0) * B[i].Get(j, j))
            
            for k := j; k < m; k++ {
                q[i].Set(k, 0, q[i].Get(k, 0) - ss1 * B[i].Get(k, j))
            }
            for k := 0; k < m; k++ {
                q[i + 1].Set(k, 0, q[i + 1].Get(k, 0) - ss1 * A[i + 1].Get(k, j))
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
                ss1 = ss1 + C[i].Get(j, k) * xc[i + 2].Get(k, 0) + D[i].Get(j, k) * xc[0].Get(k, 0)
            }
            ss1 = q[i].Get(j, 0) - ss1
            if j < m - 1 {
                for k := j + 1; k < m; k++ {
                    ss1 = ss1 - B[i].Get(j, k) * xc[i + 1].Get(k, 0)
                }
            }
            xc[i + 1].Set(j, 0, ss1 / U[i].Get(j, 0))
         }
    }
}



func sparseEye(span int) (I *matrix.SparseMatrix) {
	I = matrix.ZerosSparse(span, span)

	for i := 0; i < span; i++ {
		I.Set(i, i, 1)
	}

	return
}

func sign(x float64) float64 {
	if x < 0 {
		return -1
	}
	return 1
}

func SparseQR(A *matrix.SparseMatrix) (Q, R *matrix.SparseMatrix) {
	m := A.Rows()
	n := A.Cols()
	QR := A.Copy()
	Q = matrix.ZerosSparse(m, n)
	R = matrix.ZerosSparse(m, n)
	i, j, k := 0, 0, 0
	norm := float64(0.0)
	s := float64(0.0)

	for k = 0; k < n; k++ {
		norm = 0
		for i = k; i < m; i++ {
			norm = math.Hypot(norm, QR.Get(i, k))
		}

		if norm != 0.0 {
			if QR.Get(k, k) < 0 {
				norm = -norm
			}

			for i = k; i < m; i++ {
				QR.Set(i, k, QR.Get(i, k)/norm)
			}
			QR.Set(k, k, QR.Get(k, k)+1.0)

			for j = k + 1; j < n; j++ {
				s = 0.0
				for i = k; i < m; i++ {
					s += QR.Get(i, k) * QR.Get(i, j)
				}
				s = -s / QR.Get(k, k)
				for i = k; i < m; i++ {
					QR.Set(i, j, QR.Get(i, j)+s*QR.Get(i, k))

					if i < j {
						R.Set(i, j, QR.Get(i, j))
					}
				}

			}
		}

		R.Set(k, k, -norm)

	}

	//Q Matrix:
	i, j, k = 0, 0, 0

	for k = n - 1; k >= 0; k-- {
		Q.Set(k, k, 1.0)
		for j = k; j < n; j++ {
			if QR.Get(k, k) != 0 {
				s = 0.0
				for i = k; i < m; i++ {
					s += QR.Get(i, k) * Q.Get(i, j)
				}
				s = -s / QR.Get(k, k)
				for i = k; i < m; i++ {
					Q.Set(i, j, Q.Get(i, j)+s*QR.Get(i, k))
				}
			}
		}
	}

	return
}

func solveUpper(A matrix.MatrixRO, b matrix.MatrixRO) *matrix.DenseMatrix {
	x := make([]float64, A.Cols())
	for i := A.Rows() - 1; i >= 0; i-- {
		x[i] = b.Get(i, 0)
		for j := i + 1; j < A.Cols(); j++ {
			x[i] -= x[j] * A.Get(i, j)
		}
		x[i] /= A.Get(i, i)
	}
	return matrix.MakeDenseMatrix(x, A.Cols(), 1)
}

func solve(A *matrix.SparseMatrix, b matrix.MatrixRO) matrix.Matrix {
	Q, R := SparseQR(A)
	Qtb, _ := Q.Transpose().Times(b)
	
	return solveUpper(R, Qtb)
}