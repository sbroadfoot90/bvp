package bvp

import (
 	"github.com/skelterjohn/go.matrix"
)

type BVP struct {
	ODE ODE
	X matrix.Matrix
	T []float64
	B0, B1 matrix.Matrix
}

func (bvp *BVP) solve (bta, b matrix.Matrix) {
	
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