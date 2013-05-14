package bvp

import (
	"fmt"
)

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

type MatrixError string

func (me MatrixError) Error() string {
	return string(me)
}
