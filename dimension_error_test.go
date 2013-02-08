package bvp

import (
	"testing"
)

func TestDimensionError(t *testing.T) {

	testError := NewDimensionError("x", 4, 6, 5, 7)

	if testError.Error() != "Variable x, received dimensions (5, 7), expected dimensions (4, 6)" {
		t.Errorf("DimensionError not returning expected string representation")
	}
}
