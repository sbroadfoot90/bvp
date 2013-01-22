package bvp

import (
	"testing"
	"github.com/skelterjohn/go.matrix"
)


func TestWriteMatrix(t *testing.T) {
	
	a, err := matrix.ParseMatlab("[1 2 3; 4 2 1]")
	
	if err != nil {
		t.Errorf(err.Error());
		return
	}
	err = writeMatrix(a, "asdf.csv")
	
	if err != nil {
		t.Errorf(err.Error());
	}
	
	return
}