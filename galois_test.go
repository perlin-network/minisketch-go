package minisketch

import (
	"github.com/stretchr/testify/assert"
	"testing"
	"testing/quick"
)

func TestGaloisAxioms(t *testing.T) {
	f := func(a Element) bool {
		mul := a.Mul(a)
		exp := a.Exp(2)

		if !assert.Equal(t, mul, exp) {
			return false
		}

		div := exp.Div(a)

		if !assert.Equal(t, a, div) {
			return false
		}

		sqrt := exp.Sqrt()
		exp2 := sqrt.Exp(2)

		if !assert.Equal(t, exp, exp2) {
			return false
		}

		return true
	}

	assert.NoError(t, quick.Check(f, &quick.Config{MaxCount: 1000}))
}
