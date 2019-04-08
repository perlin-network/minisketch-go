package minisketch

import (
	"github.com/stretchr/testify/assert"
	"testing"
	"testing/quick"
)

func TestBerlekampMassey(t *testing.T) {
	sketch := Poly{8160, 22369280, 75107501056, 300239975088128, 1384206083625254912, 1535815439233325851}
	coeffs := BerlekampMassey(sketch)

	for left, right := 0, len(coeffs)-1; left < right; left, right = left+1, right-1 {
		coeffs[left], coeffs[right] = coeffs[right], coeffs[left]
	}

	assert.ElementsMatch(t, Poly{5000, 3000, 2000}, FindRoots(coeffs))
}

func TestSketch(t *testing.T) {
	f := func(a []Element, b []Element) bool {
		alice := New(len(a) + len(b)).Add(a...)
		bob := New(len(a) + len(b)).Add(b...)

		if !assert.ElementsMatch(t, a, alice.Decode()) {
			return false
		}

		if !assert.ElementsMatch(t, b, bob.Decode()) {
			return false
		}

		return assert.ElementsMatch(t, setDifference(a, b), bob.Merge(alice).Decode())
	}

	assert.NoError(t, quick.Check(f, nil))
}

func TestDecodeSketchOverCapacity(t *testing.T) {
	f := func(a []Element) bool {
		return assert.Empty(t, New(0).Add(a...).Decode())
	}

	assert.NoError(t, quick.Check(f, nil))
}

func TestDecodeSketch(t *testing.T) {
	a := []Element{2000, 4000, 5000}
	b := []Element{4000, 5000, 1000}

	alice := New(4).Add(a...)
	bob := New(12).Add(b...)

	assert.ElementsMatch(t, a, alice.Decode())
	assert.ElementsMatch(t, b, bob.Decode())

	assert.ElementsMatch(t, setDifference(a, b), bob.Merge(alice).Decode())
}

func setDifference(a, b []Element) []Element {
	visited := make(map[Element]struct{})

	for _, elem := range a {
		_, seen := visited[elem]

		if seen {
			delete(visited, elem)
			continue
		}

		visited[elem] = struct{}{}
	}

	var expected []Element

	for _, elem := range b {
		_, seen := visited[elem]

		if seen {
			delete(visited, elem)
			continue
		}

		visited[elem] = struct{}{}
	}

	for elem := range visited {
		expected = append(expected, elem)
	}

	return expected
}
