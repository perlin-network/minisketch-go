package minisketch

// FindRoots recursively factors a polynomial whose coefficients are over
// GF(2^64), and returns its roots.
//
// It implements Berlekamp's Trace algorithm, and returns nil should no roots
// be found.
func FindRoots(poly Poly) []Element {
	return findRoots([]Poly{poly}, 0, false, 0, 1)
}

func findRoots(stack []Poly, pos int, factorizable bool, depth int, a Element) []Element {
	var roots []Element

	if len(stack[pos]) == 2 {
		return append(roots, stack[pos][0])
	}

	if pos+3 > len(stack) {
		resized := make([]Poly, (pos+3)*2)
		copy(resized, stack)

		stack = resized
	}

	stack[pos+1] = Poly{}
	stack[pos+2] = Poly{}

	for i := 0; ; i, depth = i+1, depth+1 {
		stack[pos+2] = Trace(a, stack[pos])

		if i >= 1 && !factorizable {
			stack[pos+1] = stack[pos+2].Sqr()

			for i := 0; i < len(stack[pos+2]); i++ {
				stack[pos+1][i] = stack[pos+1][i].Add(stack[pos+2][i])
			}

			for len(stack[pos+1]) > 0 && stack[pos+1][len(stack[pos+1])-1] == 0 {
				stack[pos+1] = stack[pos+1][:len(stack[pos+1])-1]
			}

			stack[pos+1] = stack[pos+1].Mod(stack[pos])

			if len(stack[pos+1]) != 0 {
				return nil
			}

			factorizable = true
		}

		if factorizable && !((len(stack[pos])-2)>>uint(Dimension-depth) == 0) {
			return nil
		}

		a = a.Mul(2)

		stack[pos+1] = make(Poly, len(stack[pos]))
		copy(stack[pos+1], stack[pos])

		stack[pos+2] = GCD(stack[pos+2], stack[pos+1])
		stack[pos+1] = Poly{}

		if len(stack[pos+2]) != len(stack[pos]) && len(stack[pos+2]) > 1 {
			break
		}
	}

	stack[pos+2] = stack[pos+2].Monic()

	stack[pos+1], _ = stack[pos].Div(stack[pos+2])
	stack[pos], stack[pos+2] = stack[pos+2], stack[pos]

	roots = append(roots, findRoots(stack, pos+1, factorizable, depth, a)...)
	roots = append(roots, findRoots(stack, pos, true, depth, a)...)

	return roots
}
