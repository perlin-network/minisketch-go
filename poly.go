package minisketch

// Poly is a polynomial whose coefficienets are under the field GF(2^64).
type Poly []Element

// Monic converts an arbitrary polynomial into a monic polynomial.
func (val Poly) Monic() Poly {
	end := len(val) - 1

	if val[end] == 1 {
		return val
	}

	inv := val[end].Inv()
	val[end] = 1

	for i := 0; i < end; i++ {
		val[i] = val[i].Mul(inv)
	}

	return val
}

// Mod yields the remainder after dividing a polynomial by another polynomial.
func (val Poly) Mod(mod Poly) Poly {
	if len(val) < len(mod) {
		return val
	}

	for len(val) >= len(mod) {
		term := val[len(val)-1]
		val = val[:len(val)-1]

		if term == 0 {
			continue
		}

		for i := 0; i < len(mod)-1; i++ {
			val[len(val)-len(mod)+1+i] = val[len(val)-len(mod)+1+i].Add(term.Mul(mod[i]))
		}
	}

	for len(val) > 0 && val[len(val)-1] == 0 {
		val = val[:len(val)-1]
	}

	return val
}

// Div yields both the quotient and remainder after diving a polynomial by
// another polynomial.
func (val Poly) Div(mod Poly) (Poly, Poly) {
	if len(val) < len(mod) {
		return Poly{}, val
	}

	div := make(Poly, len(val)-len(mod)+1)

	for len(val) >= len(mod) {
		term := val[len(val)-1]
		val = val[:len(val)-1]

		div[len(val)-len(mod)+1] = term

		if term == 0 {
			continue
		}

		for i := 0; i < len(mod)-1; i++ {
			val[len(val)-len(mod)+1+i] = val[len(val)-len(mod)+1+i].Add(term.Mul(mod[i]))
		}
	}

	for len(val) > 0 && val[len(val)-1] == 0 {
		val = val[:len(val)-1]
	}

	return div, val
}

// Sqr multiplies a polynomial with itself, or otherwise computes the
// polynomial raised to the power of two.
func (val Poly) Sqr() Poly {
	if len(val) == 0 {
		return val
	}

	out := make(Poly, len(val)*2-1)
	copy(out, val)

	for x := len(out) - 1; x >= 0; x-- {
		if x&1 == 1 {
			out[x] = 0
		} else {
			out[x] = out[x/2].Exp(2)
		}
	}

	return out
}

// Trace computes a polynomial whose coefficients are over GF(2^64) divisible by
// any other polynomial whose coefficients are over GF(2^64).
//
// It is used for algorithms involved with factoring the roots of a polynomial over
// a finite field, such as Berlekamp's Trace algorithm.
func Trace(a Element, mod Poly) Poly {
	out := Poly{0, a}

	for i := 0; i < Dimension-1; i++ {
		out = out.Sqr()

		if len(out) < 2 {
			resized := make(Poly, 2)
			copy(resized, out)
			out = resized
		}

		out[1] = a
		out = out.Mod(mod)
	}

	return out
}

// GCD returns the greatest common divisor polynomial between two polynomials.
func GCD(a, b Poly) Poly {
	if len(a) < len(b) {
		a, b = b, a
	}

	for len(b) > 0 {
		if len(b) == 1 {
			return Poly{1}
		}

		b = b.Monic()
		a = a.Mod(b)

		a, b = b, a
	}

	return a
}
