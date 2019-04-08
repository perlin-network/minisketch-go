package minisketch

// GF(2^64)

const (
	// The number of times to add an element together to yield a sum of zero.
	Characteristic = 2

	// The number of dimensions of each element in the finite field.
	Dimension = 64

	// The total number of elements in the finite field.
	Order = 1 << Dimension

	// The generator modulus polynomial defining the finite field.
	// By default, the modulus specified below is x^64 + x^4 + x^3 + x + 1.
	Modulus = Element(1<<4 | 1<<3 | 1<<1 | 1<<0)
)

// Element is an element under GF(2^64).
type Element uint64

// Add adds two elements.
func (a Element) Add(b Element) Element {
	return a ^ b
}

// Sub subs two elements.
func (a Element) Sub(b Element) Element {
	return a ^ b
}

// Mul multiplies two elements.
func (a Element) Mul(b Element) (p Element) {
	for i := 0; i < Dimension; i++ {
		ab0 := a &^ (b&1 - 1)
		ra63 := Modulus &^ (a>>(Dimension-1) - 1)

		p, a, b = p^ab0, a<<1^ra63, b>>1
	}

	return
}

// Inv provides the multiplicative inverse of an element.
func (a Element) Inv() Element {
	if a == 0 {
		return a
	}

	p := Element(1)

	for i := 0; i < Dimension-1; i++ {
		a = a.Mul(a)
		p = p.Mul(a)
	}

	return p
}

// Div multiplies an element by another elements inverse. It is
// equivalent to dividing an element by another element.
func (a Element) Div(b Element) Element {
	return a.Mul(b.Inv())
}

// Exp raises an element up to a power specified by another element.
func (a Element) Exp(b Element) (p Element) {
	for p = 1; b != 0; b >>= 1 {
		if b&0x01 == 1 {
			p = p.Mul(a)
		}
		a = a.Mul(a)
	}

	return
}

// Sqrt provides the square root of an element, such that squaring
// the yielded element returns the original element.
func (a Element) Sqrt() Element {
	for i := 0; i < Dimension-1; i++ {
		a = a.Exp(2)
	}

	return a
}
