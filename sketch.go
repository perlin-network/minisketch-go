package minisketch

// Sketch may be used for deriving the symmetric set difference of an arbitrary
// set of items, by being merged with a sketch from some external third-party.
//
// It takes advantage of the isomorphism between the formal power series, and
// its infinite sum counterpart representative as an irreducible polynomial to
// yield a set representation.
//
// Numerous optimizations could be done by considering the isomorphism taken
// under a finite field, such as GF(2^64), such as taking advantage of the
// Frobenius endomorphism to only require sending every odd term coefficient
// of the formal power series before sending our sketch to a third party, as
// our external third-party can expand every odd term coefficient into its
// even term coefficient counterpart through the "freshman's dream".
type Sketch []Element

// Add immutably adds elements to the sketch. If an element is added
// which was already added to the sketch, then the element is deleted
// from the sketch.
//
func (s Sketch) Add(items ...Element) Sketch {
	for _, item := range items {
		sqr := item.Mul(item)

		for i := range s {
			s[i] = s[i].Add(item)
			item = item.Mul(sqr)
		}
	}

	return s
}

// Merge immutably merges two sketches together, returning a new sketch
// whose capacity is the minimum capacity of the two sketches.
//
// Merging two sketches allows for one to then decode for the symmetric
// set difference between the items representative of said two sketches.
func (s Sketch) Merge(other Sketch) Sketch {
	if len(other) < len(s) {
		s = s[:len(other)]
	}

	for i := range s {
		s[i] = s[i].Add(other[i])
	}

	return s
}

// Decode decodes a sketch which has less items than its predetermined
// capacity into its list of items.
func (s Sketch) Decode() (items []Element) {
	terms := Frobenius(s)
	polynomial := BerlekampMassey(terms)

	if len(polynomial) == 0 || len(polynomial) == 1 {
		return []Element{}
	}

	for left, right := 0, len(polynomial)-1; left < right; left, right = left+1, right-1 {
		polynomial[left], polynomial[right] = polynomial[right], polynomial[left]
	}

	return FindRoots(polynomial)
}

// New creates a new sketch with a specified capacity.
func New(capacity int) Sketch {
	return make(Sketch, capacity)
}

// BerlekampMassey yields a minimal polynomial of a recurrent sequence in
// any arbitrary finite field by its coefficients in descending power order.
//
// It is used in particular to solve for the error locator polynomial in
// Bose–Chaudhuri–Hocquenghem (BCH) codes.
func BerlekampMassey(syndromes []Element) []Element {
	current, prev, tmp := Poly{1}, Poly{1}, Poly{}
	b, bInv, bHaveInv := Element(1), Element(1), true

	for n := 0; n < len(syndromes); n++ {
		d := syndromes[n]

		for i := 1; i < len(current); i++ {
			d = d.Add(current[i].Mul(syndromes[n-i]))
		}

		if d == 0 {
			continue
		}

		x := n + 1 - (len(current) - 1) - (len(prev) - 1)

		if !bHaveInv {
			bInv, bHaveInv = b.Inv(), true
		}

		swap := 2*(len(current)-1) <= n

		if swap {
			resized := make(Poly, len(prev)+x)
			copy(resized, current)

			current, tmp = resized, current
		}

		mul := d.Mul(bInv)

		for i := 0; i < len(prev); i++ {
			current[i+x] = current[i+x].Add(prev[i].Mul(mul))
		}

		if swap {
			prev, tmp = tmp, prev
			b, bHaveInv = d, false
		}
	}

	return current
}

// Frobenius takes advantage of the Frobenius endomorphism of finite fields
// to expand a set of odd syndromes into its even and odd counterparts.
//
// It is used for de-compacting a sketch, which is typically denoted by its
// odd coefficients under a finite field representative of a formal power
// series' terms.
func Frobenius(syndromes []Element) []Element {
	result := make([]Element, len(syndromes)*2)

	for i := 0; i < 2*len(syndromes); i++ {
		if i&1 == 1 {
			result[i] = result[i>>1].Exp(2)
		} else {
			result[i] = syndromes[i>>1]
		}
	}

	return result
}
