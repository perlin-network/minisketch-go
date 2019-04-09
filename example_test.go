package minisketch

import "fmt"

// Example is example usage on discovering the symmetric set difference of items
// between two peers Alice and Bob, given that they produce sketches of the set
// of items they have without going over their set capacity limits.
func Example() {
	// The list of items Alice and Bob have which they are looking to reconcile.
	a := []Element{2000, 4000, 5000}
	b := []Element{4000, 5000, 1000}

	// The sketches Alice and Bob have produced of their items, with sufficient
	// capacity len(a) + len(b).
	alice := New(len(a) + len(b)).Add(a...)
	bob := New(len(a) + len(b)).Add(b...)

	// Bob merges his sketch with Alice, and decodes the symmetric set difference
	// of his items versus Alice's  items.
	merged := bob.Merge(alice).Decode()

	// Print elements which are not in both sets.
	for _, elem := range merged {
		fmt.Println(elem)
	}

	// Output:
	// 2000
	// 1000
}
