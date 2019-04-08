# A demonstration of how PinSketch works for set reconciliation
# by Kenta Iwasaki and Briney Bai.

import numpy as np

np.set_printoptions(precision=2)


# M is our set of items that we wish to convert into a compact sketch.
#
# A sketch is simply the first 2n coefficients of the formal power
# series definition of S(M), where n is the capacity (maximum amount
# of items) we are able to hold in a sketch.
#
# The formal power series definition of S(M) is: sum[ sum[m_i ^ i] * x ^ i ]
#
# For clarity, this means that our set of items can be represented
# through simply just 2n coefficient terms of S(M).
#
# To decode our sketch back into our set of items, we take advantage
# of the fact that S(M) has a polynomial representation that can be
# factorized down to our items, M.
#
# That is, S(M) also equals: product[ (1 - m_i * x) ^ (-1) ]
#
# So, we form n linear equations with n terms, solving for the coefficients
# of S(M)'s polynomial representation given the 2n coefficients of the
# formal power series definition of S(M).
#
# More explicitly, let L represent the coefficients of S(M)'s polynomial
# representation, and S represent the coefficients of S(M)'s formal power
# series representation.
#
# Solve for L such that S * L = 0.

def create_sketch_naive(M, capacity):
    return [sum(m_i ** (i + 1) for m_i in M) for i in range(2 * capacity)]


if __name__ == "__main__":
    M = [5000, 3000, 2000]
    n = 3

    print("These are our set of items M:", M, "\n")

    sketch = create_sketch_naive(M, capacity=n)
    print(
        "Our sketch, with a capacity of [n = %d] items, is the coefficients of the first [2n = %d] terms of the formal "
        "power series representation of S(M):" % (n, 2 * n), sketch, "\n"
    )

    S_m = np.zeros((n, n + 1))
    for row in range(n):
        S_m[row] = list(reversed(sketch[row: row + n + 1]))

    print(
        "We solve for the coefficients of the polynomial representation of S(M), by solving for L such that S(M)L = "
        "0.\n\nS(M)L:\n%s\n" % S_m
    )

    L = [1.0] + list(np.linalg.solve(S_m[:, 1:], -1 * S_m[:, 0]))
    L = list(reversed(L))

    print("L's polynomial coefficients in ascending order after solving for L:", L)

    print("L's polynomial factors through efficient root inversion (recovered items):", np.roots(list(reversed(L))))
    print("L's polynomial factors through manual root inversion (recovered items):", 1 / np.roots(L), "\n")

    enough_capacity = np.linalg.det(S_m[:, 1:]) > 0
    print("Could we successfully decode the sketch with a capacity of %d items, given %d items?" % (n, len(M)),
          "Yes." if enough_capacity else "No.")
