# A demonstration of how PinSketch works for set reconciliation
# by Kenta Iwasaki and Briney Bai.

import numpy as np

from galois import *

np.set_printoptions(precision=2)


def create_sketch_naive(M, capacity):
    sketch = []

    for i in range(2 * capacity):
        term = 0

        for m_i in M:
            term = gadd(term, gexp(m_i, (i + 1)))

        sketch.append(term)

    return sketch


def create_sketch(M, capacity):
    sketch = []

    for i in range(0, 2 * capacity, 2):
        term = 0

        for m_i in M:
            term = gadd(term, gexp(m_i, (i + 1)))

        sketch.append(term)

    return sketch


def frobenius(sketch):
    result = [0] * len(sketch) * 2

    for i in range(2 * len(sketch)):
        result[i] = gexp(result[i >> 1], 2) if i & 1 else sketch[i >> 1]

    return result


def berlekamp_massey(S):
    current, prev, tmp = [1], [1], []
    b, b_inv, b_have_inv = 1, 1, True

    for n in range(len(S)):
        d = S[n]

        for i in range(1, len(current), 1):
            d = gadd(d, gmul(S[n - i], current[i]))

        if d == 0: continue

        x = n + 1 - (len(current) - 1) - (len(prev) - 1)

        if not b_have_inv:
            b_inv, b_have_inv = ginv(b), True

        swap = 2 * (len(current) - 1) <= n

        if swap:
            tmp = current
            current = current + [0] * (len(prev) + x - len(current))

        mul = gmul(d, b_inv)
        for i in range(len(prev)):
            current[i + x] = gadd(current[i + x], gmul(prev[i], mul))

        if swap:
            prev, tmp = tmp, prev
            b, b_have_inv = d, False
    return current


def monic(a):
    if a[-1] == 1: return a

    inv = ginv(a[-1])
    a[-1] = 1

    for i in range(len(a) - 1):
        a[i] = gmul(inv, a[i])

    return a


def poly_mod(val, mod):
    if len(mod) > 0 and mod[-1] != 1: raise Exception("modulo must be a monic polynomial")
    if len(val) < len(mod): return val
    if val[-1] == 0: raise Exception("val must not have a leading coefficient of zero")

    while len(val) >= len(mod):
        term = val.pop(-1)

        if term == 0: continue

        for x in range(len(mod) - 1):
            val[len(val) - len(mod) + 1 + x] = gadd(val[len(val) - len(mod) + 1 + x], gmul(term, mod[x]))

    while len(val) > 0 and val[-1] == 0:
        val.pop(-1)

    return val


def div_mod(val, mod, div):
    if len(val) < len(mod):
        div.clear()
        return div, val

    div = div + [0] * ((len(val) - len(mod) + 1) - len(div))

    while len(val) >= len(mod):
        term = val.pop(-1)

        div[len(val) + 1 - len(mod)] = term

        if term == 0: continue

        for x in range(len(mod) - 1):
            val[len(val) - len(mod) + 1 + x] = gadd(val[len(val) - len(mod) + 1 + x], gmul(mod[x], term))

    while len(val) > 0 and val[-1] == 0:
        val.pop(-1)

    return div, val


def trace_mod(param, mod, bits=64):
    out = [0, param]

    for i in range(bits - 1):
        out = poly_sqr(out)

        if len(out) < 2:
            out = out + [0] * (2 - len(out))

        out[1] = param
        poly_mod(out, mod)

    return out


def poly_sqr(a):
    if len(a) == 0: return a
    a = a + [0] * ((len(a) * 2 - 1) - len(a))

    for x in range(len(a) - 1, -1, -1):
        a[x] = 0 if (x & 1) else gexp(a[x // 2], 2)

    return a


def gcd(a, b):
    if len(a) < len(b):
        a, b = b, a

    while len(b) > 0:
        if len(b) == 1:
            return [1], b

        b = monic(b)
        a = poly_mod(a, b)

        a, b = b, a

    return a, b


def find_roots(poly):
    roots = []
    stack = [poly]

    if not rec_find_roots(stack, 0, roots, False, 0, 1):
        return []

    return roots


def rec_find_roots(stack, pos, roots, fully_factorizable, depth, randv):
    if len(stack[pos]) == 2:
        roots.append(stack[pos][0])
        return True

    if pos + 3 > len(stack):
        stack = stack + ([[]] * (((pos + 3) * 2) - len(stack)))

    iter = 0

    stack[pos + 1].clear()
    stack[pos + 2].clear()

    while True:
        stack[pos + 2] = trace_mod(randv, stack[pos])

        if iter >= 1 and not fully_factorizable:
            stack[pos + 1] = stack[pos + 2].copy()
            poly_sqr(stack[pos + 1])

            for i in range(len(stack[pos + 2])):
                stack[pos + 1][i] = gadd(stack[pos + 1][i], stack[pos + 2][i])

            while len(stack[pos + 1]) > 0 and stack[pos + 1][-1] == 0: stack[pos + 1].pop(-1)

            poly_mod(stack[pos + 1], stack[pos])

            if len(stack[pos + 1]) != 0:
                return False

            fully_factorizable = True

        if fully_factorizable:
            if not ((len(stack[pos]) - 2) >> (64 - depth) == 0):
                return False

        depth += 1

        randv = gmul(randv, 2)
        stack[pos + 1] = stack[pos].copy()

        gcd(stack[pos + 2], stack[pos + 1])

        if len(stack[pos + 2]) != len(stack[pos]) and len(stack[pos + 2]) > 1:
            break

        iter += 1

    monic(stack[pos + 2])
    stack[pos + 1], _ = div_mod(stack[pos], stack[pos + 2], stack[pos + 1])

    stack[pos], stack[pos + 2] = stack[pos + 2], stack[pos]

    if not rec_find_roots(stack, pos + 1, roots, fully_factorizable, depth, randv):
        return False

    return rec_find_roots(stack, pos, roots, True, depth, randv)


if __name__ == "__main__":
    M = [5000, 3000, 2000]
    n = 3

    print("Here are all [2 * n = %d] coefficients representing our sketch:" % (2 * n), create_sketch_naive(M, n))
    print(
        "We actually only need to send every odd coefficient ([n = %d] coefficients) to represent our sketch as "
        "we can take advantage of the Frobenius endomorphism as our coefficients are elements under a Galois field:" % n,
        create_sketch(M, n)
    )

    sketch = frobenius(create_sketch(M, n))

    print(
        "With the Frobenius endomorphism, we recovered [2 * n = %d] coefficients:" % (2 * n),
        sketch
    )

    L = list(reversed(berlekamp_massey(sketch)))

    print("L's polynomial coefficients in ascending order after solving for L:", L)
    print("L's polynomial factors through root inversion via Berlekamp's Trace Algorithm (recovered items):",
          find_roots(L))
