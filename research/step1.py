# A demonstration of how PinSketch works for set reconciliation
# by Kenta Iwasaki and Briney Bai.

from functools import reduce

import numpy as np

np.set_printoptions(precision=2)


def S1(M, x, n=10):  # S(M) = sum[ sum[m_i ^ i] * x ^ i ]
    return sum((sum(m_i ** i for m_i in M) ** i) * (x ** i) for i in range(n))


def S2(M, x):  # S(M) = product[ (1 - m_i * x) ^ (-1) ]
    return reduce(lambda result, m_i: result * ((1 - m_i * x) ** -1), M, 1.0)


# If lim |mx| -> 0, prove lim S(m) -> 1.
def check_formal_power_series_convergence(M):
    for x in range(100):
        x = 1 / (x + 1)
        print("Given x = %f, S1(M, x) = %f and S2(M, x) = %f." % (x, S1(M, x), S2(M, x)))
        print("Dividing S1(M, x) / S2(M, x) yields %f." % (S1(M, x) / S2(M, x)))
        print()


if __name__ == "__main__":
    check_formal_power_series_convergence([0.25, 0.5, 0.75])
