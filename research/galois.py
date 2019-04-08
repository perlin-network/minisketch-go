# A demonstration of how PinSketch works for set reconciliation
# by Kenta Iwasaki and Briney Bai.

def gadd(a, b):
    return a ^ b


def gsub(a, b):
    return a ^ b


# x^64 + x^4 + x^3 + x + 1

def gmul(a, b, bits=64, modulus=int('10000000000000000000000000000000000000000000000000000000000011011', 2)):
    p = 0

    for i in range(bits):
        p ^= -(b & 1) & a

        mask = -((a >> bits - 1) & 1)
        a = (a << 1) ^ (modulus & mask)
        b >>= 1

    return p


def ginv(a, bits=64, modulus=int('10000000000000000000000000000000000000000000000000000000000011011', 2)):
    if a == 0: return a

    order, i, p = 1 << bits, 2, 1

    while i < order:
        a = gmul(a, a, bits=bits, modulus=modulus)
        p = gmul(p, a, bits=bits, modulus=modulus)
        i <<= 1

    return p


def gdiv(a, b, bits=64, modulus=int('10000000000000000000000000000000000000000000000000000000000011011', 2)):
    return gmul(a, ginv(b, bits=bits, modulus=modulus), bits=bits, modulus=modulus)


def gexp(a, e, bits=64, modulus=int('10000000000000000000000000000000000000000000000000000000000011011', 2)):
    accum = 1
    while e != 0:
        if e & 0x01 == 1: accum = gmul(accum, a, bits=bits, modulus=modulus)
        a = gmul(a, a, bits=bits, modulus=modulus)
        e >>= 1
    return accum


def gsqrt(a, bits=64, modulus=int('10000000000000000000000000000000000000000000000000000000000011011', 2)):
    for i in range(bits - 1):
        a = gexp(a, 2, bits=bits, modulus=modulus)
    return a


if __name__ == "__main__":
    print("Modulus:", int('10000000000000000000000000000000000000000000000000000000000011011', 2))
    print("2^64 - 1:", 2 ** 64 - 1)
    print("(2^64 - 1)(2^64 - 1):", gmul(2 ** 64 - 1, 2 ** 64 - 1))
    print("(2^64 - 1)^(-1):", ginv(gmul(2 ** 64 - 1, 2 ** 64 - 1)))
    print("((2^64 - 1)(2^64 - 1)) / (2^64 - 1):", gdiv(gmul(2 ** 64 - 1, 2 ** 64 - 1), 2 ** 64 - 1))
    print("[(2^64 - 1)(2^64 - 1)] ^ 12:", gexp(gmul(2 ** 64 - 1, 2 ** 64 - 1), 12))
    print("sqrt((2^64 - 1)(2^64 - 1)):", gsqrt(gmul(2 ** 64 - 1, 2 ** 64 - 1)))
