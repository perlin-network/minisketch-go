# A demonstration of how PinSketch works for set reconciliation
# by Kenta Iwasaki and Briney Bai.

from step3 import *

A = [2000, 3000, 5000]
B = [4000, 5000, 1000]
n = 12

A_sketch = create_sketch(A, n)
B_sketch = create_sketch(B, n)

C_sketch = frobenius([gadd(b, a) for a, b in zip(A_sketch, B_sketch)])
L = list(reversed(berlekamp_massey(C_sketch)))

print(find_roots(L))
