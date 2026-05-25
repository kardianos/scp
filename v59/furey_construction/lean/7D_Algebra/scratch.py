import sympy as sp

# 8x8 octonion mult table
# indices: 0 is scalar, 1-7 are e1-e7
table = [
  [(1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (1,7)],
  [(1,1), (-1,0), (1,3), (-1,2), (1,5), (-1,4), (1,7), (-1,6)],
  [(1,2), (-1,3), (-1,0), (1,1), (1,6), (1,7), (-1,4), (-1,5)],
  [(1,3), (1,2), (-1,1), (-1,0), (1,7), (-1,6), (1,5), (-1,4)],
  [(1,4), (-1,5), (-1,6), (-1,7), (-1,0), (1,1), (1,2), (1,3)],
  [(1,5), (1,4), (-1,7), (1,6), (-1,1), (-1,0), (-1,3), (1,2)],
  [(1,6), (1,7), (1,4), (-1,5), (-1,2), (1,3), (-1,0), (-1,1)],
  [(1,7), (-1,6), (1,5), (1,4), (-1,3), (-1,2), (1,1), (-1,0)]
]

c = sp.symbols('c0:8')

def mult(a, b):
    res = [0]*8
    for i in range(8):
        for j in range(8):
            sgn, k = table[i][j]
            res[k] += sgn * a[i] * b[j]
    return res

def conj(a):
    return [a[0]] + [-x for x in a[1:]]

def scalar(a):
    return a[0]

def norm_sq(a):
    return sum(x**2 for x in a)

c_var = list(c)

# Potential: V = 0.5 * (norm_sq(M) - v2)^2 + lambda * norm_sq(mult(M, M)) + mu * norm_sq(M)?
# Wait, the note says: V = rho_M + lambda * scalar(M * M) + mu * norm2
# But rho_M is 0.5 * (norm_sq(M) - v2). Oh, wait, in OctonionAlgebra.lean:
# rhoM(c, v2) = 0.5 * (octScalarProd c (octConj c) - v2)
# Wait, if V = rho_M + ... that is linear in norm2. That doesn't have a minimum at v2.
# Usually it's (rho_M)^2. Let's check max/octonion_sensitivity_analysis.mac or full_octonion_perturbation.mac.
