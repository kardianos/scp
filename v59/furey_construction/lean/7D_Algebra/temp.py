import itertools

octMultTable = [
  [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7)],
  [(1, 1), (-1, 0), (1, 3), (-1, 2), (1, 5), (-1, 4), (-1, 7), (1, 6)],
  [(1, 2), (-1, 3), (-1, 0), (1, 1), (1, 6), (1, 7), (-1, 4), (-1, 5)],
  [(1, 3), (1, 2), (-1, 1), (-1, 0), (1, 7), (-1, 6), (-1, 5), (1, 4)],
  [(1, 4), (-1, 5), (-1, 6), (-1, 7), (-1, 0), (1, 1), (1, 2), (1, 3)],
  [(1, 5), (1, 4), (-1, 7), (1, 6), (-1, 1), (-1, 0), (1, 3), (-1, 2)],
  [(1, 6), (1, 7), (1, 4), (1, 5), (-1, 2), (-1, 3), (-1, 0), (-1, 1)],
  [(1, 7), (-1, 6), (1, 5), (-1, 4), (-1, 3), (1, 2), (1, 1), (-1, 0)]
]

def gamma(k):
    M = [[0]*8 for _ in range(8)]
    row = octMultTable[k+1]
    for p in range(8):
        for m in range(8):
            sgn, target = row[m]
            if target == p: M[p][m] = sgn
    return M

def mat_mul(A, B):
    C = [[0]*8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            s = sum(A[i][k] * B[k][j] for k in range(8))
            C[i][j] = s
    return C

L_bivectors = []
for i in range(7):
    for j in range(i+1, 7):
        L_bivectors.append(mat_mul(gamma(i), gamma(j)))

F_fourforms = []
for i, j, k, l in itertools.combinations(range(7), 4):
    M = mat_mul(mat_mul(gamma(i), gamma(j)), mat_mul(gamma(k), gamma(l)))
    F_fourforms.append(M)

print("--- L-grade per index ---")
l_res = [sum(mat_mul(M,M)[i][i] for M in L_bivectors) for i in range(8)]
print(l_res)

print("--- F-grade per index ---")
f_res = [sum(mat_mul(M,M)[i][i] for M in F_fourforms) for i in range(8)]
print(f_res)
