import itertools

octMultTable = [
  [(1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (1,7)],
  [(1,1), (-1,0), (1,3), (-1,2), (1,5), (-1,4), (1,7), (-1,6)],
  [(1,2), (-1,3), (-1,0), (1,1), (1,6), (1,7), (-1,4), (-1,5)],
  [(1,3), (1,2), (-1,1), (-1,0), (1,7), (-1,6), (1,5), (-1,4)],
  [(1,4), (-1,5), (-1,6), (-1,7), (-1,0), (1,1), (1,2), (1,3)],
  [(1,5), (1,4), (-1,7), (1,6), (-1,1), (-1,0), (-1,3), (1,2)],
  [(1,6), (1,7), (1,4), (-1,5), (-1,2), (1,3), (-1,0), (-1,1)],
  [(1,7), (-1,6), (1,5), (1,4), (-1,3), (-1,2), (1,1), (-1,0)]
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

def analyze_squares(name, ops):
    print(f"--- {name} ---")
    squares = [mat_mul(M, M) for M in ops]
    uniq_squares = []
    for sq in squares:
        if sq not in uniq_squares:
            uniq_squares.append(sq)
    
    for sq in uniq_squares:
        diag = [sq[i][i] for i in range(8)]
        off_diags = sum(abs(sq[i][j]) for i in range(8) for j in range(8) if i!=j)
        print(f"Diag: {diag}, sum(|off_diags|): {off_diags}")

analyze_squares("L_bivectors squares", L_bivectors)
analyze_squares("F_fourforms squares", F_fourforms)
