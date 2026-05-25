import itertools

# octMultTable from SevenDAlgebra.lean
octMultTable = [
  [(1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (1,7)],   # 0 *
  [(1,1), (-1,0), (1,3), (-1,2), (1,5), (-1,4), (1,7), (-1,6)], # 1 *
  [(1,2), (-1,3), (-1,0), (1,1), (1,6), (1,7), (-1,4), (-1,5)],
  [(1,3), (1,2), (-1,1), (-1,0), (1,7), (-1,6), (1,5), (-1,4)],
  [(1,4), (-1,5), (-1,6), (-1,7), (-1,0), (1,1), (1,2), (1,3)],
  [(1,5), (1,4), (-1,7), (1,6), (-1,1), (-1,0), (-1,3), (1,2)],
  [(1,6), (1,7), (1,4), (-1,5), (-1,2), (1,3), (-1,0), (-1,1)],
  [(1,7), (-1,6), (1,5), (1,4), (-1,3), (-1,2), (1,1), (-1,0)]
]

def gamma(k):
    # k in 0..6
    M = [[0]*8 for _ in range(8)]
    row = octMultTable[k+1]
    for p in range(8):
        for m in range(8):
            sgn, target = row[m]
            if target == p:
                M[p][m] = sgn
    return M

def mat_mul(A, B):
    C = [[0]*8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            s = 0
            for k in range(8):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C

def sector_diags(M, inds):
    return [M[i][i] for i in inds]

L_bivectors = []
for i in range(7):
    for j in range(i+1, 7):
        L_bivectors.append(mat_mul(gamma(i), gamma(j)))

F_fourforms = []
for i, j, k, l in itertools.combinations(range(7), 4):
    M = mat_mul(mat_mul(gamma(i), gamma(j)), mat_mul(gamma(k), gamma(l)))
    F_fourforms.append(M)

d_inds = [1,2,3]
u_inds = [4,5,6]
lep_inds = [0,7]

def analyze(name, ops, inds):
    print(f"\n--- {name} ---")
    diags = [sector_diags(M, inds) for M in ops]
    zeros = sum(1 for d in diags if d == [0]*len(inds))
    nonzeros = [d for d in diags if d != [0]*len(inds)]
    print(f"Total ops: {len(ops)}")
    print(f"All zero diags: {zeros}")
    if nonzeros:
        print(f"Nonzero diags patterns (unique):")
        uniq = set(tuple(d) for d in nonzeros)
        for u in uniq: print(u)

analyze("L_bivectors on Leptons (N=0,3)", L_bivectors, lep_inds)
analyze("L_bivectors on d-quarks (N=1)", L_bivectors, d_inds)
analyze("F_fourforms on Leptons", F_fourforms, lep_inds)
analyze("F_fourforms on d-quarks", F_fourforms, d_inds)
analyze("F_fourforms on u-quarks", F_fourforms, u_inds)
analyze("L_bivectors on u-quarks", L_bivectors, u_inds)
