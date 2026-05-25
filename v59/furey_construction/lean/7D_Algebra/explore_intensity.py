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

def sector_trace(M, inds):
    return sum(M[i][i] for i in inds)

L_bivectors = []
for i in range(7):
    for j in range(i+1, 7):
        L_bivectors.append(mat_mul(gamma(i), gamma(j)))

F_fourforms = []
for i, j, k, l in itertools.combinations(range(7), 4):
    M = mat_mul(mat_mul(gamma(i), gamma(j)), mat_mul(gamma(k), gamma(l)))
    F_fourforms.append(M)

inds_lep = [0, 7]
inds_d = [1, 2, 3]
inds_u = [4, 5, 6]

def analyze_sector_intensity(name, ops):
    print(f"\n=== {name} ===")
    total_lep = sum(sector_trace(mat_mul(M, M), inds_lep) for M in ops)
    total_d = sum(sector_trace(mat_mul(M, M), inds_d) for M in ops)
    total_u = sum(sector_trace(mat_mul(M, M), inds_u) for M in ops)
    
    print(f"Total Intensity (Sum of M^2 traces) on Leptons (N=0,3): {total_lep}")
    print(f"Total Intensity (Sum of M^2 traces) on d-quarks (N=1): {total_d}")
    print(f"Total Intensity (Sum of M^2 traces) on u-quarks (N=2): {total_u}")

analyze_sector_intensity("L-grade (21 bivectors)", L_bivectors)
analyze_sector_intensity("F-grade (35 fourforms)", F_fourforms)
