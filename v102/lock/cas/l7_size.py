#!/usr/bin/env python3
"""L7 — size, a report (../PREREG.md §4 L7, as amended): E0(n) of S-law
rings, n = 3..10, S^z = 0, g/Delta = 0.1; n <= 8 without truncation,
n >= 9 with sum c <= 3."""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import DELTA, basis_ring, hamiltonian

t0 = time.time()
g = 0.1*DELTA
E = {}
print(f"# L7 size: E0(n) of S-law rings, g/Delta = 0.1 (energies relative to n E_f(2))")
for n in range(3, 11):
    cmax = None if n <= 8 else 3
    st, ix = basis_ring(n, cmax=cmax, sz=0)
    H = hamiltonian(n, g, st, ix)
    w = np.linalg.eigvalsh(H.toarray())[:2] if H.shape[0] < 3000 else np.sort(eigsh(H, k=2, which='SA', tol=1e-12)[0])
    E[n] = w[0]
    print(f"  n={n:2d} cmax={str(cmax):>4} dim={H.shape[0]:7d}  E0 = {w[0]:+.7f}  E0/n = {w[0]/n:+.7f}  [{time.time()-t0:.0f}s]", flush=True)
st, ix = basis_ring(8, cmax=3, sz=0)
e83 = np.sort(eigsh(hamiltonian(8, g, st, ix), k=1, which='SA', tol=1e-12)[0])[0]
print(f"  truncation scale at n=8: E0(sum c<=3) - E0(exact) = {e83 - E[8]:+.2e}")
print("  second difference E0(n+1) + E0(n-1) - 2 E0(n):")
print("   " + "  ".join(f"n={n}:{E[n+1]+E[n-1]-2*E[n]:+.2e}" for n in range(4, 10)))
print("  fission energy E0(a) + E0(b) - E0(a+b)  (> 0: the joined ring is lower):")
for N in range(6, 11):
    print(f"   N={N:2d}: " + "  ".join(f"{a}+{N-a}:{E[a]+E[N-a]-E[N]:+.5f}" for a in range(3, N//2 + 1)))
print(f"  E0/n is {'monotone decreasing' if all(E[n+1]/(n+1) < E[n]/n for n in range(3,10)) else 'not monotone'} over n = 3..10; "
      f"even n: " + " ".join(f"{E[n]/n:+.6f}" for n in (4, 6, 8, 10)) + "; odd n: " + " ".join(f"{E[n]/n:+.6f}" for n in (3, 5, 7, 9)))
