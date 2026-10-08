#!/usr/bin/env python3
"""L4(b) — one line (../PREREG.md §3 L4b).  Micro ring n = 8, sum c <= 3,
g/Delta = 0.1: lowest S^z = 2 state and lowest non-triplet S = 0
excitation against 2 Delta_8 (bar 0.98); no level with <sum c> > 0.5
below 2 Delta_8.  Reported: the first space-carrying level on the
n = 6 ring without truncation (dense)."""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import DELTA, basis_ring, hamiltonian, diag_obs

t0 = time.time()
g = 0.1*DELTA; n = 8; cmax = 3
lv = {}
for (sz, k) in ((0, 40), (2, 20), (4, 4)):
    st, ix = basis_ring(n, cmax=cmax, sz=sz)
    H = hamiltonian(n, g, st, ix)
    w, v = eigsh(H, k=k, which='SA', tol=1e-12); o = np.argsort(w); w = w[o]; v = v[:, o]
    C, _, _ = diag_obs(n, st)
    sc = (v**2*C.sum(1)[:, None]).sum(0)
    lv[sz] = (w, sc)
E0 = lv[0][0][0]; D8 = lv[2][0][0] - E0
e2 = lv[4][0][0] - E0
trip = lv[2][0]
s0 = [x for x in lv[0][0][1:] if np.min(np.abs(trip - x)) > 1e-9]
es0 = s0[0] - E0
below = lv[0][0] - E0 < 2*D8
spc = lv[0][1][below].max()
print(f"# L4(b) one line, n = {n}, sum c <= {cmax}, g/Delta = 0.1")
print(f"  gap Delta_8 = {D8:.6f};  lowest S^z=2: {e2:.6f} = {e2/(2*D8):.4f} x 2 Delta_8;  "
      f"lowest non-triplet S=0 excitation: {es0:.6f} = {es0/(2*D8):.4f} x 2 Delta_8")
print(f"  levels below 2 Delta_8 in S^z = 0: {below.sum()};  their max <sum c> = {spc:.4f} (bar < 0.5)")
ok = e2 >= 0.98*2*D8 and es0 >= 0.98*2*D8 and spc < 0.5
print(f"L4(b) {'PASS' if ok else 'MISS'}")
# report: first space-carrying level, n = 6 without truncation
st, ix = basis_ring(6, sz=0); H = hamiltonian(6, g, st, ix)
w, V = np.linalg.eigh(H.toarray()); C, _, _ = diag_obs(6, st)
sc = (V**2*C.sum(1)[:, None]).sum(0)
hit = np.where(sc > 0.5)[0][0]
st2, ix2 = basis_ring(6, sz=2); w2 = np.linalg.eigvalsh(hamiltonian(6, g, st2, ix2).toarray())
D6 = w2[0] - w[0]
print(f"  report: n = 6 (exact, S^z = 0): first level with <sum c> > 0.5 is index {hit} of {len(w)}, "
      f"E - E0 = {w[hit]-w[0]:.5f} = {(w[hit]-w[0])/DELTA:.3f} Delta = {(w[hit]-w[0])/(2*D6):.1f} x 2 Delta_6 "
      f"(spin-manifold S^z=0 states: 141)  [{time.time()-t0:.0f}s]")
