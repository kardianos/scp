#!/usr/bin/env python3
"""L1 — the S-law ring's spectrum against the spin-1 Heisenberg ring
(../PREREG.md §3 L1).  g/Delta = 0.1, t = 0, S^z = 0, first 12
excitations with multiplicity, one fitted J."""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import DELTA, basis_ring, hamiltonian, total_S, heis_ring, fit_J
from clusters import effective

t0 = time.time()
gr = 0.1; g = gr*DELTA
h2, r3, co, info = effective(g)
print(f"# L1  S-law ring vs spin-1 Heisenberg ring; g/Delta = {gr}, g = {g:.5f}, Delta = {DELTA:.5f}")
print(f"  two-cell cluster: J_h2 = {co['J']:.6f}  (2g^2/Delta = {2*g*g/DELTA:.6f});  K/J = {co['K']/co['J']:+.5f};  "
      f"three-site |r3|/J = {np.abs(np.linalg.eigvalsh(r3[(2,2,2)])).max()/co['J']:.4f}")

def heis_levels(n, k=13):
    H, Z, P = heis_ring(n)
    szd = sum(Z).diagonal(); keep = np.where(np.abs(szd) < 1e-9)[0]
    Hs = H[keep][:, keep]
    if Hs.shape[0] <= 4000: w = np.linalg.eigvalsh(Hs.toarray())[:k]
    else: w = np.sort(eigsh(Hs, k=k, which='SA', tol=1e-12)[0])
    return w

ok_all = True
rows = []
for (n, cmax) in ((6, None), (8, 3), (8, 2)):
    st, ix = basis_ring(n, cmax=cmax, sz=0)
    H = hamiltonian(n, g, st, ix)
    if H.shape[0] <= 4000:
        w, v = np.linalg.eigh(H.toarray()); w = w[:20]; v = v[:, :20]
    else:
        w, v = eigsh(H, k=20, which='SA', tol=1e-12); o = np.argsort(w); w = w[o]; v = v[:, o]
    hw = heis_levels(n)
    hx = hw[1:13] - hw[0]; ex = w[1:13] - w[0]
    J, miss = fit_J(ex, hx)
    S0 = total_S(v[:, 0], st, ix)
    ok = miss.max() < 0.05 and abs(S0) < 1e-6
    gate = cmax != 2
    if gate: ok_all &= ok
    print(f"  n={n} cmax={'none' if cmax is None else cmax:>4} dim={H.shape[0]:6d}  E0={w[0]:+.7f}  J_fit={J:.6f}  "
          f"J_fit/J_h2={J/co['J']:.4f}  max miss (12 levels)={miss.max():.3e}  S_ground={S0:.4f}  "
          f"{('PASS' if ok else 'MISS') if gate else '(convergence check, not gated)'}  [{time.time()-t0:.0f}s]", flush=True)
    rows.append((n, cmax, w[:13], J))
print("  levels (E - E0)/J_fit vs Heisenberg ring:")
for (n, cmax, w, J) in rows:
    hw = heis_levels(n)
    print(f"   n={n} cmax={cmax}: " + " ".join(f"{(a - w[0])/J:.4f}/{b - hw[0]:.4f}" for a, b in zip(w[1:13], hw[1:13])))
print(f"L1 {'PASS' if ok_all else 'MISS'}")
