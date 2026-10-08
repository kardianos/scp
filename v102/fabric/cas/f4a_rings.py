#!/usr/bin/env python3
"""F4(a) — the hole on effective rings (reading V), n = 6..14, g/Delta = 0.1.
Ring of n cells, cell 0 at cap - 1 (spin 1/2), the rest full (spin 1);
h2 + r3 including the hole-centred (2,1,2).  Levels in 2S^z = 1, 3, 5;
a level's spin is the largest S^z sector it appears in.  The no-hole
ring's gap (same n) gives the approximate one-magnon threshold.
Also the spin profile of the S^z = 1/2 ground state."""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
import eff

t0 = time.time()
gr = 0.1
h2, r3, co = eff.terms_for(gr); J = co['J']
print(f"# F4(a) one hole on effective rings, g/Delta = {gr}, J = {J:.6f}")

def low(H, k):
    if H.shape[0] <= 2000:
        w, v = np.linalg.eigh(H.toarray()); return w[:k], v[:, :k]
    w, v = eigsh(H, k=k, which='SA', tol=1e-11); o = np.argsort(w); return w[o], v[:, o]

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
for n in range(6, nmax + 1, 2):
    caps0 = [2]*n
    H0, *_ = eff.sector_hamiltonian(caps0, eff.chain_terms(caps0, h2, r3, ring=True), 0)
    H1, *_ = eff.sector_hamiltonian(caps0, eff.chain_terms(caps0, h2, r3, ring=True), 2)
    e00 = low(H0, 1)[0][0]; e01 = low(H1, 1)[0][0]
    gap0 = e01 - e00
    caps = [1] + [2]*(n - 1)
    T = eff.chain_terms(caps, h2, r3, ring=True)
    lev = {}
    for tw in (1, 3, 5):
        H, basis, dims, strides = eff.sector_hamiltonian(caps, T, tw)
        k = 10 if tw == 1 else 6
        w, v = low(H, k)
        lev[tw] = w
        if tw == 1:
            prof = eff.sz_profile(v[:, 0], caps, basis, dims, strides)
            dim1 = H.shape[0]
    E0 = lev[1][0]
    def spin_of(e, tw):
        S = tw/2
        for up in (tw + 2, tw + 4):
            if up in lev and np.any(np.abs(lev[up] - e) < 1e-8*max(1, abs(e)) + 1e-10): S = up/2
        return S
    rows = []
    for tw in (1, 3, 5):
        for e in lev[tw]:
            S = spin_of(e, tw)
            if abs(S - tw/2) < 1e-9: rows.append((e, S))
    rows.sort()
    thr = E0 + gap0
    print(f"  n={n:2d} (2S^z=1 dim {dim1}): ring without hole: gap {gap0/J:.4f} J.  With the hole: ground S = {rows[0][1]:.1f}; "
          f"levels (E-E0 in J, S): " + "  ".join(f"{(e-E0)/J:.4f}(S={S:.1f})" for e, S in rows[:8]) +
          f";  approx one-magnon threshold {gap0/J:.4f} J  [{time.time()-t0:.0f}s]", flush=True)
    r = np.array([min(i, n - i) for i in range(n)])
    half = n//2
    pr = [prof[(r == d)].sum() for d in range(half + 1)]
    print(f"        S^z profile by distance from the hole (sum over both sides): " + " ".join(f"{x:+.4f}" for x in pr) +
          f";  total {prof.sum():+.4f}", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
