#!/usr/bin/env python3
"""v99 reopen — numbers for options A, B, A+B (see ../OPTIONS.md).

E1  A alone: two-mode (Jones-vector) cells, quantized as bosons, with a
    spin-blind hop of either sign and a spin-blind on-site load term.
    Ground-state total spin.  Eisenberg-Lieb predicts maximal spin
    (ferromagnet) when the hop is sign-free (-t, or any sign on an even
    ring).  The fabric's hop is +t (H = +J(b_i^+ b_j + h.c.), from
    z_i' = cos z_i - i sin z_j), which is frustrated on odd rings.
E2  B / A+B: spin-1 and spin-1/2 antiferromagnetic Heisenberg rings.
    Gap vs n (uniform for S=1, closing for S=1/2) and ground energy vs n
    (extensive: no interior n*).
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

def lowest(H, k):
    k = min(k, H.shape[0] - 2)
    if H.shape[0] <= 400:
        w, v = np.linalg.eigh(H.toarray()); return w[:k], v[:, :k]
    w, v = eigsh(H, k=k, which='SA', tol=1e-10)
    o = np.argsort(w); return w[o], v[:, o]

# ---------------------------------------------------------------- E1
print("# E1  A alone: two-mode boson rings, spin-blind hop t and load U")
print("#     S_tot of the ground manifold; max = N/2 (ferromagnet)")
def boson_ring(n, N, t, U):
    M = 2*n                                  # mode (site i, comp a) -> 2i+a
    basis = []
    for c in itertools.combinations_with_replacement(range(M), N):
        occ = [0]*M
        for m in c: occ[m] += 1
        basis.append(tuple(occ))
    idx = {b: k for k, b in enumerate(basis)}
    D = len(basis)
    rows, cols, vals = [], [], []
    srows, scols, svals = [], [], []
    for k, b in enumerate(basis):
        # on-site spin-blind load
        e = 0.0
        for i in range(n):
            ni = b[2*i] + b[2*i+1]; e += 0.5*U*ni*(ni-1)
        rows.append(k); cols.append(k); vals.append(e)
        # hop b_{j,a}^+ b_{i,a}, both directions, ring
        for i in range(n):
            for j in ((i+1) % n, (i-1) % n):
                if n == 2 and j == (i-1) % n: continue
                for a in (0, 1):
                    src, dst = 2*i+a, 2*j+a
                    if b[src] == 0: continue
                    nb = list(b); amp = np.sqrt(nb[src]); nb[src] -= 1
                    amp *= np.sqrt(nb[dst]+1); nb[dst] += 1
                    rows.append(idx[tuple(nb)]); cols.append(k); vals.append(t*amp)
        # total S+ = sum_i b_{i,0}^+ b_{i,1}
        for i in range(n):
            src, dst = 2*i+1, 2*i
            if b[src] == 0: continue
            nb = list(b); amp = np.sqrt(nb[src]); nb[src] -= 1
            amp *= np.sqrt(nb[dst]+1); nb[dst] += 1
            srows.append(idx[tuple(nb)]); scols.append(k); svals.append(amp)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(D, D))
    Sp = sp.csr_matrix((svals, (srows, scols)), shape=(D, D))
    Sz = sp.diags([0.5*sum(b[2*i]-b[2*i+1] for i in range(n)) for b in basis])
    S2 = Sp @ Sp.T + Sz @ Sz - Sz
    return H, S2, D

print(" n  N  S/cell   t     U    dim   E0          S_tot over ground level   max")
for n, N in [(3, 3), (4, 4), (5, 5), (3, 6), (4, 8)]:
    for t in (-1.0, +1.0):
        for U in (0.0, 4.0, 20.0):
            H, S2, D = boson_ring(n, N, t, U)
            w, v = lowest(H, min(N + 4, D - 2))
            g = np.sum(w < w[0] + 1e-8)
            # S of each returned ground vector, after diagonalizing S^2 in the level
            blk = v[:, :g].T @ (S2 @ v[:, :g]); s2 = np.linalg.eigvalsh((blk + blk.T)/2)
            Ss = -0.5 + np.sqrt(0.25 + np.clip(s2, 0, None))
            rng_ = f"{Ss.min():.3f}" if Ss.max() - Ss.min() < 1e-6 else f"{Ss.min():.3f}..{Ss.max():.3f}"
            print(f" {n}  {N:2d}  {N/n/2:4.1f}  {t:+4.0f}  {U:5.1f}  {D:5d}  {w[0]:9.4f}   {rng_:>14s}           {N/2:4.1f}")

# ---------------------------------------------------------------- E2
print()
print("# E2  AF Heisenberg rings: gap and ground energy vs n")
def spin_ops(s):
    m = np.arange(s, -s-1, -1); d = len(m)
    Sz = sp.diags(m)
    Sp = sp.lil_matrix((d, d))
    for k in range(1, d):
        Sp[k-1, k] = np.sqrt(s*(s+1) - m[k]*(m[k]+1))
    return sp.csr_matrix(Sz), sp.csr_matrix(Sp)

def heis_ring(n, s):
    Sz, Sp = spin_ops(s); d = Sz.shape[0]
    def at(op, i):
        return sp.kron(sp.kron(sp.identity(d**i), op), sp.identity(d**(n-i-1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((d**n, d**n))
    for i in range(n):
        j = (i+1) % n
        H = H + Z[i] @ Z[j] + 0.5*(P[i] @ P[j].T + P[i].T @ P[j])
    Szt = sum(Z)
    return H, Szt

def sector(H, Szt, m):
    keep = np.where(np.abs(Szt.diagonal() - m) < 1e-9)[0]
    return H[keep][:, keep]

for s, ns in [(1.0, range(3, 13)), (0.5, range(4, 19, 2))]:
    print(f" S={s}:   n    E0(n)       E0/n      gap(n)     n*gap    [t]")
    E = {}
    for n in ns:
        t0 = time.time()
        H, Szt = heis_ring(n, s)
        m0 = 0.0 if (s*n) % 1 == 0 else 0.5
        e0 = lowest(sector(H, Szt, m0), 1)[0][0]
        e1 = lowest(sector(H, Szt, m0 + 1), 1)[0][0]
        E[n] = e0
        print(f"        {n:2d}  {e0:10.5f}  {e0/n:9.5f}  {e1-e0:9.5f}  {n*(e1-e0):8.3f}   [{time.time()-t0:4.1f}s]")
    if s == 1.0:
        ev = [n for n in E if n % 2 == 0]
        print("   second difference of E0 on even n (E(n+2)-2E(n)+E(n-2)):")
        for n in ev[1:-1]:
            print(f"        n={n:2d}: {E[n+2]-2*E[n]+E[n-2]:+.5f}")
        print("   E0(n) - n*e_inf with e_inf = -1.401484 (literature):")
        for n in E:
            print(f"        n={n:2d}: {E[n] - n*(-1.401484):+.5f}")

# ---------------------------------------------------------------- E3
print()
print("# E3  B alone: spin-1/2 per link end, two per cell, distinguishable.")
print("#     Alternating ring: J_link = +1 (AF) across links, J_cell inside a cell.")
print("#     J_cell -> -inf is the spin-1 chain (gap -> 0.4105/4 = 0.103);")
print("#     J_cell = 0 is isolated link singlets (gap 1); J_cell = +1 is the")
print("#     uniform spin-1/2 chain (gapless).  The gap moves continuously.")
def alt_ring(ncell, Jc):
    n = 2*ncell
    Sz, Sp = spin_ops(0.5)
    def at(op, i):
        return sp.kron(sp.kron(sp.identity(2**i), op), sp.identity(2**(n-i-1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((2**n, 2**n))
    for i in range(n):
        j = (i+1) % n
        J = Jc if i % 2 == 0 else 1.0       # (2c, 2c+1) inside cell c; (2c+1, 2c+2) across a link
        H = H + J*(Z[i] @ Z[j] + 0.5*(P[i] @ P[j].T + P[i].T @ P[j]))
    return H, sum(Z)
print("  J_cell   gap(4 cells)  gap(6 cells)  gap(8 cells)")
for Jc in (1.0, 0.5, 0.0, -0.5, -1.0, -2.0, -4.0, -8.0):
    gs = []
    for nc in (4, 6, 8):
        H, Szt = alt_ring(nc, Jc)
        e0 = lowest(sector(H, Szt, 0.0), 1)[0][0]
        e1 = lowest(sector(H, Szt, 1.0), 1)[0][0]
        gs.append(e1 - e0)
    print(f"  {Jc:+5.1f}   {gs[0]:10.4f}   {gs[1]:10.4f}   {gs[2]:10.4f}")
print("E3 verdict: with the on-site rule left free, B's gap (and so xi) is a")
print("   continuous function of J_cell/J_link. A's bosonic modes fix J_cell = -inf.")
