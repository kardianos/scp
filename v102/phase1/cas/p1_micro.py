#!/usr/bin/env python3
"""P1 and P4 — exact diagonalization.  See ../PREREG.md.

P1   microscopic ring with transport.  Cells: two bosonic modes, cap 2
     (R6: two link ends).  Links: one scalar space mode, c_l <= 2, and
     sum c <= cmax.
     H = Delta sum c^+c + g sum (c^+ A + h.c.) + t sum_a (b_ia^+ b_ja + h.c.)
     Transport above the cap is refused.  Sector: full occupation
     Q = sum N + 2 sum c = 2n, and S^z = 0.  The first 8 distinct
     excitations are mapped onto the spin-1 Heisenberg ring with one
     fitted J_eff.
P4a  lowest level carrying a real space quantum (<space> > 0.5),
     against twice the magnon gap.
P4b  effective spin-1 rings, n = 8, 10, 12: lowest S = 0 excitation and
     lowest S = 2 state, against 2*Delta_n.
"""
import itertools, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

t0 = time.time()
def stamp(): return f"[{time.time()-t0:5.0f}s]"


# ---------------------------------------------------------------- microscopic ring
def micro(n, g, t, Delta=1.0, cmax=2):
    by_N = {N: [(a, N - a) for a in range(N + 1)] for N in range(3)}
    basis, idx = [], {}
    for cs in itertools.product(range(3), repeat=n):
        sc = sum(cs)
        if sc > cmax: continue
        rem = 2*n - 2*sc
        for Ns in itertools.product(range(3), repeat=n):
            if sum(Ns) != rem: continue
            for cells in itertools.product(*[by_N[N] for N in Ns]):
                if sum(a - b for a, b in cells) != 0: continue
                idx[cs + sum(cells, ())] = len(basis); basis.append((cs, cells))
    D = len(basis)
    rows, cols, vals = [], [], []
    space = np.zeros(D)
    for k, (cs, cells) in enumerate(basis):
        space[k] = sum(cs)
        rows.append(k); cols.append(k); vals.append(Delta*sum(cs))
        for l in range(n):
            i, j = l, (l + 1) % n
            if sum(cs) < cmax and cs[l] < 2:                       # c_l^+ A_l (+ h.c.)
                for (ai, aj, sgn) in ((0, 1, +1), (1, 0, -1)):
                    if cells[i][ai] == 0 or cells[j][aj] == 0: continue
                    amp = sgn*g*np.sqrt(cells[i][ai]*cells[j][aj]*(cs[l] + 1))
                    nc = list(cells)
                    x = list(nc[i]); x[ai] -= 1; nc[i] = tuple(x)
                    y = list(nc[j]); y[aj] -= 1; nc[j] = tuple(y)
                    cs2 = list(cs); cs2[l] += 1
                    k2 = idx.get(tuple(cs2) + sum(nc, ()))
                    if k2 is None: continue
                    rows += [k2, k]; cols += [k, k2]; vals += [amp, amp]
            if t != 0:                                             # transport, cap-refused
                for (src, dst) in ((i, j), (j, i)):
                    if sum(cells[dst]) >= 2: continue
                    for a in (0, 1):
                        if cells[src][a] == 0: continue
                        amp = t*np.sqrt(cells[src][a]*(cells[dst][a] + 1))
                        nc = list(cells)
                        x = list(nc[src]); x[a] -= 1; nc[src] = tuple(x)
                        y = list(nc[dst]); y[a] += 1; nc[dst] = tuple(y)
                        k2 = idx.get(cs + sum(nc, ()))
                        if k2 is None: continue
                        rows.append(k2); cols.append(k); vals.append(amp)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(D, D))
    return H, space, basis


def S2_sz0(v, basis):
    """<S^2> in the S^z = 0 sector: S^2 = S^- S^+, so <S^2> = |S^+ v|^2."""
    acc = {}
    for k, (cs, cells) in enumerate(basis):
        if abs(v[k]) < 1e-13: continue
        for i, (a, b) in enumerate(cells):
            if b == 0: continue
            nc = list(cells); nc[i] = (a + 1, b - 1)
            key = cs + sum(nc, ())
            acc[key] = acc.get(key, 0.0) + v[k]*np.sqrt(b*(a + 1))
    s2 = sum(x*x for x in acc.values())
    return -0.5 + np.sqrt(0.25 + s2)


def lowest(H, k):
    D = H.shape[0]
    if D <= 3000:
        w, v = np.linalg.eigh(H.toarray()); return w[:k], v[:, :k]
    w, v = eigsh(H, k=k, which='SA', tol=1e-12)
    o = np.argsort(w); return w[o], v[:, o]


def distinct(w, tol):
    out = [w[0]]
    for x in w[1:]:
        if x - out[-1] > tol: out.append(x)
    return np.array(out)


# ---------------------------------------------------------------- effective rings
def spin1(n):
    Sz = sp.diags([1.0, 0.0, -1.0]); Sp = sp.csr_matrix(np.diag([np.sqrt(2), np.sqrt(2)], 1))
    at = lambda op, i: sp.kron(sp.kron(sp.identity(3**i), op), sp.identity(3**(n - i - 1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((3**n, 3**n))
    for i in range(n):
        j = (i + 1) % n
        H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
    return H, sum(Z).diagonal()

def sector(H, szd, m):
    keep = np.where(np.abs(szd - m) < 1e-9)[0]
    return H[keep][:, keep]


# ---------------------------------------------------------------- P1
if __name__ == "__main__":
    print(f"# P1 microscopic ring vs spin-1 Heisenberg ring (S^z = 0 sectors) {stamp()}", flush=True)
    heis = {}
    for n in (4, 6):
        H, szd = spin1(n)
        w = np.linalg.eigvalsh(sector(H, szd, 0).toarray())
        heis[n] = distinct(w, 1e-9)
    g = 0.1
    P1_ok = True
    P4a = {}
    for (n, cmax, ts) in ((4, 4, (0.0, 0.1, 0.2)), (6, 2, (0.0, 0.1, 0.2))):
        hx = heis[n][1:9] - heis[n][0]
        for t in ts:
            H, space, basis = micro(n, g, t, cmax=cmax)
            kk = 40 if n == 4 else 30
            w, v = lowest(H, kk)
            lv = distinct(w, 1e-9*abs(g*g))
            ex = lv[1:9] - lv[0]
            J = float(ex@hx/(hx@hx))
            miss = np.abs(ex/(J*hx) - 1)
            S0 = S2_sz0(v[:, 0], basis)
            sp0 = float(v[:, 0]**2@space)
            ok = (miss.max() < 0.05) and abs(S0) < 1e-6
            P1_ok &= ok
            print(f"  n={n} cmax={cmax} t/D={t:.1f}  dim={H.shape[0]:7d}  J_eff={J:.6f} (2g^2/D={2*g*g:.4f})"
                  f"  max miss over 8 levels={miss.max():.2e}  S_ground={S0:.4f}  <space>_ground={sp0:.4f}"
                  f"  {'PASS' if ok else 'MISS'}  {stamp()}", flush=True)
            if n == 6 and t == 0.1:
                P4a['H'] = H; P4a['space'] = space; P4a['J'] = J; P4a['gap'] = ex[0]
    print(f"P1 {'PASS' if P1_ok else 'MISS'}")

    # ---------------------------------------------------------------- P4a
    print(f"\n# P4a lowest level with a real space quantum, n=6 g/D=0.1 t/D=0.1 {stamp()}", flush=True)
    H = P4a['H']; space = P4a['space']
    w, v = eigsh(H, k=175, which='SA', tol=1e-10)
    o = np.argsort(w); w = w[o]; v = v[:, o]
    spv = (v**2*space[:, None]).sum(0)
    hit = np.where(spv > 0.5)[0]
    mag = P4a['gap']
    if len(hit):
        e = w[hit[0]] - w[0]
        print(f"  first level with <space> > 0.5: index {hit[0]}, E - E0 = {e:.5f}  (Delta = 1)")
        print(f"  magnon gap on this ring = {mag:.6f};  ratio to 2*gap = {e/(2*mag):.1f}  "
              f"{'PASS' if e > 2*mag else 'MISS'}")
    else:
        print(f"  no level with <space> > 0.5 among the lowest 175 (top at E-E0 = {w[-1]-w[0]:.4f}); "
              f"all lie above 2*gap = {2*mag:.5f}: PASS")

    # ---------------------------------------------------------------- P4b
    print(f"\n# P4b effective spin-1 rings: second line? {stamp()}", flush=True)
    print("   n   Delta_n     lowest S=0 exc   /(2 Delta_n)   lowest S=2   /(2 Delta_n)")
    P4b_ok = True
    for n in (8, 10, 12):
        H, szd = spin1(n)
        w0 = np.sort(eigsh(sector(H, szd, 0), k=14, which='SA', tol=1e-12)[0])
        w1 = np.sort(eigsh(sector(H, szd, 1), k=14, which='SA', tol=1e-12)[0])
        w2 = np.sort(eigsh(sector(H, szd, 2), k=2, which='SA', tol=1e-12)[0])
        E0 = w0[0]; Dn = w1[0] - E0
        s0 = [x for x in w0[1:] if np.min(np.abs(w1 - x)) > 1e-8]
        eS0 = s0[0] - E0 if s0 else float('nan')
        eS2 = w2[0] - E0
        ok = (eS0 >= 0.98*2*Dn) and (eS2 >= 0.98*2*Dn)
        P4b_ok &= ok
        print(f"  {n:3d}  {Dn:.6f}    {eS0:.6f}       {eS0/(2*Dn):.4f}        {eS2:.6f}   {eS2/(2*Dn):.4f}   "
              f"{'ok' if ok else 'BELOW'}  {stamp()}", flush=True)
    print(f"P4b(rings) {'PASS' if P4b_ok else 'MISS'}")
