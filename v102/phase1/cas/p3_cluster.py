#!/usr/bin/env python3
"""P3/P5 by a better-conditioned route: the exact low-energy effective
Hamiltonian of the microscopic model, then DMRG on the effective chain.

1. Open clusters of k = 2 and k = 3 cells (caps from R6), each link with
   a scalar space mode c <= 2, full occupation.  Exact diagonalization;
   the lowest prod(cap+1) states span the spin manifold (separated by
   about Delta from states carrying space).  The des Cloizeaux effective
   Hamiltonian on the model space (all cells full, no space) is
   H_eff = (X X^T)^-1/2 X E X^T (X X^T)^-1/2.
2. h2 is the k = 2 result (one link): it decomposes exactly into
   a + J S.S + K (S.S)^2.  r3 = h3 - h2(12) - h2(23) is the genuinely
   two-link part.  A linked-cluster sum through k = 3 is exact to fourth
   order (g^4, g^2 t^2).
3. Effective chain: H = sum h2 + sum r3 (degree-1 end cells use their
   own cap-1 clusters), divided by J.  DMRG (validated in P2) gives xi
   with the same window as the effective Heisenberg chain.
Usage: p3_cluster.py [L] [chi]
"""
import sys, itertools, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from dmrg import DMRG, build_mpo, spin_ops

t0 = time.time()
L = int(sys.argv[1]) if len(sys.argv) > 1 else 48
CHI = int(sys.argv[2]) if len(sys.argv) > 2 else 48


def cluster(caps, g, t, Delta=1.0, cmax_link=2):
    k = len(caps); nl = k - 1; Q = sum(caps)
    cells_by = [[(a, N - a) for N in range(cap + 1) for a in range(N + 1)] for cap in caps]
    basis, idx = [], {}
    for cs in itertools.product(range(cmax_link + 1), repeat=nl):
        for cells in itertools.product(*cells_by):
            if sum(sum(c) for c in cells) + 2*sum(cs) != Q: continue
            idx[(cs, cells)] = len(basis); basis.append((cs, cells))
    D = len(basis); H = np.zeros((D, D))
    for q, (cs, cells) in enumerate(basis):
        H[q, q] += Delta*sum(cs)
        for l in range(nl):
            i, j = l, l + 1
            if cs[l] < cmax_link:
                for (ai, aj, sgn) in ((0, 1, +1), (1, 0, -1)):
                    if cells[i][ai] == 0 or cells[j][aj] == 0: continue
                    amp = sgn*g*np.sqrt(cells[i][ai]*cells[j][aj]*(cs[l] + 1))
                    nc = list(cells)
                    x = list(nc[i]); x[ai] -= 1; nc[i] = tuple(x)
                    y = list(nc[j]); y[aj] -= 1; nc[j] = tuple(y)
                    cs2 = list(cs); cs2[l] += 1
                    q2 = idx[(tuple(cs2), tuple(nc))]
                    H[q2, q] += amp; H[q, q2] += amp
            if t != 0:
                for (src, dst) in ((i, j), (j, i)):
                    if sum(cells[dst]) >= caps[dst]: continue
                    for a in (0, 1):
                        if cells[src][a] == 0: continue
                        amp = t*np.sqrt(cells[src][a]*(cells[dst][a] + 1))
                        nc = list(cells)
                        x = list(nc[src]); x[a] -= 1; nc[src] = tuple(x)
                        y = list(nc[dst]); y[a] += 1; nc[dst] = tuple(y)
                        H[idx[(cs, tuple(nc))], q] += amp
    w, V = np.linalg.eigh(H)
    # model space, ordered as the kron of spin states m = S, S-1, ..., -S  (n1 = cap, cap-1, ..., 0)
    model = []
    for n1s in itertools.product(*[range(cap, -1, -1) for cap in caps]):
        cells = tuple((n1, cap - n1) for n1, cap in zip(n1s, caps))
        model.append(idx[((0,)*nl, cells)])
    m = len(model)
    gapsep = w[m] - w[m - 1]
    X = V[model, :m]
    S = X@X.T
    ev, U = np.linalg.eigh(S)
    Sm = U@np.diag(ev**-0.5)@U.T
    Heff = Sm@X@np.diag(w[:m])@X.T@Sm
    return 0.5*(Heff + Heff.T), gapsep, float(ev.min())


def spinvec(s):
    Sz, Sp = spin_ops(s); return Sz, Sp

def SdotS(s1, s2):
    Z1, P1 = spinvec(s1); Z2, P2 = spinvec(s2)
    return np.kron(Z1, Z2) + 0.5*(np.kron(P1, P2.T) + np.kron(P1.T, P2))

def op_svd2(h, d1, d2, tol=1e-12):
    T = h.reshape(d1, d2, d1, d2).transpose(0, 2, 1, 3).reshape(d1*d1, d2*d2)
    U, s, Vt = np.linalg.svd(T)
    out = []
    for k in range(len(s)):
        if s[k] > tol*max(s[0], 1e-300):
            out.append((s[k], U[:, k].reshape(d1, d1), Vt[k].reshape(d2, d2)))
    return out

def op_svd3(h, d1, d2, d3, tol=1e-12):
    T = h.reshape(d1, d2, d3, d1, d2, d3).transpose(0, 3, 1, 4, 2, 5).reshape(d1*d1, d2*d2*d3*d3)
    U, s, Vt = np.linalg.svd(T, full_matrices=False)
    out = []
    for k in range(len(s)):
        if s[k] <= tol*s[0]: continue
        M = Vt[k].reshape(d2*d2, d3*d3)
        U2, s2, V2 = np.linalg.svd(M, full_matrices=False)
        for l in range(len(s2)):
            if s2[l] <= tol*s2[0]: continue
            out.append((s[k]*s2[l], U[:, k].reshape(d1, d1), U2[:, l].reshape(d2, d2), V2[l].reshape(d3, d3)))
    return out


def effective(g, t):
    caps_bulk = 2
    h2 = {}; info = {}
    for pair in ((2, 2), (1, 2), (2, 1)):
        h, sep, evmin = cluster(list(pair), g, t)
        h2[pair] = h; info[pair] = (sep, evmin)
    # decompose the bulk h2
    SS = SdotS(1, 1); B = [np.eye(9), SS, SS@SS]
    A = np.stack([b.ravel() for b in B], 1)
    coef, *_ = np.linalg.lstsq(A, h2[(2, 2)].ravel(), rcond=None)
    resid = np.abs(h2[(2, 2)] - sum(c*b for c, b in zip(coef, B))).max()
    a, J, K = coef
    r3 = {}
    for trip in ((2, 2, 2), (1, 2, 2), (2, 2, 1)):
        h3, sep, evmin = cluster(list(trip), g, t)
        d = [c + 1 for c in trip]
        h12 = np.kron(h2[(trip[0], trip[1])], np.eye(d[2]))
        h23 = np.kron(np.eye(d[0]), h2[(trip[1], trip[2])])
        r3[trip] = h3 - h12 - h23
        info[trip] = (sep, evmin)
    return h2, r3, (a, J, K, resid), info


def chain_mpo(h2, r3, J, L):
    caps = [1] + [2]*(L - 2) + [1]
    dims = [c + 1 for c in caps]
    terms = []
    for i in range(L - 1):
        h = h2[(caps[i], caps[i + 1])]/J
        for (s, A, B) in op_svd2(h, dims[i], dims[i + 1]):
            terms.append((i, [A, B], s))
    for i in range(L - 2):
        trip = (caps[i], caps[i + 1], caps[i + 2])
        r = r3[trip]/J
        if np.abs(r).max() < 1e-14: continue
        for (s, A, B, C) in op_svd3(r, dims[i], dims[i + 1], dims[i + 2], tol=1e-10):
            terms.append((i, [A, B, C], s))
    W, c = build_mpo(dims, [None]*L, terms)
    return W, dims, caps


def xi_of(d, caps, rs=np.arange(2, 25), win=(6, 20)):
    mid = len(caps)//2
    ops = [spin_ops(c/2) for c in caps]
    Z = [o[0] for o in ops]; P = [o[1] for o in ops]
    C = []
    for r in rs:
        i = mid - r//2; j = i + r
        fz = [None]*len(caps); fm = [None]*len(caps); fz[j] = Z[j]; fm[j] = P[j].T
        C.append(d.corr_row(Z[i], i, fz, [j])[j] + d.corr_row(P[i], i, fm, [j])[j])
    C = np.array(C); m = (rs >= win[0]) & (rs <= win[1])
    return -1/np.polyfit(rs[m], np.log(np.abs(C[m])*np.sqrt(rs[m])), 1)[0], C


if __name__ == "__main__":
    print(f"# P3/P5 via the exact cluster effective Hamiltonian; chain L={L}, chi={CHI}")
    print("  g/D   t/D    J          2g^2/D    K/J (biquad)   SU2 resid   |r3|/J (bulk)   gap-to-space   xi(chain)   E/J")
    pts = [(0.05, 0.0), (0.1, 0.0), (0.2, 0.0), (0.1, 0.1), (0.2, 0.2), (0.3, 0.0)]
    if len(sys.argv) > 3:
        pts = [tuple(map(float, p.split(','))) for p in sys.argv[3:]]
    for (g, t) in pts:
        h2, r3, (a, J, K, resid), info = effective(g, t)
        r3n = np.abs(np.linalg.eigvalsh(r3[(2, 2, 2)])).max()/J
        W, dims, caps = chain_mpo(h2, r3, J, L)
        d = DMRG(W, dims, chi=CHI, seed=11)
        E = d.run(sweeps=14, etol=1e-11)
        xi, C = xi_of(d, caps)
        print(f"  {g:4.2f}  {t:4.2f}  {J:.6f}  {2*g*g:.6f}   {K/J:+.5f}      {resid:.1e}     {r3n:.4f}          "
              f"{info[(2,2,2)][0]:.4f}       {xi:7.4f}    {E:.6f}   [{time.time()-t0:.0f}s]", flush=True)
    # reference: pure Heisenberg chain, same L and window
    h2H = {(2, 2): SdotS(1, 1), (1, 2): SdotS(0.5, 1), (2, 1): SdotS(1, 0.5)}
    r3H = {k: np.zeros((np.prod([c + 1 for c in k]),)*2) for k in ((2, 2, 2), (1, 2, 2), (2, 2, 1))}
    W, dims, caps = chain_mpo(h2H, r3H, 1.0, L)
    d = DMRG(W, dims, chi=CHI, seed=11); E = d.run(sweeps=14, etol=1e-11)
    xi, C = xi_of(d, caps)
    print(f"  Heisenberg reference (K=0, r3=0): xi = {xi:.4f}   E = {E:.6f}")
