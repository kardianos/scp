#!/usr/bin/env python3
"""Leading-order (second-order in g) effective model of the S-law fiber for
M = 2 (3+2) or M = 4 (the widening W), with an optional mobile hole (T).

Sites: 'full' (N = 2), 'hole' (N = 1), 'end' (R6 end cell, cap 1, N = 1,
never hops) or 'mob' (N = 2 or N = 1: the mobile-hole site space).
Link exchange: -(g^2/Delta_type) A^+A,  A = sum P[mu,nu] b_i mu b_j nu,
written as sum_{(rho mu),(sigma nu)} P[rho,sigma] P[mu,nu] E^i_{rho mu} E^j_{sigma nu}
with E_{rho mu} = b^+_rho b_mu (it keeps each cell's N), and Delta_type the
first-conversion cost of the pair: (2,2) -> Delta, (2,1)/(1,2) -> 0.2275,
(1,1) -> -0.0243 (not perturbative; never used with one hole).
The (1,1) link (two adjacent N = 1 cells) is priced by P11 (../../bind/):
'zero' (default, the original: no term), 'core' (penalty P11_CORE on every
(1,1) configuration), 'score' (penalty P11_CORE on the singlet pair only, via
(P11_CORE/|P|^2) A^+A), 'upper' / 'lower' (the link's exact two-level hybrid of
the singlet pair and one space quantum, field-like / space-like branch:
E = D11/2 +- sqrt(D11^2/4 + g^2 |P|^2), applied as (E/|P|^2) A^+A).
Hop (T, T2 rule, mob sites only): t sum_mu (b^+_j mu b_i mu + h.c.), N = 2 -> 1
and 1 -> 2 (the hole moves); b on N = 1 is truncated (one hole).

Builds: (i) a term list for ../../phase1/cas/dmrg.py build_mpo (with any
number of penalties), (ii) exact sparse rings/chains in a sector of the
conserved labels (2S^z[, 2I^z], N_holes)."""
import sys, itertools
import numpy as np
import scipy.sparse as sp
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../phase1/cas')
import graphm as GM

DELTA_T = {(2, 2): GM.OMEGA - 2*(GM.EF[2] - GM.EF[1]),
           (2, 1): GM.OMEGA - (GM.EF[2] - GM.EF[1]) - GM.EF[1],
           (1, 2): GM.OMEGA - (GM.EF[2] - GM.EF[1]) - GM.EF[1],
           (1, 1): GM.OMEGA - 2*GM.EF[1]}

class Local:
    """A site space: a list of occupation tuples (N = 2 states first)."""
    def __init__(self, M, kind):
        self.M = M; self.kind = kind
        if kind == 'full': self.occ = GM.cell_states(M, 2)
        elif kind in ('hole', 'end'): self.occ = GM.cell_states(M, 1)
        elif kind == 'mob': self.occ = GM.cell_states(M, 2) + GM.cell_states(M, 1)
        self.d = len(self.occ); self.idx = {o: k for k, o in enumerate(self.occ)}
        self.N = np.array([sum(o) for o in self.occ])
        self.E = {}
        for r in range(M):
            for m in range(M):
                X = np.zeros((self.d, self.d))
                for k, o in enumerate(self.occ):
                    if o[m] == 0: continue
                    o2 = list(o); o2[m] -= 1; o2[r] += 1; o2 = tuple(o2)
                    if o2 in self.idx:
                        X[self.idx[o2], k] += np.sqrt(o[m]*(o2[r] if r != m else o[m]))
                self.E[(r, m)] = X
        self.b = []; self.bd = []
        for m in range(M):
            X = np.zeros((self.d, self.d))
            for k, o in enumerate(self.occ):
                if o[m] == 0: continue
                o2 = list(o); o2[m] -= 1; o2 = tuple(o2)
                if o2 in self.idx: X[self.idx[o2], k] = np.sqrt(o[m])
            self.b.append(X); self.bd.append(X.T.copy())
        self.P = {n: np.diag((self.N == n).astype(float)) for n in (1, 2)}
        if M == 4:
            self.lab = np.array([(sum(o[a]*GM.MODES_W[a][0] for a in range(4)), sum(o[a]*GM.MODES_W[a][1] for a in range(4))) for o in self.occ])
        else:
            self.lab = np.array([(o[0] - o[1],) for o in self.occ])
        self.Sz = np.diag(self.lab[:, 0]/2.0)
        self.Iz = np.diag(self.lab[:, 1]/2.0) if M == 4 else None
        self.nh = np.diag((self.N == 1).astype(float)) if kind == 'mob' else None

P11 = 'zero'
P11_CORE = 1.0          # rate units (about 104 J0): a hard core

def p11_coef(P, g, branch):
    """Coefficient of A^+A on a (1,1) link: the exact branch energy / |P|^2."""
    c2 = float((P*P).sum()); d = DELTA_T[(1, 1)]
    root = np.sqrt(d*d/4 + g*g*c2)
    return (d/2 + root)/c2 if branch == 'upper' else (d/2 - root)/c2

def link_terms(Li, Lj, P, g, p11=None):
    """List of (coef, op_i, op_j) for the exchange on a link, by type blocks."""
    M = P.shape[0]
    p11 = P11 if p11 is None else p11
    out = []
    for ni in (2, 1):
        if ni not in Li.N: continue
        for nj in (2, 1):
            if nj not in Lj.N: continue
            if Li.kind == 'end' and ni == 2: continue
            if Lj.kind == 'end' and nj == 2: continue
            if ni == 1 and nj == 1:
                if p11 == 'zero': continue
                if p11 == 'core':
                    out.append((P11_CORE, Li.P[1], Lj.P[1])); continue
                c = P11_CORE/float((P*P).sum()) if p11 == 'score' else p11_coef(P, g, p11)
            else:
                c = -g*g/DELTA_T[(ni, nj)]
            for (r, m) in itertools.product(range(M), repeat=2):
                for (s, n) in itertools.product(range(M), repeat=2):
                    w = P[r, s]*P[m, n]
                    if w == 0: continue
                    out.append((c*w, Li.P[ni]@Li.E[(r, m)], Lj.P[nj]@Lj.E[(s, n)]))
    return out

def hop_terms(Li, Lj, t):
    if t == 0 or Li.kind != 'mob' or Lj.kind != 'mob': return []
    out = []
    for m in range(Li.M):
        out.append((t, Li.b[m], Lj.bd[m]))     # quantum i -> j
        out.append((t, Li.bd[m], Lj.b[m]))     # quantum j -> i
    return out

def two_site_matrix(terms):
    return sum(c*np.kron(a, b) for c, a, b in terms)

def chain_terms(locs, P, g, t=0.0, ring=False):
    n = len(locs); T = []
    for i in range(n if ring else n - 1):
        j = (i + 1) % n
        T.append(((i, j), two_site_matrix(link_terms(locs[i], locs[j], P, g) + hop_terms(locs[i], locs[j], t))))
    return T

def dmrg_terms(locs, P, g, t=0.0):
    T = []
    for i in range(len(locs) - 1):
        for (c, a, b) in link_terms(locs[i], locs[i + 1], P, g) + hop_terms(locs[i], locs[i + 1], t):
            T.append((i, [a, b], c))
    return T

def sector(locs, terms, target, labfun):
    """Exact H in the sector where sum of labfun(loc) over sites == target."""
    n = len(locs); dims = np.array([L.d for L in locs])
    strides = np.ones(n, dtype=np.int64)
    for i in range(n - 2, -1, -1): strides[i] = strides[i + 1]*dims[i + 1]
    D = int(np.prod(dims)); codes = np.arange(D, dtype=np.int64)
    tot = None
    for i in range(n):
        li = labfun(locs[i])[(codes//strides[i]) % dims[i]]
        tot = li if tot is None else tot + li
    keep = np.all(tot == np.array(target)[None, :], axis=1)
    basis = codes[keep]; del codes, tot
    Db = len(basis)
    rows, cols, vals = [], [], []
    for (sites, Mx) in terms:
        ds = [dims[s] for s in sites]
        digs = [(basis//strides[s]) % dims[s] for s in sites]
        lo = np.zeros(Db, dtype=np.int64)
        for k, s in enumerate(sites): lo = lo*ds[k] + digs[k]
        for (a, b) in np.argwhere(np.abs(Mx) > 1e-14):
            sel = np.where(lo == b)[0]
            if len(sel) == 0: continue
            da = []; x = a
            for k in range(len(sites) - 1, -1, -1): da.append(x % ds[k]); x //= ds[k]
            da = da[::-1]
            new = basis[sel].copy()
            for k, s in enumerate(sites): new += (da[k] - digs[k][sel])*strides[s]
            pos = np.searchsorted(basis, new)
            ok = (pos < Db) & (basis[np.minimum(pos, Db - 1)] == new)
            rows.append(pos[ok]); cols.append(sel[ok]); vals.append(np.full(ok.sum(), Mx[a, b]))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(Db, Db)) if vals else sp.csr_matrix((Db, Db))
    return H, basis, dims, strides

def labfun_factory(M, mobile):
    def f(L):
        lab = L.lab.copy()
        if mobile:
            nh = (L.N == 1).astype(int) if L.kind == 'mob' else np.zeros(L.d, int)
            lab = np.concatenate([lab, nh[:, None]], 1)
        return lab
    return f
