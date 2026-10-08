#!/usr/bin/env python3
"""The S-law on an arbitrary graph of cells, with the option rules of
../PREREG.md §2.

State = (cells, cs): cells[i] = (a, b) quanta in Jones modes 1, 2 of cell i;
cs[l] = space quanta on link l's oscillator.  The basis is the closure, under
H, of the seed states (breadth first), inside the filters (sum c <= cmax,
at most `outmax` quanta on cells outside `core`).

    H = sum_i E_f(N_i) + sum_l Omega_l c_l^+ c_l + onsite
        + g sum_{l in conv} (c_l^+ A_l + h.c.)        A_l = sum M_ab b_ia b_jb
        + sum_{l} t_l sum_a (b_ia^+ b_ja + h.c.)      [transport, by rule]
        + sum_{(l,m)} t_s (c_l^+ c_m + h.c.)          [space flow F1]

Ceiling rules: 'fixed' (caps list, None = no ceiling), 'geo' (cap = degree),
'held' (cap = number of links whose far cell holds field; None if zero).
A state is valid iff N_i <= cap_i wherever a cap exists; every process
lands only on valid states (so the ceiling refuses symmetrically).
Transport rules (symmetric by construction: checked in both the state and
its image, so H is Hermitian):
  'T2'  ceiling only;
  'T3'  neither end cell has N = cap, before or after;
  'T4'  neither end cell has N + sum_{l at cell} c_l = cap, before or after.
T0 is T2 with t = 0 on the links it shuts; T1 is T2 with a small t.
"""
import sys
from collections import deque
import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq

W1, W2, QD, CAPL, A0 = 1.65, 2.9, 1.2, 2.5, 1.15
HBAR = A0/(2*np.pi)
OMEGA = W2

def _Ef(N):
    if N == 0: return 0.0
    return brentq(lambda e: e - N*W1/(1 + QD*HBAR*e/CAPL), 0, N*W1 + 1)

EF = np.array([_Ef(N) for N in range(16)])
DELTA = OMEGA - 2*(EF[2] - EF[1])
WEDGE = np.array([[0.0, 1.0], [-1.0, 0.0]])        # A = b_i1 b_j2 - b_i2 b_j1


def twist_matrix(theta):
    """A_l -> b_i^T eps U(theta) b_j, U = exp(-i theta sigma_y/2) (real)."""
    c, s = np.cos(theta/2), np.sin(theta/2)
    U = np.array([[c, -s], [s, c]])
    # A = b_i1 (U b_j)_2 - b_i2 (U b_j)_1  =>  M[0, b] = U[1, b], M[1, b] = -U[0, b]
    return np.array([U[1], -U[0]])


class Graph:
    def __init__(self, ncell, links, g, conv=None, t=None, rule='T2', ceil='fixed',
                 caps=None, Omega=OMEGA, twist=None, shop=None, onsite=None,
                 cmax=None, core=None, outmax=None, EFt=EF):
        self.n = ncell; self.links = list(links); self.L = len(links)
        self.g = g
        self.conv = [True]*self.L if conv is None else list(conv)
        self.t = [0.0]*self.L if t is None else list(t)
        self.rule = rule; self.ceil = ceil
        self.caps = caps
        self.Omega = [Omega]*self.L if np.isscalar(Omega) else list(Omega)
        self.M = [WEDGE if (twist is None or l not in twist) else twist[l] for l in range(self.L)]
        self.shop = [] if shop is None else list(shop)
        self.onsite = onsite                     # function(N_i, i) -> complex, or None
        self.cmax = cmax; self.core = core; self.outmax = outmax
        self.EF = EFt
        self.nbr = [[] for _ in range(ncell)]     # (link, other cell)
        for l, (i, j) in enumerate(self.links):
            self.nbr[i].append((l, j)); self.nbr[j].append((l, i))
        self.deg = [len(x) for x in self.nbr]

    # ------------------------------------------------------------ ceilings
    def cap(self, cells, k):
        if self.ceil == 'fixed':
            return None if self.caps is None else self.caps[k]
        if self.ceil == 'geo':
            return self.deg[k]
        if self.ceil == 'held':
            m = sum(1 for (_, o) in self.nbr[k] if cells[o][0] + cells[o][1] > 0)
            return m if m > 0 else None
        raise ValueError(self.ceil)

    def valid(self, cells, cs=None):
        for k in range(self.n):
            cp = self.cap(cells, k)
            if cp is not None and cells[k][0] + cells[k][1] > cp: return False
        if self.outmax is not None:
            out = sum(cells[k][0] + cells[k][1] for k in range(self.n) if k not in self.core)
            if out > self.outmax: return False
        return True

    def frozen(self, cells, cs, k):
        cp = self.cap(cells, k)
        if cp is None: return False
        N = cells[k][0] + cells[k][1]
        if self.rule == 'T3': return N >= cp
        if self.rule == 'T4': return N + sum(cs[l] for (l, _) in self.nbr[k]) >= cp
        return False

    def hop_ok(self, s, s2, i, j):
        if self.rule == 'T2': return True
        for (cells, cs) in (s, s2):
            if self.frozen(cells, cs, i) or self.frozen(cells, cs, j): return False
        return True

    # ------------------------------------------------------------ H on one state
    def diag(self, s):
        cells, cs = s
        e = sum(self.EF[a + b] for a, b in cells) + sum(o*c for o, c in zip(self.Omega, cs))
        if self.onsite is not None:
            e = e + sum(self.onsite(a + b, k) for k, (a, b) in enumerate(cells))
        return e

    def offdiag(self, s):
        cells, cs = s
        out = []
        sc = sum(cs)
        for l, (i, j) in enumerate(self.links):
            if self.conv[l] and self.g != 0:
                M = self.M[l]
                for al in (0, 1):
                    for be in (0, 1):
                        m = M[al, be]
                        if m == 0: continue
                        # c^+ A
                        if cells[i][al] > 0 and cells[j][be] > 0 and (self.cmax is None or sc + 1 <= self.cmax):
                            nc = list(cells)
                            x = list(nc[i]); x[al] -= 1; nc[i] = tuple(x)
                            y = list(nc[j]); y[be] -= 1; nc[j] = tuple(y)
                            cs2 = list(cs); cs2[l] += 1
                            nc = tuple(nc)
                            if self.valid(nc):
                                out.append(((nc, tuple(cs2)), self.g*m*np.sqrt(cells[i][al]*cells[j][be]*(cs[l] + 1))))
                        # c A^+
                        if cs[l] > 0:
                            nc = list(cells)
                            x = list(nc[i]); x[al] += 1; nc[i] = tuple(x)
                            y = list(nc[j]); y[be] += 1; nc[j] = tuple(y)
                            cs2 = list(cs); cs2[l] -= 1
                            nc = tuple(nc)
                            if self.valid(nc):
                                out.append(((nc, tuple(cs2)), self.g*m*np.sqrt((cells[i][al] + 1)*(cells[j][be] + 1)*cs[l])))
            tl = self.t[l]
            if tl != 0:
                for (src, dst) in ((i, j), (j, i)):
                    for a in (0, 1):
                        if cells[src][a] == 0: continue
                        nc = list(cells)
                        x = list(nc[src]); x[a] -= 1; nc[src] = tuple(x)
                        y = list(nc[dst]); y[a] += 1; nc[dst] = tuple(y)
                        nc = tuple(nc)
                        if not self.valid(nc): continue
                        s2 = (nc, cs)
                        if not self.hop_ok(s, s2, src, dst): continue
                        out.append((s2, tl*np.sqrt(cells[src][a]*(cells[dst][a] + 1))))
        for (l1, l2, ts) in self.shop:
            for (p, q) in ((l1, l2), (l2, l1)):
                if cs[p] == 0: continue
                cs2 = list(cs); cs2[p] -= 1; cs2[q] += 1
                out.append(((cells, tuple(cs2)), ts*np.sqrt(cs[p]*(cs[q] + 1))))
        return out

    # ------------------------------------------------------------ basis + matrix
    def build(self, seeds, keep=None, maxdim=5_000_000):
        """keep(s) -> bool: extra filter (e.g. an S^z sector)."""
        idx = {}; states = []
        dq = deque()
        for s in seeds:
            if s[1] is not None and self.cmax is not None and sum(s[1]) > self.cmax: continue
            if not self.valid(s[0]): continue
            if s not in idx and (keep is None or keep(s)):
                idx[s] = len(states); states.append(s); dq.append(s)
        rows, cols, vals = [], [], []
        while dq:
            s = dq.popleft(); k = idx[s]
            for (s2, amp) in self.offdiag(s):
                if keep is not None and not keep(s2): continue
                k2 = idx.get(s2)
                if k2 is None:
                    k2 = len(states); idx[s2] = k2; states.append(s2); dq.append(s2)
                    if len(states) > maxdim: raise RuntimeError("basis too large")
                rows.append(k2); cols.append(k); vals.append(amp)
        D = len(states)
        dvals = [self.diag(s) for s in states]
        cplx = any(isinstance(v, complex) for v in dvals)
        H = sp.csr_matrix((vals, (rows, cols)), shape=(D, D), dtype=complex if cplx else float)
        H = H + sp.diags(np.array(dvals, dtype=complex if cplx else float))
        self.states, self.idx = states, idx
        return H


# ---------------------------------------------------------------- helpers
def ring_links(n, offset=0):
    return [(offset + i, offset + (i + 1) % n) for i in range(n)]

def full_cells(n):
    return [(a, 2 - a) for a in range(3)]

def sz_of(s):
    return sum(a - b for a, b in s[0])

def ring_seeds(n, L, ncell=None, sz=None):
    """All-field ring product states (ring cells 0..n-1 at N = 2), other cells
    empty, no space.  sz: 2 S^z filter or None."""
    import itertools
    ncell = n if ncell is None else ncell
    out = []
    for combo in itertools.product([(2, 0), (1, 1), (0, 2)], repeat=n):
        if sz is not None and sum(a - b for a, b in combo) != sz: continue
        cells = list(combo) + [(0, 0)]*(ncell - n)
        out.append((tuple(cells), (0,)*L))
    return out

def spin_ops(s):
    m = np.arange(s, -s - 1, -1)
    Sz = np.diag(m)
    Sp = np.zeros((len(m), len(m)))
    for k in range(1, len(m)):
        Sp[k - 1, k] = np.sqrt(s*(s + 1) - m[k]*(m[k] + 1))
    return Sz, Sp

def heis_spectrum(n, s=1, pbc=True, sz=0):
    """Low spectrum of the spin-s Heisenberg ring in the S^z = sz sector."""
    d = int(round(2*s + 1))
    Sz, Sp = spin_ops(s)
    at = lambda op, i: sp.kron(sp.kron(sp.identity(d**i), sp.csr_matrix(op)), sp.identity(d**(n - i - 1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((d**n, d**n))
    for i in range(n if pbc else n - 1):
        j = (i + 1) % n
        H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
    tot = sum(Z).diagonal()
    keep = np.where(np.abs(tot - sz) < 1e-9)[0]
    return np.linalg.eigvalsh(H[keep][:, keep].toarray())

def fit_J(ex, hx):
    J = float(ex@hx/(hx@hx))
    return J, np.abs(ex/(J*hx) - 1)


def total_S(v, states):
    """S for a state with S^z = 0 (S^2 = S^- S^+ there); S^+ = sum_i b_i1^+ b_i2."""
    acc = {}
    for k, (cells, cs) in enumerate(states):
        if abs(v[k]) < 1e-13: continue
        for i, (a, b) in enumerate(cells):
            if b == 0: continue
            nc = list(cells); nc[i] = (a + 1, b - 1)
            key = (tuple(nc), cs)
            acc[key] = acc.get(key, 0.0) + v[k]*np.sqrt(b*(a + 1))
    s2 = sum(abs(x)**2 for x in acc.values())
    return -0.5 + np.sqrt(0.25 + s2)


def obs_arrays(states, ncell, L):
    """Diagonal observables: N (D, ncell), C (D, L)."""
    D = len(states)
    N = np.zeros((D, ncell)); C = np.zeros((D, L))
    for k, (cells, cs) in enumerate(states):
        N[k] = [a + b for a, b in cells]; C[k] = cs
    return N, C
