#!/usr/bin/env python3
"""The S-law with M field modes per cell (generalizes ../../options/cas/graph.py,
which has M = 2).  Used for the 3+4 widening (W): M = 4 modes per quantum,
labelled mu = (s, tau) with s, tau in {up, down} (spin x isospin).

State = (cells, cs): cells[i] = tuple of M occupations; cs[l] = space quanta.
    H = sum_i E_f(N_i) + sum_l Omega c_l^+ c_l + g sum_l (c_l^+ A_l + h.c.)
        + t sum_{l} sum_mu (b_i mu^+ b_j mu + h.c.)   [T2 rule: refused into a cell at cap]
A_l = sum_{mu nu} P[mu, nu] b_i mu b_j nu.
For W, P = eps (x) eps (the pair is a singlet in spin and in isospin).
For M = 2, P = eps reproduces the 3+2 S-law.
"""
import sys, itertools
from collections import deque
import numpy as np
import scipy.sparse as sp
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../options/cas')
from graph import EF, OMEGA, DELTA

EPS = np.array([[0.0, 1.0], [-1.0, 0.0]])
P_W = np.kron(EPS, EPS)            # modes ordered (s, tau): (u,u), (u,d), (d,u), (d,d)
P_2 = EPS
MODES_W = [(+1, +1), (+1, -1), (-1, +1), (-1, -1)]   # (2 s_z, 2 tau_z) of each mode

class GraphM:
    def __init__(self, M, ncell, links, g, P, caps, t=None, Omega=OMEGA, cmax=None, EFt=EF):
        self.M, self.n, self.links, self.g, self.P = M, ncell, list(links), g, P
        self.caps = caps; self.L = len(links)
        self.t = [0.0]*self.L if t is None else list(t)
        self.Omega = Omega; self.cmax = cmax; self.EF = EFt
        self.pairs = [(a, b, P[a, b]) for a in range(M) for b in range(M) if P[a, b] != 0]

    def valid(self, cells):
        return all(sum(c) <= cap for c, cap in zip(cells, self.caps))

    def diag(self, s):
        cells, cs = s
        return sum(self.EF[sum(c)] for c in cells) + self.Omega*sum(cs)

    def offdiag(self, s):
        cells, cs = s; out = []; sc = sum(cs)
        for l, (i, j) in enumerate(self.links):
            for (a, b, p) in self.pairs:
                if cells[i][a] > 0 and cells[j][b] > 0 and (self.cmax is None or sc + 1 <= self.cmax):
                    nc = list(cells)
                    x = list(nc[i]); x[a] -= 1; nc[i] = tuple(x)
                    y = list(nc[j]); y[b] -= 1; nc[j] = tuple(y)
                    cs2 = list(cs); cs2[l] += 1
                    out.append(((tuple(nc), tuple(cs2)), self.g*p*np.sqrt(cells[i][a]*cells[j][b]*(cs[l] + 1))))
                if cs[l] > 0:
                    nc = list(cells)
                    x = list(nc[i]); x[a] += 1; nc[i] = tuple(x)
                    y = list(nc[j]); y[b] += 1; nc[j] = tuple(y)
                    nc = tuple(nc)
                    if self.valid(nc):
                        cs2 = list(cs); cs2[l] -= 1
                        out.append(((nc, tuple(cs2)), self.g*p*np.sqrt((cells[i][a] + 1)*(cells[j][b] + 1)*cs[l])))
            tl = self.t[l]
            if tl != 0:
                for (src, dst) in ((i, j), (j, i)):
                    for a in range(self.M):
                        if cells[src][a] == 0: continue
                        nc = list(cells)
                        x = list(nc[src]); x[a] -= 1; nc[src] = tuple(x)
                        y = list(nc[dst]); y[a] += 1; nc[dst] = tuple(y)
                        nc = tuple(nc)
                        if not self.valid(nc): continue
                        out.append(((nc, cs), tl*np.sqrt(cells[src][a]*(cells[dst][a] + 1))))
        return out

    def build(self, seeds, keep=None):
        idx = {}; states = []; dq = deque()
        for s in seeds:
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
                rows.append(k2); cols.append(k); vals.append(amp)
        D = len(states)
        H = sp.csr_matrix((vals, (rows, cols)), shape=(D, D)) + sp.diags([self.diag(s) for s in states])
        self.states, self.idx = states, idx
        return H


def cell_states(M, N):
    """All occupation tuples of M modes with N quanta, in a fixed order."""
    out = []
    for combo in itertools.combinations_with_replacement(range(M), N):
        occ = [0]*M
        for a in combo: occ[a] += 1
        out.append(tuple(occ))
    return sorted(out, reverse=True)

def charges_W(cells):
    """(2 S^z, 2 I^z) of a W state."""
    s2 = sum(c[a]*MODES_W[a][0] for c in cells for a in range(4))
    t2 = sum(c[a]*MODES_W[a][1] for c in cells for a in range(4))
    return s2, t2
