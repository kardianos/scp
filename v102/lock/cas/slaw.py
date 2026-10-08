#!/usr/bin/env python3
"""The S-law (../PREREG.md §2) as an exact sparse Hamiltonian on rings.

State: (cs, cells).  cs[l] = space quanta on link l = (l, l+1 mod n);
cells[i] = (a, b) = quanta in Jones modes 1 and 2 of cell i.
    H = sum_i E_f(N_i) + Omega sum_l c_l^+ c_l + g sum_l (c_l^+ A_l + h.c.)
        [+ t sum_a (b_ia^+ b_ja + h.c.), cap-refused; only for L5b]
A_l = b_i1 b_j2 - b_i2 b_j1.  The constant n*E_f(2) is subtracted, so the
all-field states sit at 0.  With t = 0 the basis is built inside the local
constraint N_i + c_{i-1} + c_i = 2 (PREREG algebra item 1); with t > 0 the
basis is the full Q = 2n, cap-2 space.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq

# ---------------------------------------------------------------- table constants
W1, W2, QD, CAP, A0 = 1.65, 2.9, 1.2, 2.5, 1.15
HBAR = A0/(2*np.pi)
OMEGA = W2

def _Ef(N):
    if N == 0: return 0.0
    return brentq(lambda e: e - N*W1/(1 + QD*HBAR*e/CAP), 0, N*W1 + 1)

EF = np.array([_Ef(N) for N in range(3)])
DELTA = OMEGA - 2*(EF[2] - EF[1])          # first conversion on a ring, 0.4793


# ---------------------------------------------------------------- basis
def basis_ring(n, cmax=None, sz=None, local=True):
    """List of states and index dict.  sz: 2*S^z target (integer) or None."""
    cmax = 2*n if cmax is None else cmax
    states = []
    if local:
        for cs in itertools.product(range(3), repeat=n):
            if sum(cs) > cmax: continue
            Ns = [2 - cs[(i - 1) % n] - cs[i] for i in range(n)]
            if min(Ns) < 0: continue
            for cells in itertools.product(*[[(a, N - a) for a in range(N + 1)] for N in Ns]):
                if sz is not None and sum(a - b for a, b in cells) != sz: continue
                states.append((cs, cells))
    else:
        by_N = {N: [(a, N - a) for a in range(N + 1)] for N in range(3)}
        for cs in itertools.product(range(3), repeat=n):
            sc = sum(cs)
            if sc > cmax: continue
            rem = 2*n - 2*sc
            for Ns in itertools.product(range(3), repeat=n):
                if sum(Ns) != rem: continue
                for cells in itertools.product(*[by_N[N] for N in Ns]):
                    if sz is not None and sum(a - b for a, b in cells) != sz: continue
                    states.append((cs, cells))
    idx = {s: k for k, s in enumerate(states)}
    return states, idx


def hamiltonian(n, g, states, idx, t=0.0, Omega=OMEGA, EFt=EF):
    rows, cols, vals = [], [], []
    E2 = EFt[2]
    for k, (cs, cells) in enumerate(states):
        rows.append(k); cols.append(k)
        vals.append(sum(EFt[a + b] for a, b in cells) - n*E2 + Omega*sum(cs))
        for l in range(n):
            i, j = l, (l + 1) % n
            for (ai, aj, sgn) in ((0, 1, +1), (1, 0, -1)):          # c_l^+ A_l
                if cells[i][ai] == 0 or cells[j][aj] == 0: continue
                amp = sgn*g*np.sqrt(cells[i][ai]*cells[j][aj]*(cs[l] + 1))
                nc = list(cells)
                x = list(nc[i]); x[ai] -= 1; nc[i] = tuple(x)
                y = list(nc[j]); y[aj] -= 1; nc[j] = tuple(y)
                cs2 = list(cs); cs2[l] += 1
                k2 = idx.get((tuple(cs2), tuple(nc)))
                if k2 is None: continue
                rows += [k2, k]; cols += [k, k2]; vals += [amp, amp]
            if t != 0:
                for (src, dst) in ((i, j), (j, i)):
                    if sum(cells[dst]) >= 2: continue
                    for a in (0, 1):
                        if cells[src][a] == 0: continue
                        amp = t*np.sqrt(cells[src][a]*(cells[dst][a] + 1))
                        nc = list(cells)
                        x = list(nc[src]); x[a] -= 1; nc[src] = tuple(x)
                        y = list(nc[dst]); y[a] += 1; nc[dst] = tuple(y)
                        k2 = idx.get((cs, tuple(nc)))
                        if k2 is None: continue
                        rows.append(k2); cols.append(k); vals.append(amp)
    D = len(states)
    return sp.csr_matrix((vals, (rows, cols)), shape=(D, D))


# ---------------------------------------------------------------- diagonal observables
def diag_obs(n, states):
    """Arrays (D, n): c_l, N_i, Sz_i."""
    D = len(states)
    C = np.zeros((D, n)); N = np.zeros((D, n)); Z = np.zeros((D, n))
    for k, (cs, cells) in enumerate(states):
        C[k] = cs
        N[k] = [a + b for a, b in cells]
        Z[k] = [(a - b)/2 for a, b in cells]
    return C, N, Z


def splus(i, states_from, idx_to):
    """Matrix of S^+_i from one basis to another (S^+ = b_1^+ b_2)."""
    rows, cols, vals = [], [], []
    for k, (cs, cells) in enumerate(states_from):
        a, b = cells[i]
        if b == 0: continue
        nc = list(cells); nc[i] = (a + 1, b - 1)
        k2 = idx_to.get((cs, tuple(nc)))
        if k2 is None: continue
        rows.append(k2); cols.append(k); vals.append(np.sqrt((a + 1)*b))
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(idx_to), len(states_from)))


def total_S(v, states, idx):
    """<S^2> -> S for a state of definite S^z = 0 (S^2 = S^-S^+ there)."""
    acc = {}
    for k, (cs, cells) in enumerate(states):
        if abs(v[k]) < 1e-13: continue
        for i, (a, b) in enumerate(cells):
            if b == 0: continue
            nc = list(cells); nc[i] = (a + 1, b - 1)
            key = (cs, tuple(nc))
            acc[key] = acc.get(key, 0.0) + v[k]*np.sqrt(b*(a + 1))
    s2 = sum(abs(x)**2 for x in acc.values())
    return -0.5 + np.sqrt(0.25 + s2)


# ---------------------------------------------------------------- product initial states
def product_state(cellvecs, states, idx, n):
    """cellvecs[i] = dict {(a,b): amplitude} on N=2; all c = 0."""
    v = np.zeros(len(states), complex)
    for combo in itertools.product(*[list(cv.items()) for cv in cellvecs]):
        cells = tuple(c for c, _ in combo)
        amp = np.prod([x for _, x in combo])
        k = idx.get(((0,)*n, cells))
        if k is not None: v[k] += amp
    return v


def jones_cell(alpha, beta):
    return {(2, 0): alpha**2, (1, 1): np.sqrt(2)*alpha*beta, (0, 2): beta**2}


def neel_cells(n):
    return [{(2, 0): 1.0} if i % 2 == 0 else {(0, 2): 1.0} for i in range(n)]


def random_jones_cells(n, rng):
    out = []
    for _ in range(n):
        z = rng.normal(size=2) + 1j*rng.normal(size=2); z /= np.linalg.norm(z)
        out.append(jones_cell(z[0], z[1]))
    return out


# ---------------------------------------------------------------- A^+A on two full cells
def AdagA_two_full():
    """A^+A restricted to two cells with N = 2 each (9 x 9), via (1,1) cells."""
    full = [(a, 2 - a) for a in range(3)]
    half = [(a, 1 - a) for a in range(2)]
    F = [(x, y) for x in full for y in full]; Hh = [(x, y) for x in half for y in half]
    hidx = {s: k for k, s in enumerate(Hh)}
    A = np.zeros((len(Hh), len(F)))
    for k, (ci, cj) in enumerate(F):
        for (ai, aj, sgn) in ((0, 1, +1), (1, 0, -1)):
            if ci[ai] == 0 or cj[aj] == 0: continue
            x = list(ci); x[ai] -= 1; y = list(cj); y[aj] -= 1
            A[hidx[(tuple(x), tuple(y))], k] += sgn*np.sqrt(ci[ai]*cj[aj])
    return A.T@A


# ---------------------------------------------------------------- spin-1 Heisenberg ring
def heis_ring(n):
    Sz = sp.diags([1.0, 0.0, -1.0]); Sp = sp.csr_matrix(np.diag([np.sqrt(2), np.sqrt(2)], 1))
    at = lambda op, i: sp.kron(sp.kron(sp.identity(3**i), op), sp.identity(3**(n - i - 1)), format='csr')
    Z = [at(Sz, i) for i in range(n)]; P = [at(Sp, i) for i in range(n)]
    H = sp.csr_matrix((3**n, 3**n))
    for i in range(n):
        j = (i + 1) % n
        H = H + Z[i]@Z[j] + 0.5*(P[i]@P[j].T + P[i].T@P[j])
    return H, Z, P


def fit_J(ex, hx):
    J = float(ex@hx/(hx@hx))
    return J, np.abs(ex/(J*hx) - 1)
