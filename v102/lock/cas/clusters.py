#!/usr/bin/env python3
"""Exact S-law cluster effective Hamiltonians (des Cloizeaux), as in
../../phase1/cas/p3_cluster.py, but with the S-law: loaded field energies
E_f(N), link oscillators at Omega, t = 0, and the local constraint
N_i + sum_{l at i} c_l = cap_i (R6 at t = 0).

Open cluster of k cells with caps (R6: cap = degree in the chain).
Model space: all c = 0, every cell full (spin cap/2).  H_eff on it,
ordered as the kron of m = S, S-1, ..., -S (n1 = cap, cap-1, ..., 0).
"""
import itertools, sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slaw import EF, OMEGA

def cluster(caps, g, Omega=OMEGA, EFt=EF):
    k = len(caps); nl = k - 1
    basis, idx = [], {}
    for cs in itertools.product(range(3), repeat=nl):
        Ns = []
        for i in range(k):
            used = (cs[i - 1] if i > 0 else 0) + (cs[i] if i < nl else 0)
            Ns.append(caps[i] - used)
        if min(Ns) < 0: continue
        for cells in itertools.product(*[[(a, N - a) for a in range(N + 1)] for N in Ns]):
            idx[(cs, cells)] = len(basis); basis.append((cs, cells))
    D = len(basis); H = np.zeros((D, D))
    E0 = sum(EFt[c] for c in caps)
    for q, (cs, cells) in enumerate(basis):
        H[q, q] = sum(EFt[a + b] for a, b in cells) - E0 + Omega*sum(cs)
        for l in range(nl):
            i, j = l, l + 1
            for (ai, aj, sgn) in ((0, 1, +1), (1, 0, -1)):
                if cells[i][ai] == 0 or cells[j][aj] == 0: continue
                amp = sgn*g*np.sqrt(cells[i][ai]*cells[j][aj]*(cs[l] + 1))
                nc = list(cells)
                x = list(nc[i]); x[ai] -= 1; nc[i] = tuple(x)
                y = list(nc[j]); y[aj] -= 1; nc[j] = tuple(y)
                cs2 = list(cs); cs2[l] += 1
                q2 = idx.get((tuple(cs2), tuple(nc)))
                if q2 is None: continue
                H[q2, q] += amp; H[q, q2] += amp
    w, V = np.linalg.eigh(H)
    model = []
    for n1s in itertools.product(*[range(cap, -1, -1) for cap in caps]):
        cells = tuple((n1, cap - n1) for n1, cap in zip(n1s, caps))
        model.append(idx[((0,)*nl, cells)])
    m = len(model)
    gapsep = w[m] - w[m - 1]
    X = V[model, :m]
    ev, U = np.linalg.eigh(X@X.T)
    Sm = U@np.diag(ev**-0.5)@U.T
    Heff = Sm@X@np.diag(w[:m])@X.T@Sm
    return 0.5*(Heff + Heff.T), gapsep, float(ev.min())


def spin_ops(s):
    m = np.arange(s, -s - 1, -1)
    Sz = np.diag(m)
    Sp = np.zeros((len(m), len(m)))
    for k in range(1, len(m)):
        Sp[k - 1, k] = np.sqrt(s*(s + 1) - m[k]*(m[k] + 1))
    return Sz, Sp

def SdotS(s1, s2):
    Z1, P1 = spin_ops(s1); Z2, P2 = spin_ops(s2)
    return np.kron(Z1, Z2) + 0.5*(np.kron(P1, P2.T) + np.kron(P1.T, P2))


# op_svd2 / op_svd3 copied verbatim from ../../phase1/cas/p3_cluster.py
# (that module parses argv at import)
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


def effective(g):
    """h2 for (2,2), (1,2), (2,1); r3 for (2,2,2), (1,2,2), (2,2,1)."""
    h2, info = {}, {}
    for pair in ((2, 2), (1, 2), (2, 1)):
        h, sep, evmin = cluster(list(pair), g)
        h2[pair] = h; info[pair] = (sep, evmin)
    SS = SdotS(1, 1); B = [np.eye(9), SS, SS@SS]
    A = np.stack([b.ravel() for b in B], 1)
    coef, *_ = np.linalg.lstsq(A, h2[(2, 2)].ravel(), rcond=None)
    resid = np.abs(h2[(2, 2)] - sum(c*b for c, b in zip(coef, B))).max()
    r3 = {}
    for trip in ((2, 2, 2), (1, 2, 2), (2, 2, 1)):
        h3, sep, evmin = cluster(list(trip), g)
        d = [c + 1 for c in trip]
        h12 = np.kron(h2[(trip[0], trip[1])], np.eye(d[2]))
        h23 = np.kron(np.eye(d[0]), h2[(trip[1], trip[2])])
        r3[trip] = h3 - h12 - h23
        info[trip] = (sep, evmin)
    a, J, K = coef
    # end bond: h2(1,2) = a' + J_end S(1/2).S(1)
    SSe = SdotS(0.5, 1); Be = [np.eye(6), SSe]
    Ae = np.stack([b.ravel() for b in Be], 1)
    ce, *_ = np.linalg.lstsq(Ae, h2[(1, 2)].ravel(), rcond=None)
    return h2, r3, dict(a=a, J=J, K=K, resid=resid, J_end=ce[1]), info
