#!/usr/bin/env python3
"""Compact exact MPO for a chain Hamiltonian of two-site terms h2[i] on
(i, i+1) and three-site terms r3[i] on (i, i+1, i+2), plus the quadratic
S^z penalty of ../../phase1/cas/dmrg.py build_mpo.

Channels on bond b (between sites b and b+1):
  0 start, 1 done,
  L1: d_b^2 channels  "a term began with E_pq on site b",
  L2: d_{b+1}^2 channels "a three-site term began at b-1; it waits for
      E_tu on site b+1",
  P : penalty carry.
So w <= 2 + 9 + 9 + 1 = 21, against ~156 for one channel per SVD term.
W[i]: (wl, wr, d_out, d_in), the dmrg.py convention."""
import numpy as np


def mpo_chain(dims, h2, r3, qops=None, lam=0.0, target=0.0):
    """h2[i]: (d_i d_{i+1})^2 matrix or None; r3[i]: (d_i d_{i+1} d_{i+2})^2 or None."""
    N = len(dims)
    def E(d, p, q):
        m = np.zeros((d, d)); m[p, q] = 1.0; return m
    # channel layout per bond
    lay = []
    for b in range(N - 1):
        d1 = dims[b]; d2 = dims[b + 1]
        off = 2
        L1 = off; off += d1*d1
        L2 = off; off += d2*d2 if b >= 1 else 0
        P = off; off += 1 if qops is not None else 0
        lay.append(dict(L1=L1, L2=L2, P=P, w=off))
    W = []
    for i in range(N):
        d = dims[i]; I = np.eye(d)
        wl = 1 if i == 0 else lay[i - 1]['w']
        wr = 1 if i == N - 1 else lay[i]['w']
        Wi = np.zeros((wl, wr, d, d))
        Ls, Ld = 0, (None if i == 0 else 1)
        Rs, Rd = (None if i == N - 1 else 0), (0 if i == N - 1 else 1)
        if Rs is not None: Wi[Ls, Rs] += I
        if Ld is not None: Wi[Ld, Rd] += I
        h = np.zeros((d, d))
        if qops is not None:
            q = qops[i]; h += lam*(q@q) - 2*lam*target*q
        Wi[Ls, Rd] += h
        # penalty carry
        if qops is not None:
            q = qops[i]
            if Rs is not None: Wi[Ls, lay[i]['P']] += 2*lam*q
            if i > 0 and Rs is not None: Wi[lay[i - 1]['P'], lay[i]['P']] += I
            if i > 0: Wi[lay[i - 1]['P'], Rd] += q
        # open a term at site i: E_pq into L1 of bond i
        if Rs is not None and (h2[i] is not None or (i < N - 2 and r3[i] is not None)):
            for p in range(d):
                for qq in range(d):
                    Wi[Ls, lay[i]['L1'] + p*d + qq] += E(d, p, qq)
        if i >= 1:
            dl = dims[i - 1]
            # close a two-site term (i-1, i)
            if h2[i - 1] is not None:
                h4 = h2[i - 1].reshape(dl, d, dl, d)
                for p in range(dl):
                    for qq in range(dl):
                        Wi[lay[i - 1]['L1'] + p*dl + qq, Rd] += h4[p, :, qq, :]
            # middle of a three-site term (i-1, i, i+1)
            if i <= N - 2 and r3[i - 1] is not None:
                dr = dims[i + 1]
                r6 = r3[i - 1].reshape(dl, d, dr, dl, d, dr)
                for p in range(dl):
                    for qq in range(dl):
                        for t in range(dr):
                            for u in range(dr):
                                Wi[lay[i - 1]['L1'] + p*dl + qq, lay[i]['L2'] + t*dr + u] += r6[p, :, t, qq, :, u]
        if i >= 2 and r3[i - 2] is not None:
            # close a three-site term (i-2, i-1, i)
            for t in range(d):
                for u in range(d):
                    Wi[lay[i - 1]['L2'] + t*d + u, Rd] += E(d, t, u)
        W.append(Wi)
    const = lam*target*target if qops is not None else 0.0
    return W, const


def mpo_to_dense(W):
    M = W[0][0]                                  # (wr, d, d)
    for Wi in W[1:]:
        wr = Wi.shape[1]
        # M: (w, D, D), Wi: (w, wr, d, d) -> (wr, D d, D d)
        T = np.einsum('aij,abkl->bikjl', M, Wi)
        D = M.shape[1]*Wi.shape[2]
        M = T.reshape(wr, D, D)
    return M[0]
