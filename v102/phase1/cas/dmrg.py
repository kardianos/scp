#!/usr/bin/env python3
"""v102 Phase 1 — a small two-site DMRG (numpy/scipy only).

MPS  A[i]: (Dl, d, Dr).  MPO W[i]: (wl, wr, d_out, d_in).
build_mpo() takes on-site operators, terms on consecutive sites, and
quadratic penalties lam*(sum_i q_i - target)^2 (used to pin S^z_tot).
Real arithmetic throughout (every model here is real).
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


def build_mpo(dims, onsite, terms, penalties=()):
    """dims: local dims.  onsite: list of (d,d) or None.
    terms: list of (start, [op0, op1, ...], coef) on consecutive sites.
    penalties: list of (qops, lam, target); qops: list of (d,d) or None.
    Returns (W list, constant)."""
    N = len(dims)
    # channels per bond b (between site b and b+1): 0=start, 1=done, then extras
    chans = [dict() for _ in range(N - 1)]
    nch = [2]*(N - 1)
    def ch(b, key):
        if key not in chans[b]:
            chans[b][key] = nch[b]; nch[b] += 1
        return chans[b][key]
    for k, (s, ops, coef) in enumerate(terms):
        for m in range(len(ops) - 1):
            ch(s + m, ('t', k, m))
    for p in range(len(penalties)):
        for b in range(N - 1):
            ch(b, ('p', p))
    const = 0.0
    W = []
    for i in range(N):
        d = dims[i]
        wl = 1 if i == 0 else nch[i - 1]
        wr = 1 if i == N - 1 else nch[i]
        Wi = np.zeros((wl, wr, d, d))
        L_start = 0
        L_done = None if i == 0 else 1
        R_start = None if i == N - 1 else 0
        R_done = 0 if i == N - 1 else 1
        I = np.eye(d)
        if R_start is not None: Wi[L_start, R_start] += I
        if L_done is not None: Wi[L_done, R_done] += I
        h = np.zeros((d, d)) if onsite[i] is None else onsite[i].copy()
        for (qops, lam, target) in penalties:
            q = qops[i]
            if q is not None:
                h += lam*(q@q) - 2*lam*target*q
        Wi[L_start, R_done] += h
        for p, (qops, lam, target) in enumerate(penalties):
            q = qops[i] if qops[i] is not None else np.zeros((d, d))
            if R_start is not None:
                Wi[L_start, chans[i][('p', p)]] += 2*lam*q
            if i > 0 and R_start is not None:
                Wi[chans[i - 1][('p', p)], chans[i][('p', p)]] += I
            if i > 0:
                Wi[chans[i - 1][('p', p)], R_done] += q
        for k, (s, ops, coef) in enumerate(terms):
            n = len(ops)
            if not (s <= i < s + n): continue
            m = i - s
            if n == 1:
                Wi[L_start, R_done] += coef*ops[0]
            elif m == 0:
                Wi[L_start, chans[i][('t', k, 0)]] += coef*ops[0]
            elif m < n - 1:
                Wi[chans[i - 1][('t', k, m - 1)], chans[i][('t', k, m)]] += ops[m]
            else:
                Wi[chans[i - 1][('t', k, m - 1)], R_done] += ops[m]
        W.append(Wi)
    for (qops, lam, target) in penalties:
        const += lam*target*target
    return W, const


def _lenv(L, A, W):
    t = np.tensordot(L, A, axes=([0], [0]))            # w a' t b
    t = np.tensordot(t, W, axes=([0, 2], [0, 3]))      # a' b v s
    return np.tensordot(t, A, axes=([0, 3], [0, 1]))   # b v b'


def _renv(R, A, W):
    t = np.tensordot(A, R, axes=([2], [0]))            # a t v b'
    t = np.tensordot(t, W, axes=([1, 2], [3, 1]))      # a b' w s
    return np.tensordot(t, A, axes=([3, 1], [1, 2]))   # a w a'


def _heff(L, W1, W2, R, shape):
    def mv(x):
        x = x.reshape(shape)
        t = np.tensordot(L, x, axes=([0], [0]))        # w a' s1 s2 b
        t = np.tensordot(t, W1, axes=([0, 2], [0, 3])) # a' s2 b u s1'
        t = np.tensordot(t, W2, axes=([3, 1], [0, 3])) # a' b s1' v s2'
        t = np.tensordot(t, R, axes=([1, 3], [0, 1]))  # a' s1' s2' b'
        return t.ravel()
    return mv


class DMRG:
    def __init__(self, W, dims, chi=64, seed=0, init=None):
        self.W, self.dims, self.N = W, dims, len(dims)
        self.chi = chi
        rng = np.random.default_rng(seed)
        A = []
        Dl = 1
        for i, d in enumerate(dims):
            Dr = 1 if i == self.N - 1 else min(chi, 4)
            if init is not None:
                T = np.zeros((Dl, d, Dr)); T[0, init[i], 0] = 1.0
                T += 1e-3*rng.normal(size=T.shape)
            else:
                T = rng.normal(size=(Dl, d, Dr))
            A.append(T); Dl = Dr
        self.A = A
        self._right_canon()
        self.L = [None]*(self.N + 1); self.R = [None]*(self.N + 1)
        self.L[0] = np.ones((1, 1, 1)); self.R[self.N] = np.ones((1, 1, 1))
        for i in range(self.N - 1, 0, -1):
            self.R[i] = _renv(self.R[i + 1], self.A[i], self.W[i])

    def _right_canon(self):
        for i in range(self.N - 1, 0, -1):
            Dl, d, Dr = self.A[i].shape
            M = self.A[i].reshape(Dl, d*Dr)
            Q, Rm = np.linalg.qr(M.T)
            self.A[i] = Q.T.reshape(-1, d, Dr)
            self.A[i - 1] = np.tensordot(self.A[i - 1], Rm.T, axes=([2], [0]))
        self.A[0] /= np.linalg.norm(self.A[0])

    def _solve(self, i, tol):
        th = np.tensordot(self.A[i], self.A[i + 1], axes=([2], [0]))
        shape = th.shape
        n = th.size
        mv = _heff(self.L[i], self.W[i], self.W[i + 1], self.R[i + 2], shape)
        if n <= 400:
            H = np.column_stack([mv(e) for e in np.eye(n)])
            w, v = np.linalg.eigh(0.5*(H + H.T))
            return w[0], v[:, 0].reshape(shape)
        op = LinearOperator((n, n), matvec=mv, dtype=float)
        w, v = eigsh(op, k=1, which='SA', v0=th.ravel(), tol=tol, ncv=min(n, 20))
        return w[0], v[:, 0].reshape(shape)

    def sweep(self, tol=1e-10, cutoff=1e-12):
        E = None; trunc = 0.0
        for i in range(self.N - 1):                      # left -> right
            E, th = self._solve(i, tol)
            Dl, d1, d2, Dr = th.shape
            U, S, V = np.linalg.svd(th.reshape(Dl*d1, d2*Dr), full_matrices=False)
            keep = max(1, min(self.chi, int(np.sum(S > cutoff*S[0]))))
            trunc = max(trunc, float(np.sum(S[keep:]**2)))
            U, S, V = U[:, :keep], S[:keep]/np.linalg.norm(S[:keep]), V[:keep]
            self.A[i] = U.reshape(Dl, d1, keep)
            self.A[i + 1] = (np.diag(S)@V).reshape(keep, d2, Dr)
            self.L[i + 1] = _lenv(self.L[i], self.A[i], self.W[i])
        for i in range(self.N - 2, -1, -1):              # right -> left
            E, th = self._solve(i, tol)
            Dl, d1, d2, Dr = th.shape
            U, S, V = np.linalg.svd(th.reshape(Dl*d1, d2*Dr), full_matrices=False)
            keep = max(1, min(self.chi, int(np.sum(S > cutoff*S[0]))))
            trunc = max(trunc, float(np.sum(S[keep:]**2)))
            U, S, V = U[:, :keep], S[:keep]/np.linalg.norm(S[:keep]), V[:keep]
            self.A[i] = (U@np.diag(S)).reshape(Dl, d1, keep)
            self.A[i + 1] = V.reshape(keep, d2, Dr)
            self.R[i + 1] = _renv(self.R[i + 2], self.A[i + 1], self.W[i + 1])
        return E, trunc

    def run(self, sweeps=8, etol=1e-10, verbose=False):
        Eold = None
        for s in range(sweeps):
            E, tr = self.sweep()
            if verbose: print(f"    sweep {s}: E={E:.12f} trunc={tr:.1e}", flush=True)
            if Eold is not None and abs(E - Eold) < etol: break
            Eold = E
        return E

    # ---- measurements (MPS is right-canonical with the centre at site 0 after run())
    def _left_upto(self, i):
        X = np.ones((1, 1))
        for k in range(i):
            X = np.tensordot(X, self.A[k], axes=([0], [0]))          # a' s b
            X = np.tensordot(X, self.A[k], axes=([0, 1], [0, 1]))    # b b'
        return X

    def expect1(self, op, i):
        X = self._left_upto(i)
        t = np.tensordot(X, self.A[i], axes=([0], [0]))              # a' s b
        t = np.tensordot(t, op, axes=([1], [1]))                     # a' b s'
        return float(np.tensordot(t, self.A[i], axes=([0, 2, 1], [0, 1, 2])))

    def corr_row(self, opA, i, opB_list, js):
        """<opA_i opB_j> for every j in js (all > i), one left pass."""
        X = self._left_upto(i)
        t = np.tensordot(X, self.A[i], axes=([0], [0]))
        t = np.tensordot(t, opA, axes=([1], [1]))
        X = np.tensordot(t, self.A[i], axes=([0, 2], [0, 1]))        # b b'
        out = {}
        jmax = max(js)
        for k in range(i + 1, jmax + 1):
            if k in js:
                t = np.tensordot(X, self.A[k], axes=([0], [0]))
                t = np.tensordot(t, opB_list[k], axes=([1], [1]))
                out[k] = float(np.tensordot(t, self.A[k], axes=([0, 2, 1], [0, 1, 2])))
            X = np.tensordot(X, self.A[k], axes=([0], [0]))
            X = np.tensordot(X, self.A[k], axes=([0, 1], [0, 1]))
        return out


# ---------------------------------------------------------------- operators
def spin_ops(s):
    m = np.arange(s, -s - 1, -1); d = len(m)
    Sz = np.diag(m); Sp = np.zeros((d, d))
    for q in range(1, d):
        Sp[q - 1, q] = np.sqrt(s*(s + 1) - m[q]*(m[q] + 1))
    return Sz, Sp


def heisenberg_chain(spins, J=1.0, Sz_target=None, lam=1.0):
    """Open chain, nearest-neighbour J S.S, mixed spins allowed."""
    dims = [int(round(2*s + 1)) for s in spins]
    ops = [spin_ops(s) for s in spins]
    terms = []
    for i in range(len(spins) - 1):
        (Zi, Pi), (Zj, Pj) = ops[i], ops[i + 1]
        terms += [(i, [Zi, Zj], J), (i, [Pi, Pj.T], 0.5*J), (i, [Pi.T, Pj], 0.5*J)]
    pen = []
    if Sz_target is not None:
        pen = [([o[0] for o in ops], lam, Sz_target)]
    W, c = build_mpo(dims, [None]*len(dims), terms, pen)
    return W, c, dims, ops
