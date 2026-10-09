#!/usr/bin/env python3
"""Effective chains for M-mode cells (3+2 with M = 2, the widening W with M = 4),
by des Cloizeaux on 2- and 3-cell open clusters with the local law
N_i + sum c = cap_i (t = 0), as ../../lock/cas/clusters.py does for M = 2.

Local model-space basis of a cell at cap c: graphm.cell_states(M, c)
(all c-quanta occupations, sorted descending).  Local labels: the
conserved charges of each basis state (for W: 2S^z, 2I^z; for M = 2: 2S^z).

sector_hamiltonian: H of a chain/ring in one sector of the conserved
charges, by vectorized digit arithmetic (as ../../fabric/cas/eff.py)."""
import sys, itertools
import numpy as np
import scipy.sparse as sp
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
import graphm as GM

def labels_of(M, occ):
    if M == 4:
        return (sum(occ[a]*GM.MODES_W[a][0] for a in range(4)), sum(occ[a]*GM.MODES_W[a][1] for a in range(4)))
    return (occ[0] - occ[1],)

def cluster(M, P, caps, g):
    """Exact effective Hamiltonian on the model space (all c = 0, every cell at cap)."""
    k = len(caps); nl = k - 1
    gm = GM.GraphM(M, k, [(i, i + 1) for i in range(nl)], g, P, caps)
    seeds = []
    for combo in itertools.product(*[GM.cell_states(M, c) for c in caps]):
        seeds.append((tuple(combo), (0,)*nl))
    H = gm.build(seeds).toarray()
    E0 = sum(GM.EF[c] for c in caps)
    H = H - E0*np.eye(len(gm.states))
    w, V = np.linalg.eigh(H)
    model = [gm.idx[(tuple(combo), (0,)*nl)] for combo in itertools.product(*[GM.cell_states(M, c) for c in caps])]
    m = len(model)
    X = V[model, :m]
    ev, U = np.linalg.eigh(X@X.T)
    Sm = U@np.diag(ev**-0.5)@U.T
    Heff = Sm@X@np.diag(w[:m])@X.T@Sm
    return 0.5*(Heff + Heff.T), w[m] - w[m - 1], float(ev.min())

def effective(M, P, g, pairs=((2, 2), (1, 2), (2, 1)), triples=((2, 2, 2), (1, 2, 2), (2, 2, 1), (2, 1, 2))):
    h2, info = {}, {}
    for pr in pairs:
        h, sep, evm = cluster(M, P, list(pr), g); h2[pr] = h; info[pr] = (sep, evm)
    r3 = {}
    for tr in triples:
        h3, sep, evm = cluster(M, P, list(tr), g)
        d = [len(GM.cell_states(M, c)) for c in tr]
        r3[tr] = h3 - np.kron(h2[(tr[0], tr[1])], np.eye(d[2])) - np.kron(np.eye(d[0]), h2[(tr[1], tr[2])])
        info[tr] = (sep, evm)
    return h2, r3, info

def chain_terms(caps, h2, r3, ring=False, with_r3=True):
    n = len(caps); T = []
    for i in range(n if ring else n - 1):
        j = (i + 1) % n
        T.append(((i, j), h2[(caps[i], caps[j])]))
    if with_r3:
        for i in range(n if ring else n - 2):
            j, k = (i + 1) % n, (i + 2) % n
            T.append(((i, j, k), r3[(caps[i], caps[j], caps[k])]))
    return T

def sector_hamiltonian(M, caps, terms, target):
    """target: tuple of conserved charges (2S^z[, 2I^z])."""
    loc = [GM.cell_states(M, c) for c in caps]
    dims = np.array([len(x) for x in loc]); n = len(caps)
    lab = [np.array([labels_of(M, o) for o in x]) for x in loc]      # (d_i, k)
    D = int(np.prod(dims))
    strides = np.ones(n, dtype=np.int64)
    for i in range(n - 2, -1, -1): strides[i] = strides[i + 1]*dims[i + 1]
    codes = np.arange(D, dtype=np.int64)
    tot = np.zeros((D, len(target)))
    for i in range(n):
        tot += lab[i][(codes//strides[i]) % dims[i]]
    keep = np.all(tot == np.array(target)[None, :], axis=1)
    basis = codes[keep]; Db = len(basis)
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
            for k in range(len(sites) - 1, -1, -1):
                da.append(x % ds[k]); x //= ds[k]
            da = da[::-1]
            new = basis[sel].copy()
            for k, s in enumerate(sites): new += (da[k] - digs[k][sel])*strides[s]
            pos = np.searchsorted(basis, new)
            ok = (pos < Db) & (basis[np.minimum(pos, Db - 1)] == new)
            rows.append(pos[ok]); cols.append(sel[ok]); vals.append(np.full(ok.sum(), Mx[a, b]))
    if not vals: return sp.csr_matrix((Db, Db)), basis, dims, strides, lab
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(Db, Db))
    return H, basis, dims, strides, lab

def profile(v, caps, basis, dims, strides, lab, comp=0):
    out = []
    for i in range(len(caps)):
        dig = (basis//strides[i]) % dims[i]
        out.append(float(np.abs(v)**2@(lab[i][dig, comp]/2)))
    return np.array(out)
