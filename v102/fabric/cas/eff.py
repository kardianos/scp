#!/usr/bin/env python3
"""Reading V: the fiber as the S-law effective chain (no real space states).

Sites carry spin cap/2 (cap 2: a full cell; cap 1: a hole or an R6 end
cell).  Terms are the des Cloizeaux h2 (pairs) and r3 (consecutive
triples) of ../../lock/cas/clusters.py, exact to fourth order in g, plus
the hole-centred triple (2,1,2) computed here the same way.

sector_hamiltonian builds H in one S^z sector by vectorized digit
arithmetic (basis = sorted integer codes; local site index = digit).
"""
import sys, itertools
import numpy as np
import scipy.sparse as sp
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../lock/cas')
from clusters import cluster, effective, SdotS
from slaw import DELTA

def terms_for(gr):
    """h2[(c1,c2)] and r3[(c1,c2,c3)] at g/Delta = gr, plus J."""
    g = gr*DELTA
    h2, r3, co, info = effective(g)
    h3, sep, ev = cluster([2, 1, 2], g)
    r3[(2, 1, 2)] = h3 - np.kron(h2[(2, 1)], np.eye(3)) - np.kron(np.eye(3), h2[(1, 2)])
    return h2, r3, co

def local_m(cap):
    """m values in the clusters' ordering: n1 = cap..0, m = (2 n1 - cap)/2."""
    return np.array([(2*n1 - cap)/2 for n1 in range(cap, -1, -1)])

def chain_terms(caps, h2, r3, ring=False, with_r3=True):
    """List of (sites, matrix) for a chain or ring of the given caps."""
    n = len(caps); T = []
    nb = n if ring else n - 1
    for i in range(nb):
        j = (i + 1) % n
        T.append(((i, j), h2[(caps[i], caps[j])]))
    if with_r3:
        nt = n if ring else n - 2
        for i in range(nt):
            j, k = (i + 1) % n, (i + 2) % n
            key = (caps[i], caps[j], caps[k])
            if key in r3: T.append(((i, j, k), r3[key]))
            else: raise KeyError(key)
    return T

def sector_basis(caps, twoSz):
    dims = np.array([c + 1 for c in caps]); n = len(caps)
    D = int(np.prod(dims))
    codes = np.arange(D, dtype=np.int64)
    strides = np.ones(n, dtype=np.int64)
    for i in range(n - 2, -1, -1): strides[i] = strides[i + 1]*dims[i + 1]
    tot = np.zeros(D)
    for i in range(n):
        dig = (codes//strides[i]) % dims[i]
        tot += local_m(caps[i])[dig]
    keep = codes[np.abs(2*tot - twoSz) < 1e-9]
    return keep, dims, strides

def sector_hamiltonian(caps, terms, twoSz):
    basis, dims, strides = sector_basis(caps, twoSz)
    D = len(basis)
    rows, cols, vals = [], [], []
    for (sites, M) in terms:
        ds = [dims[s] for s in sites]
        loc = np.zeros(D, dtype=np.int64); mult = 1
        digs = []
        for s in sites:
            digs.append((basis//strides[s]) % dims[s])
        # local index in kron order (first site most significant)
        for k, s in enumerate(sites):
            loc = loc*ds[k] + digs[k]
        nz = np.argwhere(np.abs(M) > 1e-15)
        for (a, b) in nz:
            sel = np.where(loc == b)[0]
            if len(sel) == 0: continue
            # digits of a
            da = []; x = a
            for k in range(len(sites) - 1, -1, -1):
                da.append(x % ds[k]); x //= ds[k]
            da = da[::-1]
            new = basis[sel].copy()
            for k, s in enumerate(sites):
                new += (da[k] - digs[k][sel])*strides[s]
            pos = np.searchsorted(basis, new)
            ok = (pos < D) & (basis[np.minimum(pos, D - 1)] == new)
            rows.append(pos[ok]); cols.append(sel[ok]); vals.append(np.full(ok.sum(), M[a, b]))
    if not vals:                     # e.g. the fully polarized sector: no conversion, H = 0 exactly
        return sp.csr_matrix((D, D)), basis, dims, strides
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
    return H, basis, dims, strides

def sz_profile(v, caps, basis, dims, strides):
    out = []
    for i, c in enumerate(caps):
        dig = (basis//strides[i]) % dims[i]
        out.append(float(np.abs(v)**2@local_m(c)[dig]))
    return np.array(out)

def total_S_from_levels(levels_by_sector):
    """levels_by_sector: {twoSz: sorted energies}. Assign S by multiplet counting."""
    raise NotImplementedError
