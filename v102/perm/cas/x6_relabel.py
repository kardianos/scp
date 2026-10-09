#!/usr/bin/env python3
"""X6 — the relabelling checks.
P1 (I): the 3+2 hole ring's doublet: the two members (2S^z = +1, -1) are
    degenerate (that is the p/n pair under I) and their profiles mirror.
P3 (I + W = W): W's spectrum is symmetric under spin <-> isospin: compare
    the lowest levels of sector pairs (2,0)/(0,2) without a hole, and
    (3,1)/(1,3) with one hole, on effective rings n = 4 (all orders)."""
import sys, time
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../fabric/cas')
import graphm as GM, weff, eff
t0 = time.time(); g = 0.1*GM.DELTA
h2, r3, co = eff.terms_for(0.1)
caps = [1] + [2]*7
T = eff.chain_terms(caps, h2, r3, ring=True)
lev = {}
for tw in (1, -1):
    H, basis, dims, strides = eff.sector_hamiltonian(caps, T, tw)
    w, v = np.linalg.eigh(H.toarray())
    lev[tw] = (w[:6], eff.sz_profile(v[:, 0], caps, basis, dims, strides))
print(f"# X6 P1: 3+2 hole ring n=8: lowest levels 2S^z=+1: {np.round(lev[1][0],8)}; 2S^z=-1: {np.round(lev[-1][0],8)}; "
      f"max diff {np.abs(lev[1][0]-lev[-1][0]).max():.1e}; profiles mirror: max |p(+)+p(-)| = {np.abs(lev[1][1]+lev[-1][1]).max():.1e}")
h2w, r3w, info = weff.effective(4, GM.P_W, g)
for (name, caps, pairs) in (("no hole", [2]*4, [((2, 0), (0, 2)), ((4, 0), (0, 4))]), ("one hole", [1, 2, 2, 2], [((3, 1), (1, 3)), ((5, 1), (1, 5))])):
    Tw = weff.chain_terms(caps, h2w, r3w, ring=True)
    for a, b in pairs:
        wa = np.linalg.eigvalsh(weff.sector_hamiltonian(4, caps, Tw, a)[0].toarray())[:8]
        wb = np.linalg.eigvalsh(weff.sector_hamiltonian(4, caps, Tw, b)[0].toarray())[:8]
        print(f"# X6 P3: W ring n=4 {name}: sectors {a} vs {b}: max |diff| of the lowest 8 = {np.abs(wa-wb).max():.1e}")
print(f"done [{time.time()-t0:.0f}s]")
