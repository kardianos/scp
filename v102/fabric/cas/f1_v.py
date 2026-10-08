#!/usr/bin/env python3
"""F1 — reading V against the microscopic S-law (../PREREG.md §4 F1).
Rings n = 6, 8 (no hole, S^z = 0) and the one-hole ring n = 6 (S^z = 1/2):
the effective chain (h2 + r3) against the exact S-law (link holds real
quanta).  g/Delta = 0.1."""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../options/cas')
import eff
import graph as G

t0 = time.time()
gr = 0.1; g = gr*G.DELTA
h2, r3, co = eff.terms_for(gr)
J = co['J']
print(f"# F1 reading V (effective, no real link states) against the S-law, g/Delta = {gr}, J_h2 = {J:.6f}")

def micro(n, hole, cmax):
    caps = [2]*n
    import itertools
    seeds = []
    for combo in itertools.product([(2, 0), (1, 1), (0, 2)], repeat=n):
        cells = list(combo)
        if hole: cells[0] = (1, 0)
        seeds.append((tuple(cells), (0,)*n))
    target = 1 if hole else 0
    gm = G.Graph(n, G.ring_links(n), g, ceil='fixed', caps=caps, cmax=cmax)
    H = gm.build(seeds, keep=lambda s: G.sz_of(s) == target)
    E2 = sum(G.EF[2] for _ in range(n)) - (G.EF[2] - G.EF[1] if hole else 0)
    H = H - E2*sp.identity(H.shape[0])
    return H, gm

for (n, hole, cmax) in ((6, False, None), (8, False, 3), (6, True, None), (8, True, 3)):
    Hm, gm = micro(n, hole, cmax)
    caps = [1 if (hole and i == 0) else 2 for i in range(n)]
    He, basis, dims, strides = eff.sector_hamiltonian(caps, eff.chain_terms(caps, h2, r3, ring=True), 1 if hole else 0)
    k = 30
    wm = np.sort(eigsh(Hm, k=k, which='SA', tol=1e-12)[0]) if Hm.shape[0] > 3000 else np.linalg.eigvalsh(Hm.toarray())[:k]
    we = np.sort(eigsh(He, k=k, which='SA', tol=1e-12)[0]) if He.shape[0] > 3000 else np.linalg.eigvalsh(He.toarray())[:k]
    # the two spectra share a constant only up to the effective a-term; compare excitations
    dm = wm - wm[0]; de = we - we[0]
    nspin = He.shape[0]
    print(f"  n={n} hole={hole} cmax={cmax}: micro dim {Hm.shape[0]} vs V dim {nspin} (V removes {Hm.shape[0]-nspin} space-carrying states); "
          f"E0 micro {wm[0]:+.6f} V {we[0]:+.6f} (diff {we[0]-wm[0]:+.2e} = {(we[0]-wm[0])/J:+.3f} J); "
          f"lowest {k} excitations: max |diff| = {np.abs(dm-de).max():.2e} = {np.abs(dm-de).max()/J:.4f} J, "
          f"relative to the first gap {np.abs(dm-de).max()/dm[1]:.4f}  [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
