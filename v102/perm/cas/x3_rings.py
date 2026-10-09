#!/usr/bin/env python3
"""X3 (rings) — T: the mobile hole on exact leading-order rings.
M = 2 (3+2): n = 6, 8, 10;  M = 4 (W): n = 4, 6.  Every site 'mob'; one hole.
Per t (units of J0 = 2 g^2/Delta): the lowest levels in the hole sector by
multiplet; the band width of the n positional states of the ground multiplet;
the lowest higher multiplet against the threshold E0(hole) + gap(no hole).
usage: x3_rings.py M n"""
import sys, time
import numpy as np
from scipy.sparse.linalg import eigsh
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here)
import graphm as GM, lo

t0 = time.time()
M, n = int(sys.argv[1]), int(sys.argv[2])
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA
locs = [lo.Local(M, 'mob') for _ in range(n)]
lf = lo.labfun_factory(M, True)
def low(H, k):
    """Lowest k eigenvalues, robust to exact degeneracy: dense when the sector
    fits, else LOBPCG with a block of 2k (block methods keep degenerate copies)."""
    k = min(k, H.shape[0] - 1)
    if H.shape[0] <= 12000:
        return np.linalg.eigvalsh(H.toarray())[:k]
    from scipy.sparse.linalg import lobpcg
    rng = np.random.default_rng(3)
    X = rng.normal(size=(H.shape[0], 2*k))
    w, _ = lobpcg(H, X, largest=False, tol=1e-9, maxiter=2000)
    return np.sort(w)[:k]
print(f"# X3 rings: M={M} n={n}, one mobile hole, leading order, units J0 = {J0:.6f}")
for tJ in (0.0, 0.1, 1.0, 3.0):
    T = lo.chain_terms(locs, P, g, t=tJ*J0, ring=True)
    if M == 2:
        secs = {'1/2': (1, 1), '3/2': (3, 1)}; nos = {'0': (0, 0), '1': (2, 0)}
    else:
        secs = {'(1/2,1/2)': (1, 1, 1), '(3/2,1/2)': (3, 1, 1), '(1/2,3/2)': (1, 3, 1), '(3/2,3/2)': (3, 3, 1)}
        nos = {'0': (0, 0, 0), '1': (2, 0, 0)}
    lev = {}
    for name, tg in secs.items():
        H, *_ = lo.sector(locs, T, tg, lf)
        lev[name] = low(H, n + 4)/J0
    if tJ == 0.0:
        e0n = low(lo.sector(locs, T, nos['0'], lf)[0], 1)[0]/J0
        e1n = low(lo.sector(locs, T, nos['1'], lf)[0], 1)[0]/J0
        gapn = e1n - e0n
    names = list(secs)
    base = lev[names[0]]
    E0 = base[0]
    # levels of the base sector that recur in a higher sector belong to a higher multiplet
    hi = np.concatenate([lev[k] for k in names[1:]])
    own = [e for e in base if not np.any(np.abs(hi - e) < 1e-7)]
    band = own[:n]
    W = band[-1] - band[0] if len(band) == n else float('nan')
    print(f"  t={tJ:4.1f} J0: ground {E0:+.5f} J0 in {names[0]}; own-multiplet levels (E-E0): " +
          " ".join(f"{e-E0:.4f}" for e in own[:n + 1]) + f";  band width of {n} positional states W = {W:.5f} J0")
    for k in names[1:]:
        print(f"           lowest {k}: {lev[k][0]-E0:+.4f} J0   (threshold hole + one quantum: {gapn:+.4f} J0;  "
              f"{'BOUND' if lev[k][0]-E0 < gapn - 1e-6 else 'above'})")
    print(f"           [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]")
