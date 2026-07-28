#!/usr/bin/env python3
"""Fair race harness: the NAIVE algorithm, identical to phase_bell.c's, in the
natural numpy idiom (vectorised over samples, Python loop over settings)."""
import sys, time
import numpy as np
TAU = 2*np.pi
N    = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
grid = int(sys.argv[2]) if len(sys.argv) > 2 else 13
arm  = int(sys.argv[3]) if len(sys.argv) > 3 else 1
rng = np.random.default_rng(20260726)
lam = rng.uniform(0, TAU, N)
t0 = time.time()
best = 0.0
g = TAU*np.arange(grid)/grid
for ap in g:
    s0 = np.sign(np.cos(lam)); sa = np.sign(np.cos(lam-ap))
    for b in g:
        sb = -np.sign(np.cos(lam-b)); wb = np.abs(np.cos(lam-b)) if arm==2 else 1.0
        for bp in g:
            sbp = -np.sign(np.cos(lam-bp)); wbp = np.abs(np.cos(lam-bp)) if arm==2 else 1.0
            if arm==2:
                Wb, Wbp = wb.sum(), wbp.sum()
                S = ((wb*s0*sb).sum()/Wb + (wbp*s0*sbp).sum()/Wbp
                     + (wb*sa*sb).sum()/Wb - (wbp*sa*sbp).sum()/Wbp)
            else:
                S = ((s0*sb).mean() + (s0*sbp).mean() + (sa*sb).mean() - (sa*sbp).mean())
            if abs(S) > best: best = abs(S)
dt = time.time()-t0
print("PY  NAIVE  N=%d grid=%d arm=%d :  max|S| = %.5f   %8.3f s" % (N,grid,arm,best,dt))
