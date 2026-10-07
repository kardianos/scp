#!/usr/bin/env python3
"""P1 diagnostic (post-registration, labelled as such in DECIDE):
compare levels WITH multiplicity, and give the spin and <space> of each
level, on the n = 4 ring.  It separates a split accidental degeneracy
of the 4-ring from a real change of ground state."""
import sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from p1_micro import micro, S2_sz0, spin1, sector
n = 4
H, szd = spin1(n)
wh, vh = np.linalg.eigh(sector(H, szd, 0).toarray())
print("Heisenberg n=4 ring, S^z=0 sector (J=1): E - E0 =", np.round(wh - wh[0], 4).tolist())
for t in (0.0, 0.05, 0.1, 0.15, 0.2):
    Hm, space, basis = micro(n, 0.1, t, cmax=4)
    w, v = np.linalg.eigh(Hm.toarray())
    ex = w[:19] - w[0]
    hx = wh[:19] - wh[0]
    J = float(ex[1:]@hx[1:]/(hx[1:]@hx[1:]))
    Ss = [S2_sz0(v[:, k], basis) for k in range(6)]
    sps = [float(v[:, k]**2@space) for k in range(6)]
    miss = np.abs(ex[1:]/(J*hx[1:]) - 1)
    print(f"t/D={t:.2f}: J_fit={J:.6f}  max miss (19 levels, with multiplicity)={miss.max():.3e}")
    print("    lowest 6:  E-E0 =", np.round(ex[:6], 6).tolist())
    print("               S    =", np.round(Ss, 3).tolist())
    print("            <space> =", np.round(sps, 3).tolist())
