#!/usr/bin/env python3
"""Post-registration diagnostic: where does transport of conversion holes
destroy the antiferromagnetic lock?  Ground-state spin and <space> vs t/Delta
on the n = 4 ring (full capped space), for three g/Delta.  Estimate to
test: in a ferromagnetic background each conversion hole hops with
amplitude 2t, so the converted state sits near Delta - 8t (two holes,
band bottom -4t each).  That predicts t_c/Delta ~ 1/8, independent of g
at leading order.  Sensitivity is measured; nothing is chosen."""
import sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from p1_micro import micro, S2_sz0
print("  g/D    t_c/D (first t with S_ground > 0)   <space> just below / at t_c")
for g in (0.05, 0.1, 0.2):
    prev = None; found = None
    for t in np.round(np.arange(0.0, 0.30001, 0.01), 3):
        H, space, basis = micro(4, g, t, cmax=4)
        w, v = np.linalg.eigh(H.toarray())
        S = S2_sz0(v[:, 0], basis); sp0 = float(v[:, 0]**2@space)
        if S > 0.5 and found is None:
            found = (t, sp0, prev); break
        prev = sp0
    if found:
        print(f"  {g:4.2f}   {found[0]:.2f}                                 {found[2]:.3f} / {found[1]:.3f}")
    else:
        print(f"  {g:4.2f}   > 0.30")
