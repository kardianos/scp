#!/usr/bin/env python3
"""L4(a) — one DMRG run of the effective chain (../PREREG.md §3 L4).
Usage: l4_mass.py {slaw|heis} L_bulk Sz chi
slaw: h2 + r3 from the S-law clusters at g/Delta = 0.1 (clusters.py),
      in units of the bulk J_h2.  heis: spin-1 Heisenberg, J = 1.
Chain: [cap 1] + [cap 2]*L_bulk + [cap 1]  (R6 end cells), so the
spin-1/2 ends screen the edge states (P2's chain).
Prints one line: mode L_bulk Sz chi E J."""
import sys, time
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../phase1/cas')
from dmrg import DMRG, spin_ops
from clusters import effective, SdotS
from mpo3 import mpo_chain
from slaw import DELTA

mode, Lb, Sz, chi = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
gr = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
t0 = time.time()
if mode == 'slaw':
    h2, r3, co, info = effective(gr*DELTA); J = co['J']
else:
    h2 = {(2, 2): SdotS(1, 1), (1, 2): SdotS(0.5, 1), (2, 1): SdotS(1, 0.5)}
    r3 = {}; J = 1.0
caps = [1] + [2]*Lb + [1]; L = len(caps); dims = [c + 1 for c in caps]
H2 = [h2[(caps[i], caps[i + 1])]/J for i in range(L - 1)] + [None]
R3 = [(r3[(caps[i], caps[i + 1], caps[i + 2])]/J if (caps[i], caps[i + 1], caps[i + 2]) in r3 else None)
      for i in range(L - 2)] + [None, None]
qops = [spin_ops(c/2)[0] for c in caps]
W, const = mpo_chain(dims, H2, R3, qops, 1.0, Sz)          # compact exact MPO (mpo3.py)
d = DMRG(W, dims, chi=chi, seed=7 + int(Sz))
E = d.run(sweeps=14, etol=1e-11) + const
print(f"{mode} {Lb} {Sz:.0f} {chi} {E:.10f} {J:.8f} {gr} [{time.time()-t0:.0f}s]", flush=True)
