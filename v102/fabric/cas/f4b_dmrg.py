#!/usr/bin/env python3
"""F4(b) — the hole in a long effective open chain (reading V), DMRG.
Chain: [cap 1 end] + [cap 2]*Lh + [cap 1 hole] + [cap 2]*Lh + [cap 1 end].
R6 end cells screen the Haldane edge states, so the only free spin is the
hole's.  usage: f4b_dmrg.py L_bulk twoSz chi [hole|nohole]
Prints the energy and the S^z profile."""
import sys, time
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../phase1/cas'); sys.path.insert(0, here + '/../../lock/cas')
import eff
from dmrg import DMRG, spin_ops
from mpo3 import mpo_chain

t0 = time.time()
Lb, tw, chi = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
mode = sys.argv[4] if len(sys.argv) > 4 else 'hole'
h2, r3, co = eff.terms_for(0.1); J = co['J']
Lh = Lb//2
caps = [1] + [2]*Lh + ([1] if mode == 'hole' else [2]) + [2]*Lh + [1]
L = len(caps); dims = [c + 1 for c in caps]
H2 = [h2[(caps[i], caps[i + 1])]/J for i in range(L - 1)] + [None]
R3 = [r3[(caps[i], caps[i + 1], caps[i + 2])]/J for i in range(L - 2)] + [None, None]
qops = [spin_ops(c/2)[0] for c in caps]
W, const = mpo_chain(dims, H2, R3, qops, 1.0, tw/2)
d = DMRG(W, dims, chi=chi, seed=17 + tw)
E = d.run(sweeps=16, etol=1e-11) + const
prof = [d.expect1(qops[i], i) for i in range(L)]
print(f"{mode} Lb={Lb} twoSz={tw} chi={chi} E/J={E:.10f} sumSz={sum(prof):.6f} [{time.time()-t0:.0f}s]")
print("profile " + " ".join(f"{x:+.5f}" for x in prof), flush=True)
