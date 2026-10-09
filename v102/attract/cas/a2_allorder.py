#!/usr/bin/env python3
"""A2 — 3+2 chains with the all-order (fourth-order-exact) S-law effective
terms h2 + r3 at a given g/Delta, own compact-MPO DMRG (as
../../fabric/cas/f4b_dmrg.py, which is fixed at g/Delta = 0.1).
usage: a2_allorder.py gr L_bulk twoSz chi hole|nohole
Chain: [cap 1 end] + [2]*Lh + [1 hole] + [2]*Lh + [cap 1 end] (hole), or
[cap 1] + [2]*L_bulk + [cap 1] (nohole).  Energies in units of the bulk J."""
import sys, time
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here + '/../../fabric/cas'); sys.path.insert(0, here + '/../../phase1/cas'); sys.path.insert(0, here + '/../../lock/cas')
import eff
from dmrg import DMRG, spin_ops
from mpo3 import mpo_chain
t0 = time.time()
gr, Lb, tw, chi, mode = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
h2, r3, co = eff.terms_for(gr); J = co['J']
Lh = Lb//2
caps = [1] + [2]*Lh + ([1] if mode == 'hole' else [2]) + [2]*Lh + [1] if mode == 'hole' else [1] + [2]*Lb + [1]
L = len(caps); dims = [c + 1 for c in caps]
H2 = [h2[(caps[i], caps[i + 1])]/J for i in range(L - 1)] + [None]
R3 = [r3[(caps[i], caps[i + 1], caps[i + 2])]/J for i in range(L - 2)] + [None, None]
qops = [spin_ops(c/2)[0] for c in caps]
W, const = mpo_chain(dims, H2, R3, qops, 1.0, tw/2)
d = DMRG(W, dims, chi=chi, seed=17 + tw)
E = d.run(sweeps=16, etol=1e-11) + const
prof = [d.expect1(qops[i], i) for i in range(L)]
print(f"gr={gr} mode={mode} Lb={Lb} twoSz={tw} chi={chi} E/J={E:.10f} J={J:.8f} sumSz={sum(prof):.6f} [{time.time()-t0:.0f}s]")
print("profile " + " ".join(f"{x:+.5f}" for x in prof), flush=True)
