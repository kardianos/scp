#!/usr/bin/env python3
"""Q2 (../PREREG.md §3): the W background's correlation from the middle cell
of a long leading-order open chain, by the TeNPy model of
../../perm/cas/tp_chain.py (a thin driver: TeNPy does the work over numpy,
BLAS and LAPACK).  Run with ~/.venvs/scp-tenpy/bin/python -I.
usage: vq2_tp.py M L_bulk chi
Prints E/J0 and <Sz_c Sz_c+r>, <Iz_c Iz_c+r> for r = 0 .. L-1-c."""
import os, sys, time, warnings
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here + '/../../perm/cas')
import graphm as GM, lo, tp_chain as TC
from tenpy.networks.mps import MPS
from tenpy.algorithms import dmrg
warnings.filterwarnings("ignore")

t0 = time.time()
M, Lb, chi = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA
kinds = ['end'] + ['full']*Lb + ['end']
locs = [lo.Local(M, k) for k in kinds]
model = TC.Chain(locs, P, g, 0.0, J0, False)
sites = model.lat.mps_sites()
L = len(locs); sgn = 1; init = []
target = (0, 0) if M == 4 else (0,)
for Lc in locs:
    w = np.array([(1 if Lc.kind == 'end' else 2)*sgn]*Lc.lab.shape[1]); sgn = -sgn
    init.append(int(np.argmin(np.abs(Lc.lab - w[None, :]).sum(1))))
def tot(init): return np.sum([locs[i].lab[init[i]] for i in range(L)], 0)
for _ in range(4*L):
    diff = np.array(target) - tot(init)
    if not np.any(diff): break
    best = None
    for i in range(L):
        for k in range(locs[i].d):
            nd = np.abs(diff - (locs[i].lab[k] - locs[i].lab[init[i]])).sum()
            if best is None or nd < best[0]: best = (nd, i, k)
    init[best[1]] = best[2]
assert not np.any(np.array(target) - tot(init))
psi = MPS.from_product_state(sites, [str(k) for k in init], bc='finite')
opts = {'max_N_sites_per_ring': None, 'max_trunc_err': 1.0, 'mixer': True, 'max_E_err': 1e-10, 'max_S_err': 1e-7,
        'min_sweeps': 6, 'max_sweeps': 20, 'trunc_params': {'chi_max': chi, 'svd_min': 1e-12}}
info = dmrg.run(psi, model, opts)
psi.canonical_form()
c = L//2
js = list(range(c, L))
czz = psi.correlation_function('Sz', 'Sz', sites1=[c], sites2=js)[0]
print(f"M={M} nohole Lb={Lb} chi={chi} E/J0={info['E']:.10f} trunc={float(np.max(info['sweep_statistics']['max_trunc_err'][-1:])):.1e} centre={c} [{time.time()-t0:.0f}s]")
print("SzSz " + " ".join(f"{x:+.10f}" for x in czz))
if M == 4:
    cii = psi.correlation_function('Iz', 'Iz', sites1=[c], sites2=js)[0]
    print("IzIz " + " ".join(f"{x:+.10f}" for x in cii))
