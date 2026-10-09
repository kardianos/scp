#!/usr/bin/env python3
"""DMRG on leading-order chains (lo.py), energies in units of J0 = 2 g^2/Delta.
usage: lo_dmrg.py M mode L_bulk target t_over_J chi
  M: 2 (3+2) or 4 (W);  mode: nohole | hole | mob
  target: comma list of conserved charges, 2S^z[,2I^z]
Chain: [end] + bulk + [end]; hole: a static hole in the middle;
mob: every bulk site is 'mob' and one hole is pinned in number by a penalty.
Prints E and the profiles of S^z (, I^z) (, hole density)."""
import sys, time
import numpy as np
here = __file__.rsplit('/', 1)[0]
sys.path.insert(0, here); sys.path.insert(0, here + '/../../phase1/cas')
import graphm as GM, lo
from dmrg import DMRG, build_mpo
from dmrg_sym import DMRGSym

t0 = time.time()
M, mode, Lb = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
target = tuple(int(x) for x in sys.argv[4].split(','))
tJ, chi = float(sys.argv[5]), int(sys.argv[6])
P = GM.P_W if M == 4 else GM.P_2
g = 0.1*GM.DELTA; J0 = 2*g*g/GM.DELTA
Lh = Lb//2
if mode == 'nohole':
    kinds = ['end'] + ['full']*Lb + ['end']
elif mode == 'hole':
    kinds = ['end'] + ['full']*Lh + ['hole'] + ['full']*Lh + ['end']
else:
    kinds = ['end'] + ['mob']*(Lb + 1) + ['end']
locs = [lo.Local(M, k) for k in kinds]
L = len(locs)
terms = [(i, ops, c/J0) for (i, ops, c) in lo.dmrg_terms(locs, P, g, t=tJ*J0)]
pens = [([Lc.Sz for Lc in locs], 1.0, target[0]/2)]
if M == 4: pens.append(([Lc.Iz for Lc in locs], 1.0, target[1]/2))
if mode == 'mob': pens.append(([Lc.nh if Lc.nh is not None else np.zeros((Lc.d, Lc.d)) for Lc in locs], 1.0, 1.0))
W, const = build_mpo([Lc.d for Lc in locs], [None]*L, terms, pens)
# product initial state in the target sector
init = []
sgn = 1
for i, Lc in enumerate(locs):
    labs = Lc.lab
    if Lc.kind == 'end' or (Lc.kind == 'hole'):
        want = np.array([sgn]*labs.shape[1]); sgn = -sgn
    else:
        want = np.array([2*sgn]*labs.shape[1]); sgn = -sgn
    k = int(np.argmin(np.abs(labs - want[None, :]).sum(1) + (0 if Lc.kind != 'mob' else 10*(Lc.N != 2))))
    init.append(k)
if mode == 'mob':
    mid = L//2
    Lc = locs[mid]
    cand = [k for k in range(Lc.d) if Lc.N[k] == 1]
    init[mid] = cand[0]
# fix the sector of the product state by adjusting the middle sites greedily
def tot(init):
    return np.sum([locs[i].lab[init[i]] for i in range(L)], 0)
for _ in range(4*L):
    diff = np.array(target) - tot(init)
    if not np.any(diff): break
    best = None
    for i in range(L):
        Lc = locs[i]
        for k in range(Lc.d):
            if mode == 'mob' and Lc.kind == 'mob' and Lc.N[k] != Lc.N[init[i]]: continue
            nd = np.abs(diff - (Lc.lab[k] - Lc.lab[init[i]])).sum()
            if best is None or nd < best[0]: best = (nd, i, k)
    init[best[1]] = best[2]
d = (DMRGSym if (len(sys.argv) > 7 and sys.argv[7] == 'sym') else DMRG)(W, [Lc.d for Lc in locs], chi=chi, seed=5, init=init)
E = d.run(sweeps=14, etol=1e-9) + const
prof = [d.expect1(Lc.Sz, i) for i, Lc in enumerate(locs)]
line = f"M={M} mode={mode} Lb={Lb} target={target} t/J0={tJ} chi={chi} {sys.argv[7] if len(sys.argv) > 7 else 'plain'} E/J0={E:.10f} sumSz={sum(prof):.5f} trunc={getattr(d, 'last_trunc', float('nan')):.1e}"
out = ["Sz " + " ".join(f"{x:+.5f}" for x in prof)]
if M == 4:
    pi = [d.expect1(Lc.Iz, i) for i, Lc in enumerate(locs)]
    line += f" sumIz={sum(pi):.5f}"; out.append("Iz " + " ".join(f"{x:+.5f}" for x in pi))
if mode == 'mob':
    ph = [d.expect1(Lc.nh, i) if Lc.nh is not None else 0.0 for i, Lc in enumerate(locs)]
    line += f" holes={sum(ph):.5f}"; out.append("nh " + " ".join(f"{x:+.5f}" for x in ph))
    # hole-centred spin profile C(d) = sum_i <n_h(i) S^z(i+d)> + <S^z(i-d) n_h(i)>  (d >= 1), C(0) = sum_i <n_h S^z>_i
    nh_ops = [Lc.nh if Lc.nh is not None else np.zeros((Lc.d, Lc.d)) for Lc in locs]
    sz_ops = [Lc.Sz for Lc in locs]
    C = np.zeros(L)
    for i in range(1, L - 1):
        C[0] += d.expect1(nh_ops[i]@sz_ops[i], i)
        js = list(range(i + 1, L))
        if js:
            row = d.corr_row(nh_ops[i], i, sz_ops, js)
            for j, v in row.items(): C[j - i] += v
        js2 = list(range(i + 1, L))
        for k in range(0, i):
            pass
    for k in range(0, L - 1):
        js = [j for j in range(k + 1, L - 1)]
        if js:
            row = d.corr_row(sz_ops[k], k, nh_ops, js)
            for j, v in row.items(): C[j - k] += v
    out.append("C " + " ".join(f"{x:+.5f}" for x in C))
print(line + f" [{time.time()-t0:.0f}s]")
print("\n".join(out), flush=True)
