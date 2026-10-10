#!/usr/bin/env python3
"""K3 check (reported, not gated). usage: k3_ker8.py kernel.out (a 3+2 n=10 run ending with 7 holes)"""
# Is there any R12-allowed state with every bulk site holed (3+2 n=10, 8 holes)? Lowest eigenvalues of
# H_c = sum_l Pi_l in each 8-hole charge sector next to the run's final sector; and the lowest exact level there.
import os, sys, runpy, numpy as np
sys.argv = ['k2d_replay.py', sys.argv[1], 'score', 'proj']
G = runpy.run_path(os.path.dirname(os.path.abspath(__file__)) + '/k2d_replay.py', run_name='k2d')
locs, n, sect, lo, pi_link = G['locs'], G['n'], G['sect'], G['lo'], G['pi_link']
tgt = G['tgt']
for q in (1, -1):
    nt = (tgt[0] - q, tgt[1] + 1)
    H, b, _, _ = sect(nt)
    Hc = sum(lo.sector(locs, [((l, l + 1), pi_link(locs[l], locs[l + 1]))], nt, G['lf'])[0] for l in range(n - 1))
    ev, vec = np.linalg.eigh(Hc.toarray())
    k0 = np.sum(ev < 1e-10)
    print(f"sector {nt}: dim {len(b)}, kernel of H_c: {k0} states; lowest H_c eigenvalues {np.round(ev[:4], 6)}")
    if k0:
        K = vec[:, :k0]; Hk = K.conj().T @ H.toarray() @ K
        print(f"   energies of the allowed states (rate units): {np.round(np.linalg.eigvalsh(Hk), 6)}")
