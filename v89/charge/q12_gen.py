#!/usr/bin/env python3
"""Q12 arm 1 — the kink launcher: a Z3 branch wall on a fifth chain.

A chain (open path, no closure constraint) of alternating light/heavy
voices (x = 0.195 / 0.709, the tri2 fifth pair) with every 3:2 gate
chained exactly open. The KINK: one mid-chain edge carries an extra
+2pi/3 — a real Z3 domain wall (on 3:2 edges the chart branches are
Z3; unlike a 2pi offset this is NOT gauge). Question: does the wall
propagate under the gate dynamics (mobile dressed kink — the Q12
preparation) or pin like cycle strain?

Emits: nets/q12_chain.net (kinked), nets/q12_chain_ctrl.net (no kink).
Deterministic; picks cells from the production foam TSV.
"""
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FOAM = os.path.join(HERE, "..", "prestress", "foam", "foam_s20260727.tsv")
OUT = os.path.join(HERE, "nets")
os.makedirs(OUT, exist_ok=True)
TWO_PI = 2.0 * math.pi
W2, Q, CAP, L_BOX = 2.9, 1.2, 2.5, 24.0
XA, XB = 0.195035, 0.709220        # the tri2 fifth pair (3:2 exactly)
NCHAIN = 14
MARGIN = 2.8

rows = np.loadtxt(FOAM, skiprows=1, usecols=(0, 1, 2, 3, 4))
P, R = rows[:, 1:4], rows[:, 4]
N = len(P)
D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
LINK = (D < 1.15 * (R[:, None] + R[None, :])) & ~np.eye(N, dtype=bool)
interior = np.all((P > MARGIN) & (P < L_BOX - MARGIN), axis=1)


def walk_chain(n):
    """Straightest interior linked path of n cells: greedy, scored by
    total turning; deterministic (start cells in id order)."""
    best = None
    for start in np.where(interior)[0][:400]:
        for d0 in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            path = [start]
            dirv = np.array(d0, float)
            turn = 0.0
            ok = True
            for _ in range(n - 1):
                cands = [j for j in np.where(LINK[path[-1]])[0]
                         if interior[j] and j not in path
                         and 1.30 <= D[path[-1], j] <= 1.70]
                if not cands:
                    ok = False
                    break
                scores = []
                for j in cands:
                    v = P[j] - P[path[-1]]
                    v = v / np.linalg.norm(v)
                    scores.append((1.0 - float(v @ dirv), j, v))
                scores.sort(key=lambda s: (s[0], s[1]))
                c, j, v = scores[0]
                turn += c
                path.append(int(j))
                dirv = v
            if ok and (best is None or turn < best[0]):
                best = (turn, path)
    return best[1]


chain = walk_chain(NCHAIN)
dl = [float(np.linalg.norm(P[chain[k + 1]] - P[chain[k]])) for k in range(NCHAIN - 1)]
wA = W2 / (1 + Q * XA)
wB = W2 / (1 + Q * XB)
assert abs(wA / wB - 1.5) < 1e-3, (wA, wB, wA / wB)

for fn, kink in (("q12_chain.net", True), ("q12_chain_ctrl.net", False)):
    xs = [XA if k % 2 == 0 else XB for k in range(NCHAIN)]
    ws = [wA if k % 2 == 0 else wB for k in range(NCHAIN)]
    th = [0.0]
    mid = NCHAIN // 2 - 1
    for k in range(NCHAIN - 1):
        # forward gate of edge k (i=k, j=k+1) open: th_j = th_i - w_i*d
        t = th[-1] - ws[k] * dl[k]
        if kink and k == mid:
            t += TWO_PI / 3.0
        th.append(t % TWO_PI)
    with open(os.path.join(OUT, fn), "w") as f:
        f.write(f"# Q12 kink launcher: fifth chain n={NCHAIN} "
                f"wA={wA:.6f} wB={wB:.6f} ratio={wA/wB:.4f}\n")
        f.write(f"# kink={int(kink)} at edge {mid} (+2pi/3 Z3 branch wall); "
                f"lengths=[{min(dl):.4f},{max(dl):.4f}]\n")
        for k, c in enumerate(chain):
            f.write(f"V {P[c,0]:.4f} {P[c,1]:.4f} {P[c,2]:.4f} "
                    f"{xs[k]:.6f} {th[k]:.6f}\n")
        for k in range(NCHAIN - 1):
            # 3:2 edges: p:q = w_i:w_j reduced — light->heavy is 3:2,
            # heavy->light is 2:3
            p, q = (3, 2) if k % 2 == 0 else (2, 3)
            f.write(f"E {k} {k + 1} {p} {q}\n")
    print(f"{fn}: chain={chain[:4]}...{chain[-2:]} kink_edge={mid if kink else '-'} "
          f"span={np.linalg.norm(P[chain[-1]]-P[chain[0]]):.2f}")
