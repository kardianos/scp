#!/usr/bin/env python3
# =========================================================================
# formfind_v2s.py — W1 re-fight on the ANNEALED substrate (task #4 arm 2).
#
# F3: addresses are substrate-dependent. On V2s (geom_dtar=1.52, link
# ceiling ~1.70) the strut length pi/omega at desert loads (x~0.6 ->
# d=1.86) exceeds the link ceiling: NO strut structure can ride the
# desert. Cable rings CAN — the winding integer m sets the address:
#     omega = 2*pi*m / L_ring,   x = (w2/omega - 1)/q_detune
# ring12 @ d~1.52: m=5 -> x~0.57 (mid-desert, at the M-B4 minimum);
#                  m=6 -> x=0.3356 (the strut/rung address — control);
# ring8:           m=3 -> x~0.73 (upper desert).
#
# Emits (candidates/):
#   v2s_ring12_m5.net       comp12-class desert rider
#   v2s_ring12_m5_ctrl.net  same cells/loads, phases mean-tuned BFS
#                           ignoring lengths (W1 ctrl convention)
#   v2s_ring12_m6.net       same cells, rung winding (address A/B control)
#   v2s_ring8_m3.net        upper-desert 8-ring
#   v2s_free1.net           lone voice at x=0.60 (load-line baseline)
#
# Substrate: prestress/foam/anneal_snaps/cells_000000.tsv — the annealed
# dump (NC=5039), byte-identical geometry to any run with
# seed=20260727 dmin=1.25 geom_relax=800 geom_relax_k=0.3 geom_dtar=1.52
# geom_lmax=1.70 L=24. Picks are exact coordinates -> worst_pick=0.
# Deterministic: no RNG anywhere (exhaustive scored search).
# =========================================================================
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TSV = os.path.join(HERE, "foam", "anneal_snaps", "cells_000000.tsv")
OUT = os.path.join(HERE, "candidates")
TWO_PI = 2.0 * math.pi
W2, Q, CAP, L_BOX = 2.9, 1.2, 2.5, 24.0
MARGIN = 2.8            # ring cells stay this far from box faces (W1 rule)
LINK_F = 1.15           # kernel link rule: d < 1.15*(r_i+r_j)

rows = np.loadtxt(TSV, skiprows=1, usecols=(0, 1, 2, 3, 4))
ids, P, R = rows[:, 0].astype(int), rows[:, 1:4], rows[:, 4]
N = len(ids)

GEOM_LMAX = 1.70        # V2s channel ceiling (cellfab.c:968 — cut = min(1.15*(ri+rj), geom_lmax))
D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
CUT = np.minimum(LINK_F * (R[:, None] + R[None, :]), GEOM_LMAX)
LINK = (D < CUT) & ~np.eye(N, dtype=bool)
deg = LINK.sum(1)
ld = D[LINK][::2] if False else D[np.triu(LINK)]
print(f"substrate: NC={N} degree={deg.mean():.2f} "
      f"dbar={ld.mean():.4f} sigma_d={100*ld.std()/ld.mean():.2f}%")

interior = np.all((P > MARGIN) & (P < L_BOX - MARGIN), axis=1)

# 13 deterministic plane normals: axes, face diagonals, body diagonals
NORMALS = []
for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1),
          (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, -1, 0), (1, 0, -1), (0, 1, -1),
          (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]:
    v = np.array(v, float)
    NORMALS.append(v / np.linalg.norm(v))


def find_ring(nring):
    """Best chord-free cycle of nring linked cells near an ideal circle."""
    dbar = ld.mean()
    R_ring = nring * dbar / TWO_PI
    best = None
    for ci in np.where(interior)[0]:
        c = P[ci]
        if np.any((c < MARGIN + R_ring - 0.8) | (c > L_BOX - MARGIN - R_ring + 0.8)):
            pass  # ideal circle may poke near faces; per-cell margin re-checked below
        for nrm in NORMALS:
            u = np.cross(nrm, [1.0, 0.0, 0.0])
            if np.linalg.norm(u) < 1e-6:
                u = np.cross(nrm, [0.0, 1.0, 0.0])
            u /= np.linalg.norm(u)
            v = np.cross(nrm, u)
            picks = []
            ok = True
            for k in range(nring):
                ph = TWO_PI * k / nring
                tgt = c + R_ring * (math.cos(ph) * u + math.sin(ph) * v)
                d2 = np.einsum("ij,ij->i", P - tgt, P - tgt)
                for pi_ in picks:
                    d2[pi_] = 1e18
                j = int(np.argmin(d2))
                if d2[j] > 1.0**2 or not interior[j]:
                    ok = False
                    break
                picks.append(j)
            if not ok or len(set(picks)) != nring:
                continue
            ok = all(LINK[picks[k], picks[(k + 1) % nring]] for k in range(nring))
            if not ok:
                continue
            chords = sum(1 for a in range(nring) for b in range(a + 2, nring)
                         if (a, b) != (0, nring - 1) and LINK[picks[a], picks[b]])
            if chords:
                continue
            dl = [D[picks[k], picks[(k + 1) % nring]] for k in range(nring)]
            spread = float(np.std(dl))
            if best is None or spread < best[0]:
                best = (spread, picks, dl)
    return best


def x_of_omega(w):
    return (W2 / w - 1.0) / Q


def emit_ring(fn, picks, dl, m, ctrl=False):
    nring = len(picks)
    L_ring = sum(dl)
    w = TWO_PI * m / L_ring
    x = x_of_omega(w)
    th = [0.0]
    step = (L_ring / nring) if ctrl else None
    for k in range(nring - 1):
        d = step if ctrl else dl[k]
        th.append((th[-1] - w * d) % TWO_PI)
    closure = (th[-1] - w * (step if ctrl else dl[-1]) - th[0]) % TWO_PI
    closure = min(closure, TWO_PI - closure)
    gerr = (np.array(dl) - L_ring / nring) * w if ctrl else np.zeros(nring)
    path = os.path.join(OUT, fn)
    with open(path, "w") as f:
        f.write(f"# formfind_v2s — annealed-substrate W1 re-fight candidate\n")
        f.write(f"# nring={nring} m={m} L_ring={L_ring:.4f} omega={w:.6f} "
                f"x={x:.6f} Em/voice={x*CAP:.4f} mass={nring*x*CAP:.4f}\n")
        f.write(f"# edge spread sigma={np.std(dl):.4f} "
                f"lengths=[{min(dl):.4f},{max(dl):.4f}] "
                f"closure_defect={closure:.2e} ctrl={int(ctrl)} "
                f"ctrl_gate_rms={float(np.sqrt((gerr**2).mean())):.4f}\n")
        for k, pi_ in enumerate(picks):
            f.write(f"V {P[pi_,0]:.4f} {P[pi_,1]:.4f} {P[pi_,2]:.4f} "
                    f"{x:.6f} {th[k]:.6f}\n")
        for k in range(nring):
            f.write(f"E {k} {(k + 1) % nring}\n")
    print(f"{fn}: nring={nring} m={m} L={L_ring:.3f} omega={w:.4f} x={x:.4f} "
          f"spread={np.std(dl):.4f} closure={closure:.1e}")


r12 = find_ring(12)
r8 = find_ring(8)
if r12 is None or r8 is None:
    raise SystemExit("no valid ring found")

emit_ring("v2s_ring12_m5.net", r12[1], r12[2], 5)
emit_ring("v2s_ring12_m5_ctrl.net", r12[1], r12[2], 5, ctrl=True)
emit_ring("v2s_ring12_m6.net", r12[1], r12[2], 6)
emit_ring("v2s_ring8_m3.net", r8[1], r8[2], 3)

# free1: interior cell nearest box center (each net runs in its own box;
# centrality keeps the footprint clear of the sink rim on all sides)
ctr = np.array([L_BOX / 2] * 3)
far = int(np.argmin(np.linalg.norm(P - ctr, axis=1)))
x_free = 0.60
with open(os.path.join(OUT, "v2s_free1.net"), "w") as f:
    f.write("# formfind_v2s — lone desert voice (load-line baseline)\n")
    f.write(f"# x={x_free} Em={x_free*CAP:.4f} cell={ids[far]} "
            f"dist_to_center={np.linalg.norm(P[far]-ctr):.2f}\n")
    f.write(f"V {P[far,0]:.4f} {P[far,1]:.4f} {P[far,2]:.4f} {x_free:.6f} 0.000000\n")
print(f"v2s_free1.net: cell {ids[far]} x={x_free} "
      f"center_dist={np.linalg.norm(P[far]-ctr):.2f} pos={P[far].round(2)}")
