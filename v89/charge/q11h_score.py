#!/usr/bin/env python3
"""Q11h scorer — holonomy current vs harmonic (Tenney height) and size.

Pre-registered (before scoring, recorded in the session + §7.10 draft):
  P1 frustration: fifth TRIANGLES (perimeter-quantized defect, cannot
     close in-band) circulate persistently, site-independently;
     closable cycles (oct_tri defect 0.0005, fifth_hex −0.089) can
     truly lock (fix_locked=1, fix_nu ~ 0).
  P2 Tenney fragility: octave tongue (γ/2) holds vs fifth (γ/6):
     flip rate collapses vs the fifth's ~867/600 t.u.
  P3 C-symmetry: circulation antisymmetry around every cycle
     (Kirchhoff closure of fix_slips).
Anchor row: the P17 fifth triangle (q11_pinned_v2.log).
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = {
    "fifth_tri(P17)": "runs/q11_pinned_v2.log",
    "fifth_tri2":     "runs/q11h_fifth_tri2.log",
    "oct_tri":        "runs/q11h_oct_tri.log",
    "fifth_hex":      "runs/q11h_fifth_hex.log",
}
TEN = {"3:2": 6, "2:3": 6, "2:1": 2, "1:2": 2, "1:1": 1}

print(f"{'net':>15} {'e':>2} {'fixpq':>5} {'Ten':>3} {'fix_nu':>10} "
      f"{'fix_dwb2pi':>11} {'ratio':>6} {'fix_slips':>9} {'flk':>3} "
      f"{'flips':>6} {'livepq':>6}")
for name, rel in RUNS.items():
    path = os.path.join(HERE, rel)
    if not os.path.exists(path):
        print(f"{name:>15}  -- missing --")
        continue
    ksum = 0.0
    for line in open(path):
        if not line.startswith("# RESULT slip"):
            continue
        g = lambda k, s=line: re.search(k + r"=(\S+)", s).group(1)
        fixpq = g("fix_pq")
        nu, db = float(g("fix_nu")), float(g("fix_dwb2pi"))
        sl = float(g("fix_slips"))
        ksum += sl
        ratio = nu / db if abs(db) > 1e-9 else float("nan")
        print(f"{name:>15} {g('e'):>2} {fixpq:>5} {TEN.get(fixpq,0):>3} "
              f"{nu:>10.6f} {db:>11.6f} {ratio:>6.2f} {sl:>9.2f} "
              f"{g('fix_locked'):>3} {g('flips'):>6} {g('pq'):>6}")
    print(f"{'':>15}  cycle fix_slips sum (Kirchhoff) = {ksum:+.2f}")
