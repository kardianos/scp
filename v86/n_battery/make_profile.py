#!/usr/bin/env python3
"""Emit a gauged-shooter radial profile in the v69 gprofile format at a chosen
(omega, g), for N7 seeding. Same columns and same temporal-gauge convention as
v69/theory/gauged_shooter.py writes, so downstream tooling is unchanged."""
import os, sys
import numpy as np
ROOT = "/home/d/code/scp"
sys.path.insert(0, os.path.join(ROOT, "v69/theory"))
import gauged_shooter as G          # production grid H=0.004 RMAX=150

w_t = float(sys.argv[1]); g_t = float(sys.argv[2]); out = sys.argv[3]
f0g = G.load_v66_profile(os.path.join(ROOT, "v66/results/profile_omega1.4500.txt"))
f, chi, ok, _, rn = G.solve(1.45, 0.0, f0g, np.zeros(G.N))
assert ok, "seed failed"
# march omega at g=0
w = 1.45; step = -0.004 if w_t < w else 0.004
while abs(w - w_t) > 1e-9:
    wn = w + step
    if (step < 0 and wn < w_t) or (step > 0 and wn > w_t): wn = w_t
    fn, cn, okk, _, _ = G.solve(wn, 0.0, f, chi)
    if not okk:
        step *= 0.5; assert abs(step) > 1e-6, "omega continuation stalled"; continue
    f, chi, w = fn, cn, wn
# continue in g
gc, chig = 0.0, 0.0
while gc < g_t - 1e-12:
    dg = g_t - gc
    while True:
        gg = min(gc + dg, g_t)
        cg = chi * (gg / chig) ** 2 if chig > 0 else chi
        fn, cn, okk, _, rn = G.solve(w_t, gg, f, cg)
        if okk: f, chi, gc, chig = fn, cn, gg, gg; break
        dg *= 0.5
        assert dg > 1e-7, "g continuation stalled"
o = G.observables(f, chi, w_t, g_t)
r, H = G.r, G.H
wt = w_t - chi
cpe = np.gradient(np.append(chi, chi[-1]*G.RFAC), H, edge_order=2)[:len(r)]
Er = -cpe / g_t
print("w=%.4f g=%.3f: f0=%.6f Q=%.2f E=%.2f r_half=%.3f weff0=%.5f resid=%.1e"
      % (w_t, g_t, o["f0"], o["Q"], o["Et"], o["rhalf"], o["weff0"], rn))
step_i = int(round(0.02 / H)); nmax = int(round(60.0 / 0.02))
with open(out, "w") as fp:
    fp.write("# v86 N7 profile (gauged_shooter): omega=%.6f g=%.6f m2=%.6f mu=%.6f kappa=%.6f\n"
             % (w_t, g_t, G.M2, G.MU, G.KAP))
    fp.write("# f0=%.6f Q=%.2f E_matter=%.2f E_field=%.3f E_total=%.2f chi0=%.6f weff0=%.6f r_half=%.4f\n"
             % (o["f0"], o["Q"], o["Em"], o["Ef"], o["Et"], o["chi0"], o["weff0"], o["rhalf"]))
    fp.write("# r f Er weff\n")
    for k in range(nmax + 1):
        i = k * step_i
        fp.write("%.6f %.9f %.9e %.9f\n" % (r[i], max(f[i], 0.0), Er[i], wt[i]))
print("wrote", out)
