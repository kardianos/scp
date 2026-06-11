#!/usr/bin/env python3
"""make_w146_profile.py — generate the on-branch gauged profile at
omega=1.46, g=0.05 by omega-continuation from the existing w142 profile,
reusing the v69 gauged shooter's solver. Output:
v70/results/gprofile_w146_g005.txt (same 4-column contract as v69 §5.1).
"""
import importlib.util
import os
import numpy as np

ROOT = "/home/d/code/scp"
spec = importlib.util.spec_from_file_location(
    "gsf", os.path.join(ROOT, "v69/theory/gauged_shooter_fast.py"))
gsf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gsf)

H, N, RFAC = gsf.H, gsf.N, gsf.RFAC
r = np.arange(N) * H

# seed from the w142 g=0.05 profile (r f Er weff)
src = np.loadtxt(os.path.join(ROOT, "v69/theory/gprofile_w142_g005.txt"))
f_seed = np.interp(r, src[:, 0], src[:, 1], right=0.0)
chi_seed = np.interp(r, src[:, 0], 1.42 - src[:, 3], right=0.0)

g = 0.05
w, f, chi = 1.42, f_seed.copy(), chi_seed.copy()
f, chi, ok, its, rn = gsf.solve(w, g, f, chi)
print(f"re-solve at w=1.42: ok={ok} iters={its} resid={rn:.2e}")
assert ok

for wt in np.arange(1.4225, 1.46 + 1e-9, 0.0025):
    f, chi, ok, its, rn = gsf.solve(wt, g, f, chi)
    print(f"w={wt:.4f}: ok={ok} iters={its} resid={rn:.2e} f0={f[0]:.5f}")
    if not ok:
        raise SystemExit("continuation failed")
    w = wt

o = gsf.observables(f, chi, w, g)
wt_arr = w - chi
integ = np.concatenate(([0.0], np.cumsum(
    0.5 * H * (3 * wt_arr * f * f * r * r)[:-1]
    + 0.5 * H * (3 * wt_arr * f * f * r * r)[1:])))
Er = np.zeros(N)
Er[1:] = g * integ[1:] / (r[1:] ** 2)

out = os.path.join(ROOT, "v70/results/gprofile_w146_g005.txt")
step = int(round(0.02 / H))
nmax = int(round(60.0 / 0.02))
with open(out, "w") as fp:
    fp.write("# v70 gauged_shooter profile (omega-continued from w142): "
             "omega=%.6f g=%.6f m2=%.6f mu=%.6f kappa=%.6f\n"
             % (w, g, gsf.M2, gsf.MU, gsf.KAP))
    fp.write("# f0=%.6f Q=%.2f E_matter=%.2f E_field=%.3f E_total=%.2f "
             "chi0=%.6f a0_GD(0)=%.6f weff0=%.6f\n"
             % (o["f0"], o["Q"], o["Em"], o["Ef"], o["Et"], o["chi0"],
                o["a00"], o["weff0"]))
    fp.write("# r f Er weff\n")
    for k in range(nmax + 1):
        i = k * step
        fp.write("%.6f %.9f %.9e %.9f\n" % (r[i], max(f[i], 0.0), Er[i], wt_arr[i]))
print(f"wrote {out}")
print(f"w=1.46 ball: Q={o['Q']:.2f} E_total={o['Et']:.2f} weff0={o['weff0']:.5f} "
      f"r_half={o['rhalf']:.3f}")
