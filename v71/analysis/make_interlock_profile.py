#!/usr/bin/env python3
"""make_interlock_profile.py — solve the TWO-LOW/ONE-HIGH flavored baryon
(omega = 1.38, 1.38, 1.42) for the phase-interlock molecule experiment.

The interlock needs the two anti-phased (repulsive) channels to be the
SHORT-range ones (low omega -> larger mu = sqrt(m^2-w^2)) with combined
amplitude beating the single aligned long-range channel at contact:
rep/attr(D) = [2 f_low^2 e^{-mu_low D}] / [f_high^2 e^{-mu_high D}].

Output: v71/results/interlock_profile.txt (r f0 f1 f2; f0=f1 low, f2 high).
"""
import importlib.util
import os
import numpy as np

ROOT = "/home/d/code/scp"
spec = importlib.util.spec_from_file_location(
    "fq", os.path.join(ROOT, "v71/analysis/flavored_qball.py"))
fq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fq)

r, H, N = fq.r, fq.H, fq.N

# symmetric start
src = np.loadtxt(os.path.join(ROOT, "v69/theory/gprofile_w142_g005.txt"))
f0 = np.interp(r, src[:, 0], src[:, 1], right=0.0)
F = np.vstack([f0, f0, f0])
w = np.array([1.42, 1.42, 1.42])
F, ok, it, rn = fq.solve(F, w ** 2)
print(f"symmetric: ok={ok} resid={rn:.1e}")
assert ok

# continuation: lower components 0 AND 1 together
for dw in np.arange(0.005, 0.040 + 1e-9, 0.005):
    w = np.array([1.42 - dw, 1.42 - dw, 1.42])
    F, ok, it, rn = fq.solve(F, w ** 2)
    print(f"dw={dw:.3f}: ok={ok} iters={it} f=({F[0,0]:.4f},{F[1,0]:.4f},{F[2,0]:.4f})")
    if not ok:
        raise SystemExit("continuation failed")

Q, E = fq.observables(F, w)
mu_low, mu_high = np.sqrt(2.25 - w[0] ** 2), np.sqrt(2.25 - w[2] ** 2)
ratio0 = 2 * F[0, 0] ** 2 / F[2, 0] ** 2
Dstar = np.log(ratio0) / (mu_low - mu_high)
print(f"\nw = ({w[0]:.4f}, {w[1]:.4f}, {w[2]:.4f})")
print(f"Q_a = ({Q[0]:.2f}, {Q[1]:.2f}, {Q[2]:.2f})  E = {E:.2f}")
print(f"mu_low = {mu_low:.4f}  mu_high = {mu_high:.4f}")
print(f"contact rep/attr = {ratio0:.3f}  -> predicted interlock D* ~ {Dstar:.1f}")

out = os.path.join(ROOT, "v71/results/interlock_profile.txt")
np.savetxt(out, np.column_stack([r, F[0], F[1], F[2]]),
           header=f"r f0 f1 f2   (w = {w[0]:.4f} {w[1]:.4f} {w[2]:.4f}; "
                  f"interlock: anti on comps 0,1 / aligned on comp 2)")
print(f"wrote {out}")
