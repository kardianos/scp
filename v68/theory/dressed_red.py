#!/usr/bin/env python3
"""Supplementary scan: map the strong-Gaussian RED regime of the matched-Q shift
(d_omega(n) turned red at n ~ 0.04-0.08 in the main run) and find its maximum
depth before the branch is lost. Also report the implied amplification over
vacuum drive for the e=0.2 and e=0.8 rungs."""
import numpy as np
src = open('/home/d/code/scp/v68/theory/dressed.py').read()
exec(src.split("# ================================================================ CHECK A")[0])

tb0 = make_tables('vac')
w0_300, _ = omega_at_Q(tb0, 300.0)
print(f"vacuum omega(Q=300) = {w0_300:.5f}")
wg = w0_300
best = (None, 0.0)
for n in (0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.14, 0.16):
    tb = make_tables('M2', n_override=n, use_enh=True, H=6)
    w, o = omega_at_Q(tb, 300.0, w_guess=wg)
    if w is None:
        print(f"n={n:5.3f}: BRANCH LOST (omega_min={tb.wmin:.4f})")
        continue
    wg = w
    d = w - w0_300
    if d < best[1]:
        best = (n, d)
    print(f"n={n:5.3f}: d_omega = {d:+.5f}  omega_min = {tb.wmin:.5f}  f0 = {o['f0']:.4f}")
print(f"\ndeepest red: d_omega = {best[1]:+.5f} at n = {best[0]}")
print("targets: -0.025 (e=0.2, n0e=0.0029), -0.120 (e=0.8, n0e=0.0117)")
for target in (-0.025,):
    pass
