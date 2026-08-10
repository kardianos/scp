#!/usr/bin/env python3
"""fig_fold: (a) the dressing fold (theory); (b) M3b measured door_in vs field energy."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.6))

# (a) u = beta (s + q u)^2 : sweep the lower branch to the fold
q, s = 1.2, 1.0
u = np.linspace(0, 1.0 / (2 * q * 0.12) - s / q, 400)  # generous range
betas = u / (s + q * u) ** 2
uc_star = []
b_grid = np.linspace(0.005, 1 / (4 * q * s), 200)
low, high = [], []
for b in b_grid:
    r = np.roots([b * q * q, 2 * b * q * s - 1, b * s * s])
    r = np.sort(r[np.abs(r.imag) < 1e-9].real)
    if len(r) == 0:
        low.append(np.nan); high.append(np.nan)
    else:
        low.append(r[0]); high.append(r[1] if len(r) > 1 else np.nan)
ax1.plot(b_grid, low, color=INK, lw=1.6)
ax1.plot(b_grid, high, color=GRAY, lw=1.2, linestyle=(0, (2, 2)))
bc = 1 / (4 * q * s)
uc = 1 / (2 * bc * q * q) - s / q
ax1.plot([bc], [np.interp(bc, b_grid, low)], marker="o", ms=5, color=ORANGE)
ax1.annotate("fold: no bounded slaved\nsolution beyond here",
             xy=(bc, np.interp(bc, b_grid, low)), xytext=(0.075, 1.15),
             fontsize=7.5,
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8))
ax1.text(0.02, 0.17, "stable (reversible)\ndressing branch", fontsize=7.5)
ax1.text(0.10, 0.68, "unstable branch", color=GRAY, fontsize=7.5)
ax1.set_xlabel(r"drive  $\beta = g^2 E_e/(\Delta_0^2\, \mathrm{cap})$")
ax1.set_ylabel(r"dressing  $u=|\chi|^2/\mathrm{cap}$")
ax1.set_xlim(0, 0.24); ax1.set_ylim(0, 1.5)
ax1.set_title("(a) predicted dressing fold", fontsize=8)

# (b) measured M3b (Go kernel, door_g = 0.5, blob+bath, T = 500)
na = np.array([0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0])
Ee = na ** 2
din = np.array([0, 0, 2.11, 8.46, 19.93, 50.84, 153.74])
dout = np.array([0, 0, 0, 0, 0.86, 5.30, 22.37])
ax2.plot(Ee, din, marker="o", ms=4.5, color=ORANGE, lw=1.3)
ax2.plot(Ee, dout, marker="s", ms=4, color=BLUE, lw=1.3)
ax2.text(3.1, 60, "door_in\n(intake)", color=ORANGE, fontsize=8, ha="left")
ax2.text(3.1, 4.5, "door_out\n(outtake)", color=BLUE, fontsize=8, ha="left")
ax2.axvline(0.36, color=GRAY, lw=0.8, linestyle=(0, (2, 2)))
ax2.annotate("onset = one atom\n($E_e \\approx 0.36 \\approx \\varepsilon$)",
             xy=(0.36, 1.2), xytext=(0.05, 25), fontsize=7,
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.7))
ax2.axvline(1.0, color=GRAY, lw=0.8, linestyle=(0, (4, 2)))
ax2.annotate("knee $\\approx$ fold\n($E_e \\approx 1.0$)",
             xy=(1.0, 0.05), xytext=(1.35, -0.35), fontsize=7,
             arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.7))
ax2.set_ylim(-0.45, 260)
ax2.set_xlabel(r"per-cell field energy  $E_e \approx n_a^2$")
ax2.set_ylabel("door traffic (total)")
ax2.set_yscale("symlog", linthresh=0.5)
ax2.set_xlim(0, 4.3)
ax2.set_title("(b) measured coherent-drive ramp (M3b)", fontsize=8)
for ax in (ax1, ax2):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig("fig_fold.pdf")
print("wrote fig_fold.pdf")
