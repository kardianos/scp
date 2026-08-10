#!/usr/bin/env python3
"""fig_doubleslit: two-path phase sum on the roster -> fringes; click histogram rebuilds them."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.4, 3.3), sharex=True,
                               gridspec_kw={"height_ratios": [1, 1]})

y = np.linspace(-6, 6, 800)
d, L, lam = 2.0, 10.0, 1.0          # slit spacing, screen distance, wavelength
dphi = 2 * np.pi * d * y / (lam * np.sqrt(L ** 2 + y ** 2))
env = np.exp(-(y ** 2) / 14)
I = env * np.cos(dphi / 2) ** 2
I = I / I.max()

ax1.plot(y, I, color=BLUE, lw=1.4)
ax1.fill_between(y, 0, I, color=BLUE, alpha=0.12, lw=0)
ax1.set_ylabel("wave intensity $|\\psi_e|^2$")
ax1.set_yticks([])
ax1.text(0, 1.06, "amplitudes of the two paths added, then squared",
         ha="center", fontsize=7.5)

rng = np.random.default_rng(20260802)
p = I / I.sum()
clicks = rng.choice(y, size=362, p=p)
ax2.hist(clicks, bins=48, range=(-6, 6), color=ORANGE, edgecolor="white",
         linewidth=0.4)
ax2.set_ylabel("clicks per bin")
ax2.set_xlabel("position on the far screen")
ax2.text(0, ax2.get_ylim()[1] * 1.04,
         "362 one-atom clicks, arriving one at a time",
         ha="center", fontsize=7.5)
for ax in (ax1, ax2):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig("fig_doubleslit.pdf")
print("wrote fig_doubleslit.pdf")
