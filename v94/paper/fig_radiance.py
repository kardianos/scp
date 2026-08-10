#!/usr/bin/env python3
"""fig_radiance: intake vs graded outtake O = k cap x^4 and the restoring fixed point."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

x = np.linspace(0, 1.0, 400)
k_rad, cap = 0.05, 2.5
O = k_rad * cap * x ** 4
# illustrative intake: rises then saturates (door opens with load, supply-limited)
I = 0.021 * (1 - np.exp(-4.5 * x))

fig, ax = plt.subplots(figsize=(3.4, 2.4))
ax.plot(x, O, color=ORANGE, lw=1.6)
ax.plot(x, I, color=BLUE, lw=1.6)
xs = x[np.argmin(np.abs(O - I)[100:]) + 100]
ys = k_rad * cap * xs ** 4
ax.plot([xs], [ys], marker="o", ms=5, color=INK, zorder=5)
ax.annotate(r"fixed point $x^*$: intake = outtake," "\n"
            "and restoring (outtake steeper)",
            xy=(xs, ys), xytext=(0.16, 0.078), fontsize=7.5,
            arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8))
ax.annotate("", xy=(xs - 0.015, ys * 0.4), xytext=(xs - 0.14, ys * 0.4),
            arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))
ax.annotate("", xy=(xs + 0.02, ys * 1.75), xytext=(xs + 0.15, ys * 1.75),
            arrowprops=dict(arrowstyle="->", color=GRAY, lw=1.0))
ax.text(0.62, 0.101, r"outtake $O = k_{\mathrm{rad}}\,\mathrm{cap}\,x^{4}$",
        color=ORANGE, fontsize=8)
ax.text(0.60, 0.028, "intake $I(x)$ (door opens\nwith load, supply-limited)",
        color=BLUE, fontsize=8)
ax.set_xlabel(r"load $x$ (how full the cell is)")
ax.set_ylabel("energy per unit time")
ax.set_xlim(0, 1); ax.set_ylim(0, 0.13)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig("fig_radiance.pdf")
print("wrote fig_radiance.pdf")
