#!/usr/bin/env python3
"""fig_gravity: the space-sector footprint around a mass; measured facts annotated."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

r = np.linspace(0, 10, 500)
Es0 = 1.0
core = 1.2
# graded depression: deep at the mass, healing outward over a few cells
Es = Es0 - 0.55 * np.exp(-r / 1.8)

fig, ax = plt.subplots(figsize=(3.4, 2.5))
ax.plot(r, Es, color=BLUE, lw=1.7)
ax.axhline(Es0, color=GRAY, lw=0.8, linestyle=(0, (2, 2)))
ax.text(6.4, Es0 + 0.015, r"undisturbed space level $E_{s0}$", color=GRAY,
        fontsize=7.5)
ax.axvspan(0, core, color=ORANGE, alpha=0.14, lw=0)
ax.text(core / 2, 0.40, "mass\n(converted\nspace)", ha="center", color=ORANGE,
        fontsize=7.5)
for rr in (2.2, 3.2, 4.4):
    ax.annotate("", xy=(rr - 0.55, np.interp(rr - 0.55, r, Es) + 0.02),
                xytext=(rr, np.interp(rr, r, Es) + 0.02),
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.0))
ax.text(2.5, 0.62, "pressure pushes inward\n(nothing 'reaches out')",
        fontsize=7.5)
ax.annotate("graded footprint (g1),\nmaintained while the mass lives",
            xy=(1.7, np.interp(1.7, r, Es)), xytext=(2.6, 0.40), fontsize=7.5,
            arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8))
ax.text(9.7, 0.76, "far field: no steady $1/r$ tail\nmeasured (g4) — open",
        fontsize=7.5, ha="right")
ax.set_xlabel(r"distance from the mass (cells)")
ax.set_ylabel(r"space energy $E_s$")
ax.set_xlim(0, 10); ax.set_ylim(0.3, 1.12)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig("fig_gravity.pdf")
print("wrote fig_gravity.pdf")
