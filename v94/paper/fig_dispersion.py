#!/usr/bin/env python3
"""fig_dispersion: polariton avoided crossing of the door pair."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY,
                     "axes.linewidth": 0.7})

w1, w2, Kbar, g = 1.65, 2.9, 0.45, 0.30
k = np.linspace(0, np.pi, 400)
eps = 2 * Kbar * (1 - np.cos(k))          # field band
Of = w1 + eps                             # bare field branch
Om = np.full_like(k, w2)                  # bare register line (flat)
mid = 0.5 * (Of + Om)
gap = np.sqrt(0.25 * (Of - Om) ** 2 + g ** 2)
wp, wm = mid + gap, mid - gap

fig, ax = plt.subplots(figsize=(3.4, 2.7))
ax.plot(k, Of, color=BLUE, lw=0.9, linestyle=(0, (2, 2)))
ax.plot(k, Om, color=GRAY, lw=0.9, linestyle=(0, (2, 2)))
ax.plot(k, wp, color=INK, lw=1.6)
ax.plot(k, wm, color=INK, lw=1.6)
kc = np.arccos(1 - (w2 - w1) / (2 * Kbar))
ax.annotate("", xy=(kc, w2 + g), xytext=(kc, w2 - g),
            arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.2))
ax.text(kc + 0.06, w2 + 0.02, r"gap $2g$", color=ORANGE, fontsize=8)
ax.text(0.35, 1.78, r"bare field  $\Omega_1{+}\epsilon(k)$", color=BLUE,
        fontsize=8, rotation=14)
ax.text(0.012, 2.96, r"bare register  $\Omega_2$ (flat: $\chi$ does not hop)",
        color=GRAY, fontsize=8)
ax.text(2.35, 3.62, r"$\omega_+$", fontsize=9)
ax.text(2.55, 2.71, r"$\omega_-$", fontsize=9)
ax.annotate("slope near $k{=}0$\n$= c_{\\mathrm{sub}}$", xy=(0.42, 1.85),
            xytext=(1.15, 1.28), fontsize=7.5,
            arrowprops=dict(arrowstyle="->", color=GRAY, lw=0.8))
ax.text(2.0, 2.30, "mixed branch:\npart field, part register\n(slower — refraction)",
        fontsize=7.5, ha="center")
ax.set_xlabel(r"wavenumber $k$ (per link)")
ax.set_ylabel(r"frequency $\omega$")
ax.set_xlim(0, np.pi); ax.set_ylim(1.3, 4.0)
ax.set_xticks([0, np.pi / 2, np.pi]); ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("fig_dispersion.pdf")
print("wrote fig_dispersion.pdf")
