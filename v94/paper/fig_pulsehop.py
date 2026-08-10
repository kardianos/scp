#!/usr/bin/env python3
"""fig_pulsehop: a pulse moving along a 12-cell chain over successive time steps
(actual simulation of the precession+exchange update)."""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, INK = "#2E6FB0", "#D9621E", "#8a8a8a", "#222222"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans",
                     "text.color": INK, "axes.edgecolor": GRAY})

# simulate the actual primitive on a chain: psi' = e^{-i O dt} then Givens hops
N, O1, K, dt = 40, 1.65, 0.8, 0.05
x = np.arange(N)
k0 = 1.1
psi = np.exp(-((x - 8.0) ** 2) / 8.0) * np.exp(1j * k0 * x)   # moving packet
snaps, times = [], []
steps_per_snap, nsnap = 60, 4
t = 0.0
for s in range(nsnap):
    snaps.append(np.abs(psi) ** 2); times.append(t)
    for _ in range(steps_per_snap):
        psi = psi * np.exp(-1j * O1 * dt)
        for i in range(N - 1):                      # pairwise Givens sweep
            tau = K * dt
            a, b = psi[i], psi[i + 1]
            psi[i] = np.cos(tau) * a - 1j * np.sin(tau) * b
            psi[i + 1] = np.cos(tau) * b - 1j * np.sin(tau) * a
        t += dt

fig, axes = plt.subplots(nsnap, 1, figsize=(3.4, 2.9), sharex=True)
for s, ax in enumerate(axes):
    E = snaps[s]
    ax.bar(x, E, width=0.85, color=BLUE, edgecolor="white", linewidth=0.3)
    ax.set_ylim(0, 1.15); ax.set_yticks([])
    ax.text(39.3, 0.62, f"$t={times[s]:.0f}$", ha="right", fontsize=7.5)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    pk = np.argmax(E)
    ax.annotate("", xy=(pk, 1.08), xytext=(pk - 0.001, 1.08))
axes[0].text(1, 0.85, "energy per cell $|\\psi|^2$", fontsize=7.5)
axes[-1].set_xlabel("cell number along the chain", fontsize=8)
axes[-1].set_xticks([0, 10, 20, 30, 39])
fig.tight_layout(h_pad=0.3)
fig.savefig("fig_pulsehop.pdf")
print("wrote fig_pulsehop.pdf")
