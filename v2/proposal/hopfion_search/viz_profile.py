#!/usr/bin/env python3
"""
Visualize the Skyrmion hedgehog profile f(r) and energy densities
from the shooting method solver output.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

def load_profile(filename="profile.dat"):
    data = np.loadtxt(filename)
    return {
        'r': data[:, 0],
        'f': data[:, 1],
        'fp': data[:, 2],
        'bdens': data[:, 3],
        'ed2': data[:, 4],
        'ed4': data[:, 5],
    }

def main():
    fname = sys.argv[1] if len(sys.argv) > 1 else "profile.dat"
    d = load_profile(fname)
    r, f, fp = d['r'], d['f'], d['fp']
    bdens, ed2, ed4 = d['bdens'], d['ed2'], d['ed4']

    # Trim to where profile is significant
    rmax_plot = r[np.searchsorted(-f, -0.01)]  # where f < 0.01
    rmax_plot = min(rmax_plot * 1.3, r[-1])
    mask = r <= rmax_plot

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("B=1 Hedgehog Skyrmion Profile (CHPT)", fontsize=14, fontweight='bold')

    # 1) Profile f(r)
    ax = axes[0, 0]
    ax.plot(r[mask], f[mask], 'b-', linewidth=2)
    ax.axhline(y=np.pi, color='gray', linestyle='--', alpha=0.5, label=r'$\pi$')
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('r')
    ax.set_ylabel('f(r)')
    ax.set_title('Profile function f(r)')
    ax.set_ylim(-0.2, np.pi + 0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2) Derivative f'(r)
    ax = axes[0, 1]
    ax.plot(r[mask], fp[mask], 'r-', linewidth=2)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('r')
    ax.set_ylabel("f'(r)")
    ax.set_title("Profile derivative f'(r)")
    ax.grid(True, alpha=0.3)

    # 3) Energy densities (radial, weighted by r^2 for 4pi*r^2*dr measure)
    ax = axes[1, 0]
    # The profile.dat stores the 1D integrand; multiply by r^2 for the shell contribution
    r_safe = np.where(r > 0, r, 1e-30)
    shell_e2 = ed2[mask] * r_safe[mask]**0  # already in 1D form
    shell_e4 = ed4[mask] * r_safe[mask]**0
    ax.plot(r[mask], ed2[mask], 'b-', linewidth=1.5, label=r'$\varepsilon_2$ (gradient)')
    ax.plot(r[mask], ed4[mask], 'r-', linewidth=1.5, label=r'$\varepsilon_4$ (Skyrme)')
    ax.plot(r[mask], ed2[mask] + ed4[mask], 'k--', linewidth=1, label=r'$\varepsilon_{total}$')
    ax.set_xlabel('r')
    ax.set_ylabel('Energy density')
    ax.set_title('Radial energy densities')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4) Baryon density
    ax = axes[1, 1]
    ax.plot(r[mask], bdens[mask], 'g-', linewidth=2)
    ax.fill_between(r[mask], 0, bdens[mask], alpha=0.2, color='green')
    ax.set_xlabel('r')
    ax.set_ylabel(r'$\rho_B(r)$')
    ax.set_title('Baryon charge density')
    ax.grid(True, alpha=0.3)

    # Compute integrated quantities for annotation
    dr = r[1] - r[0]
    E2 = 2 * np.pi * np.trapz(ed2, r)
    E4_coeff = 4 * np.pi / 16.0  # for e=4; this is approximate
    Q = (2.0 / np.pi) * np.trapz(-fp * np.sin(f)**2, r)

    fig.tight_layout()
    outfile = fname.replace('.dat', '.png')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")
    plt.close()

if __name__ == "__main__":
    main()
