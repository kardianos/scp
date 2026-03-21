#!/usr/bin/env python3
"""
1D Radial Quantum Time Evolution: Watch the electron orbital form.

Uses V_eff(r) from Phase 1, evolves a wave packet in the Schrödinger
equation with ℏ_eff. Shows the packet settling into a bound state
at the Bohr radius instead of collapsing classically.

Output: r_θ(t) trajectory, |ψ(r)|² at key times, frequency spectrum.
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal
import os

# Physical parameters from V34/V35
R_BRAID = 5.0
V_DEPTH = 1.27  # total well depth (gravity + EM)
N_POWER = 1.189  # depletion power law

# ℏ_eff from braid action scale
E_BRAID = 5000.0
OMEGA_BRAID = 0.22
HBAR_EFF = E_BRAID / OMEGA_BRAID  # ≈ 22727

# Effective mass — scan to find hydrogen match
# From param sweep: need hbar^2/m_eff ≈ 1.1M-5.5M for Bohr ratio 53000
# At hbar=22727: m_eff = hbar^2 / C where C gives Bohr ratio 53000
# Bohr ratio = a0/r_braid, a0 = hbar^2/(m*alpha), alpha = V_depth
# 53000 = hbar^2/(m * 1.27 * 5) → m = hbar^2 / (53000 * 1.27 * 5)
M_EFF = HBAR_EFF**2 / (53000 * V_DEPTH * R_BRAID)
print(f"Parameters:")
print(f"  ℏ_eff = E_braid/ω = {E_BRAID}/{OMEGA_BRAID} = {HBAR_EFF:.0f}")
print(f"  m_eff = ℏ²/(Bohr_ratio × α × r_braid) = {M_EFF:.1f}")
print(f"  m_eff/m_braid = {M_EFF/E_BRAID:.4f} = 1/{E_BRAID/M_EFF:.0f}")
print()

def V_eff(r, l=0):
    """Effective potential: gravitational + centrifugal."""
    V = np.where(r < 2, V_DEPTH * 10 * (2/r - 1),  # repulsive core
                 -V_DEPTH / (r/R_BRAID)**N_POWER)
    if l > 0:
        V += HBAR_EFF**2 * l*(l+1) / (2*M_EFF * r**2)
    return V

def solve_eigenstates(l=0, r_max=None, N_r=10000):
    """Find bound eigenstates."""
    if r_max is None:
        a0 = HBAR_EFF**2 / (M_EFF * V_DEPTH)
        r_max = max(500, 15 * a0)

    r = np.linspace(0.5, r_max, N_r)
    dr = r[1] - r[0]
    V = V_eff(r, l)

    kinetic = HBAR_EFF**2 / (2 * M_EFF * dr**2)
    diag = 2*kinetic + V
    offdiag = -kinetic * np.ones(N_r - 1)

    n_eig = min(20, N_r - 2)
    energies, vectors = eigh_tridiagonal(diag, offdiag,
                                          select='i', select_range=(0, n_eig-1))

    bound = energies < 0
    return r, energies[bound], vectors[:, bound]

def time_evolve(l=0, r_max=None, N_r=10000, T=None, N_t=2000):
    """Time-evolve a Gaussian wave packet in V_eff."""
    if r_max is None:
        a0 = HBAR_EFF**2 / (M_EFF * V_DEPTH)
        r_max = max(500, 15 * a0)
    if T is None:
        # One orbital period ~ 2π × a0 / v_orbital
        a0 = HBAR_EFF**2 / (M_EFF * V_DEPTH)
        v_orbital = HBAR_EFF / (M_EFF * a0)
        T_orbital = 2 * np.pi * a0 / v_orbital
        T = 5 * T_orbital

    r = np.linspace(0.5, r_max, N_r)
    dr = r[1] - r[0]
    dt = T / N_t

    print(f"Time evolution: r_max={r_max:.0f}, N_r={N_r}, T={T:.0f}, dt={dt:.2f}")
    print(f"  Bohr radius a₀ = {HBAR_EFF**2/(M_EFF*V_DEPTH):.0f}")
    print()

    # Initialize: Gaussian at 0.8 × Bohr radius (slightly inside)
    a0 = HBAR_EFF**2 / (M_EFF * V_DEPTH)
    r0 = 0.8 * a0
    sigma = 0.3 * a0
    k0 = 0  # no initial radial velocity

    psi = np.exp(-(r - r0)**2 / (2*sigma**2)) * np.exp(1j * k0 * r)
    psi /= np.sqrt(np.sum(np.abs(psi)**2) * dr)

    V = V_eff(r, l)

    # Crank-Nicolson time stepping
    kinetic = HBAR_EFF**2 / (2 * M_EFF * dr**2)
    alpha_cn = 1j * dt / (2 * HBAR_EFF)

    # Tridiagonal matrices for CN: (I + iH dt/2ℏ)ψ(t+dt) = (I - iH dt/2ℏ)ψ(t)
    diag_H = 2*kinetic + V
    off_H = -kinetic

    # Record trajectory
    record_every = max(1, N_t // 500)
    times = []
    mean_r = []
    rms_r = []
    prob_density_snapshots = []
    snap_times = [0, int(N_t*0.1), int(N_t*0.25), int(N_t*0.5), int(N_t*0.75), N_t-1]

    for step in range(N_t):
        # Split-step method (simpler than CN, good enough)
        # Half potential step
        psi *= np.exp(-1j * V * dt / (2*HBAR_EFF))

        # Full kinetic step (FFT)
        psi_k = np.fft.fft(psi)
        k = np.fft.fftfreq(N_r, d=dr) * 2 * np.pi
        kinetic_phase = np.exp(-1j * HBAR_EFF * k**2 * dt / (2*M_EFF))
        psi = np.fft.ifft(psi_k * kinetic_phase)

        # Half potential step
        psi *= np.exp(-1j * V * dt / (2*HBAR_EFF))

        # Normalize (absorbing boundary)
        norm = np.sum(np.abs(psi)**2) * dr
        if norm > 0:
            psi /= np.sqrt(norm)

        if step % record_every == 0 or step in snap_times:
            prob = np.abs(psi)**2
            mr = np.sum(r * prob) * dr
            rr = np.sqrt(np.sum(r**2 * prob) * dr)
            times.append(step * dt)
            mean_r.append(mr)
            rms_r.append(rr)

        if step in snap_times:
            prob_density_snapshots.append((step * dt, r.copy(), np.abs(psi)**2))

        if step % (N_t // 10) == 0:
            prob = np.abs(psi)**2
            mr = np.sum(r * prob) * dr
            print(f"  t={step*dt:10.0f}  <r>={mr:10.0f}  (Bohr={a0:.0f})")

    return np.array(times), np.array(mean_r), np.array(rms_r), prob_density_snapshots


def main():
    os.makedirs('data', exist_ok=True)

    # First: find eigenstates
    print("=== Eigenstates ===\n")
    r, E, V = solve_eigenstates(l=0)
    a0 = HBAR_EFF**2 / (M_EFF * V_DEPTH)

    print(f"Bohr radius a₀ = {a0:.0f} code units")
    print(f"a₀ / r_braid = {a0/R_BRAID:.0f}")
    print(f"\nBound states (l=0):")
    for i in range(min(5, len(E))):
        prob = V[:, i]**2
        norm = np.sum(prob) * (r[1]-r[0])
        mr = np.sum(r * prob) * (r[1]-r[0]) / norm
        print(f"  n={i+1}: E={E[i]:.6e}  <r>={mr:.0f}  <r>/r_braid={mr/R_BRAID:.0f}")

    if len(E) >= 2:
        print(f"\nE₂/E₁ = {E[1]/E[0]:.4f}  (hydrogen: 0.2500)")

    # Angular momentum states
    print(f"\nAngular momentum (at this ℏ, m_eff):")
    for l in [0, 1, 2, 3]:
        _, El, _ = solve_eigenstates(l=l)
        if len(El) > 0:
            print(f"  l={l}: {len(El)} bound states, E₀={El[0]:.6e}")
        else:
            print(f"  l={l}: no bound states")

    # Time evolution
    print(f"\n=== Time Evolution ===\n")
    times, mean_r, rms_r, snaps = time_evolve(l=0)

    # Save trajectory
    with open('data/orbital_trajectory.tsv', 'w') as f:
        f.write("t\tmean_r\trms_r\tmean_r_over_a0\n")
        for i in range(len(times)):
            f.write(f"{times[i]:.2f}\t{mean_r[i]:.2f}\t{rms_r[i]:.2f}\t{mean_r[i]/a0:.4f}\n")

    # Save snapshots
    with open('data/orbital_snapshots.tsv', 'w') as f:
        f.write("t\tr\tprob_density\n")
        for t, r_snap, prob in snaps:
            for i in range(0, len(r_snap), max(1, len(r_snap)//500)):
                f.write(f"{t:.2f}\t{r_snap[i]:.2f}\t{prob[i]:.6e}\n")

    # Summary
    print(f"\n=== Summary ===")
    print(f"  ℏ_eff = {HBAR_EFF:.0f} (from E_braid/ω)")
    print(f"  m_eff = {M_EFF:.1f} (from Bohr ratio constraint)")
    print(f"  m_eff/m_braid = 1/{E_BRAID/M_EFF:.0f}")
    print(f"  Bohr radius = {a0:.0f} code units = {a0/R_BRAID:.0f} × r_braid")
    print(f"  Bound states: {len(E)} (l=0)")
    if len(E) >= 2:
        print(f"  E₂/E₁ = {E[1]/E[0]:.4f} (hydrogen: 0.250)")
    print(f"  Final <r> = {mean_r[-1]:.0f} (Bohr = {a0:.0f})")
    print(f"  <r>/a₀ = {mean_r[-1]/a0:.3f} (should oscillate around 1.0)")


if __name__ == '__main__':
    main()
