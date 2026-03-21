#!/usr/bin/env python3
"""
V35 Phase 2: Radial eigenvalue problem for bound states.

Solves:  -hbar^2/(2m) u''(r) + [V_eff(r) + hbar^2 l(l+1)/(2m r^2)] u = E u

where u(r) = r * R(r) is the reduced radial wavefunction satisfying u(0) = u(inf) = 0.

Two approaches:
  A. Uniform r-grid (simple, good for moderate range)
  B. Log-grid with generalized eigenvalue (for large range)

Uses scipy.linalg.eigh_tridiagonal for efficiency.
"""

import numpy as np
import os
import sys
import time
from scipy.linalg import eigh_tridiagonal

# ── Potential from Phase 1 ─────────────────────────────────────────────────

B_grav = 0.3628
n_grav = 1.1891
V0_grav = B_grav / 5.0**n_grav
em_fraction = 0.27
total_factor = 1.0 + em_fraction  # 1.27

r_braid = 5.0
r_core = 2.0

def V_well(r):
    """Total effective potential (gravity + EM). V(5) = -1.27."""
    V = np.zeros_like(r, dtype=float)

    outer = r >= r_braid
    V[outer] = -total_factor * B_grav / (r[outer]**n_grav * V0_grav)

    inner = (r < r_braid) & (r > 0)
    C_core = 1.0 / (r_braid/r_core - 1.0)**2
    V[inner] = total_factor * (-1.0 + C_core * (r_braid/r[inner] - 1.0)**2)

    return V

# ── Approach A: Uniform r-grid ─────────────────────────────────────────────

def solve_uniform(N_r, r_min, r_max, hbar, m_eff, l_ang, n_eig=30):
    """
    Solve on uniform r-grid.
    u(r) = r*R(r), equation: -hbar^2/(2m) u'' + V_eff(r) u = E u
    V_eff = V_well(r) + hbar^2 l(l+1)/(2m r^2)
    Boundary: u(r_min) = u(r_max) = 0
    """
    r = np.linspace(r_min, r_max, N_r)
    dr = r[1] - r[0]
    h2_2m = hbar**2 / (2.0 * m_eff)

    # Effective potential
    V = V_well(r)
    V_cent = h2_2m * l_ang * (l_ang + 1) / r**2
    V_eff = V + V_cent

    # Kinetic: -hbar^2/(2m) d^2/dr^2 => tridiagonal with coefficient h2_2m/dr^2
    kin_coeff = h2_2m / dr**2

    # Diagonal: 2*kin_coeff + V_eff
    diag = 2.0 * kin_coeff + V_eff

    # Off-diagonal: -kin_coeff
    off = np.full(N_r - 1, -kin_coeff)

    # Solve for lowest eigenvalues
    n_req = min(n_eig, N_r - 2)
    eigenvalues, eigenvectors = eigh_tridiagonal(
        diag, off, select='i', select_range=(0, n_req - 1)
    )

    return eigenvalues, eigenvectors, r, dr

# ── Approach B: Log-grid with generalized eigenvalue ───────────────────────

def solve_loggrid(N_r, r_min, r_max, hbar, m_eff, l_ang, n_eig=30):
    """
    Solve on logarithmic grid using substitution.

    Let x = ln(r/r_min), so r = r_min * e^x.
    Let u(r) = r*R(r) = g(x) * r^{1/2}  (to remove first-derivative)

    Actually, simplest correct approach:
    - Uniform grid in x = ln(r), dx = const
    - u(r) = r*R(r) satisfies:
        -hbar^2/(2m) u'' + V_eff(r) u = E u
    - In x-coordinates with u(r(x)):
        du/dr = du/dx * dx/dr = (1/r) du/dx
        d^2u/dr^2 = (1/r^2)(d^2u/dx^2 - du/dx)

    So the equation becomes:
        -hbar^2/(2m) * (1/r^2) * (u_xx - u_x) + V_eff u = E u

    This has both u_xx and u_x terms. Substitute f(x) = u(r(x)):
        -hbar^2/(2m r^2) (f'' - f') + V_eff f = E f

    Multiply through by r^2:
        -hbar^2/(2m) (f'' - f') + r^2 V_eff f = E r^2 f

    This is a generalized eigenvalue problem with a first-derivative term.
    Handle the first-derivative by further substitution: f = g * e^{x/2} = g * sqrt(r)

    Then f' = g'*sqrt(r) + g/(2*sqrt(r)) = (g' + g/2) * e^{x/2}
    f'' = (g'' + g') * e^{x/2} + (g' + g/2)/2 * e^{x/2}
         Wait, let me be more careful.

    f(x) = g(x) * e^{x/2}
    f' = (g' + g/2) e^{x/2}
    f'' = (g'' + g' + g'/2 + g/4) e^{x/2} = (g'' + 3g'/2 + g/4) e^{x/2}

    Hmm, this doesn't cleanly remove f'. Let me use the standard approach instead.

    Standard approach: multiply the original equation by r^2 to get generalized EV,
    and discretize f'', f' on the uniform x-grid.
    """
    x_min = np.log(r_min)
    x_max = np.log(r_max)
    x = np.linspace(x_min, x_max, N_r)
    dx = x[1] - x[0]
    r = np.exp(x)

    h2_2m = hbar**2 / (2.0 * m_eff)

    V = V_well(r)
    V_cent = h2_2m * l_ang * (l_ang + 1) / r**2
    V_eff = V + V_cent

    # Equation in x: -h2_2m/r^2 (f'' - f') + V_eff f = E f
    # Multiply by r^2: -h2_2m (f'' - f') + r^2 V_eff f = E r^2 f
    #
    # H f = E S f where
    # H = -h2_2m (D2 - D1) + diag(r^2 * V_eff)
    # S = diag(r^2)
    #
    # D2 = second derivative: (f_{i+1} - 2f_i + f_{i-1}) / dx^2
    # D1 = first derivative: (f_{i+1} - f_{i-1}) / (2*dx)
    #
    # H is NOT tridiagonal-symmetric because of D1. But we can make it symmetric
    # with the substitution f(x) = g(x) * e^{-x/2}:
    # f' = (g' - g/2) e^{-x/2}
    # f'' = (g'' - g' + g'/2 - g/4) ... NO this gets messy.
    #
    # Better: use the substitution that makes the operator self-adjoint.
    # The original radial equation for u(r) is already self-adjoint.
    # On a UNIFORM r-grid, it's symmetric tridiagonal. The log grid
    # introduces asymmetry. For a quick solution, just use a uniform grid
    # with enough points.
    #
    # For large r_max, we need adaptive: uniform up to some r_cut, then sparse.
    # OR: use Numerov's method on log grid (shooting).

    # Fall back to uniform for correctness
    return None

# ── Shooting method for verification ──────────────────────────────────────

def shoot(hbar, m_eff, l_ang, E_trial, r_min=0.5, r_max=1000, N=10000):
    """
    Numerov method to integrate the radial equation at given E.
    Returns the wavefunction value at the matching point.
    """
    r = np.linspace(r_min, r_max, N)
    dr = r[1] - r[0]
    h2_2m = hbar**2 / (2.0 * m_eff)

    V = V_well(r)
    V_cent = h2_2m * l_ang * (l_ang + 1) / r**2

    # f(r) = 2m/hbar^2 * (V_eff - E)
    f = (V + V_cent - E_trial) / h2_2m if h2_2m > 0 else np.zeros_like(r)

    # Numerov: u_{n+1} = [2(1 - 5/12 dr^2 f_n) u_n - (1 + 1/12 dr^2 f_{n-1}) u_{n-1}]
    #                     / (1 + 1/12 dr^2 f_{n+1})

    # Outward integration from r_min
    u = np.zeros(N)
    u[0] = 0.0
    u[1] = dr  # small starting value

    for i in range(1, N-1):
        num = 2*(1 - 5.0/12.0 * dr**2 * f[i]) * u[i] \
            - (1 + 1.0/12.0 * dr**2 * f[i-1]) * u[i-1]
        den = 1 + 1.0/12.0 * dr**2 * f[i+1]
        if abs(den) < 1e-30:
            break
        u[i+1] = num / den
        # Prevent overflow
        if abs(u[i+1]) > 1e30:
            u[i+1:] = np.sign(u[i+1]) * 1e30
            break

    return u, r

# ── Multi-range uniform grid solver ───────────────────────────────────────

def solve_adaptive(hbar, m_eff, l_ang, n_eig=30):
    """
    Solve using uniform grid with adaptive range.
    For small hbar: wavefunction is near the well, use small r_max.
    For large hbar: wavefunction extends far, use large r_max.

    The Bohr radius scales as hbar^2/(m*alpha), so r_max should scale similarly.
    """
    alpha_eff = 1.27  # well depth

    # Estimate Bohr radius
    a0_est = hbar**2 / (m_eff * alpha_eff)

    # Grid should extend to ~10 * a0 (or at least past the well)
    r_max = max(50.0, 20.0 * a0_est)
    r_max = min(r_max, 5e6)  # cap to prevent memory issues

    # Grid spacing: need dr << hbar / sqrt(2*m*|V_min|) (de Broglie wavelength at well bottom)
    lambda_dB = hbar / np.sqrt(2.0 * m_eff * abs(alpha_eff)) if alpha_eff > 0 else hbar
    dr_target = min(lambda_dB / 10.0, 0.1)
    dr_target = max(dr_target, r_max / 500000)  # limit N to 500k

    N_r = int(r_max / dr_target) + 1
    N_r = max(N_r, 1000)
    N_r = min(N_r, 500000)

    r_min = dr_target  # small but nonzero

    return solve_uniform(N_r, r_min, r_max, hbar, m_eff, l_ang, n_eig), r_max, N_r

# ── Main scan ──────────────────────────────────────────────────────────────

def main():
    m_eff = 1.0
    l_ang = 0

    hbar_values = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0,
                   10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0]

    print("=" * 100)
    print("V35 PHASE 2: Radial Eigenvalue Problem — Bound State Search")
    print("=" * 100)
    print(f"  m_eff = {m_eff}")
    print(f"  l = {l_ang}")
    print(f"  Potential: V(5) = {V_well(np.array([5.0]))[0]:.4f} (well minimum)")
    print(f"  Power law: V ~ -1.27 / r^{n_grav:.3f} for r > 5")
    print(f"  Core: repulsive for r < {r_core}")
    print()

    # Hydrogen reference
    alpha = 1.27
    print("  Hydrogen-like reference (V = -alpha/r, alpha=1.27):")
    print(f"    a_0 = hbar^2/(m*alpha)")
    print(f"    E_0 = -m*alpha^2/(2*hbar^2)")
    print()

    results = []

    print(f"{'hbar':>10s}  {'r_max':>10s}  {'N_r':>8s}  {'dr':>8s}  "
          f"{'N_bnd':>5s}  {'E_0':>14s}  {'<r>_0':>12s}  {'<r>/r_b':>8s}  "
          f"{'E_1':>14s}  {'<r>_1':>12s}  {'t':>5s}")
    print("-" * 120)

    for hbar in hbar_values:
        t0 = time.time()

        sol, r_max_used, N_r_used = solve_adaptive(hbar, m_eff, l_ang)
        eigenvalues, eigenvectors, r, dr = sol

        dt = time.time() - t0

        # Count bound states
        bound_mask = eigenvalues < 0
        E_bound = eigenvalues[bound_mask]
        n_bound = len(E_bound)

        if n_bound > 0:
            psi0 = eigenvectors[:, np.where(bound_mask)[0][0]]
            # <r> = integral |u(r)|^2 r dr / integral |u(r)|^2 dr
            # (u = r*R, so |u|^2 dr is already the radial probability density on r)
            # Wait: for the reduced wavefunction u(r), |u(r)|^2 dr gives the
            # probability of finding particle between r and r+dr.
            # So <r> = int r |u|^2 dr / int |u|^2 dr
            prob0 = psi0**2
            norm0 = np.sum(prob0) * dr
            r0_mean = np.sum(prob0 * r) * dr / norm0 if norm0 > 0 else 0
            r0_rms = np.sqrt(np.sum(prob0 * r**2) * dr / norm0) if norm0 > 0 else 0
            E0 = E_bound[0]

            if n_bound > 1:
                psi1 = eigenvectors[:, np.where(bound_mask)[0][1]]
                prob1 = psi1**2
                norm1 = np.sum(prob1) * dr
                r1_mean = np.sum(prob1 * r) * dr / norm1 if norm1 > 0 else 0
                E1 = E_bound[1]
                E1_str = f"{E1:14.6e}"
                r1_str = f"{r1_mean:12.2f}"
            else:
                E1, r1_mean = None, None
                E1_str = "---"
                r1_str = "---"

            print(f"{hbar:10.3f}  {r_max_used:10.0f}  {N_r_used:8d}  {dr:8.4f}  "
                  f"{n_bound:5d}  {E0:14.6e}  {r0_mean:12.2f}  {r0_mean/r_braid:8.1f}  "
                  f"{E1_str:>14s}  {r1_str:>12s}  {dt:4.1f}s")

            results.append({
                'hbar': hbar, 'n_bound': n_bound,
                'E0': E0, 'r0_mean': r0_mean, 'r0_rms': r0_rms,
                'E1': E1, 'r1_mean': r1_mean,
                'r_max': r_max_used, 'N_r': N_r_used, 'dr': dr
            })
        else:
            print(f"{hbar:10.3f}  {r_max_used:10.0f}  {N_r_used:8d}  {dr:8.4f}  "
                  f"{0:5d}  {'---':>14s}  {'---':>12s}  {'---':>8s}  "
                  f"{'---':>14s}  {'---':>12s}  {dt:4.1f}s")
            results.append({
                'hbar': hbar, 'n_bound': 0,
                'E0': None, 'r0_mean': None, 'r0_rms': None,
                'E1': None, 'r1_mean': None,
                'r_max': r_max_used, 'N_r': N_r_used, 'dr': dr
            })

    # ── Verification with pure Coulomb ─────────────────────────────────────

    print()
    print("=" * 100)
    print("VERIFICATION: Pure Coulomb V = -1.27/r (should give E_n = -0.8065/(hbar*n)^2)")
    print("=" * 100)
    print()

    def V_coulomb(r):
        """Pure Coulomb for verification."""
        return -1.27 / r

    hbar_test = [1.0, 10.0, 100.0]
    for hbar in hbar_test:
        a0 = hbar**2 / (m_eff * 1.27)
        r_max_c = max(50, 20 * a0)
        r_max_c = min(r_max_c, 5e6)
        lambda_dB = hbar / np.sqrt(2.0 * 1.27)
        dr_c = min(lambda_dB / 10.0, 0.1)
        dr_c = max(dr_c, r_max_c / 500000)
        N_c = int(r_max_c / dr_c) + 1
        N_c = max(N_c, 1000)
        N_c = min(N_c, 500000)

        r_c = np.linspace(dr_c, r_max_c, N_c)
        dr_c = r_c[1] - r_c[0]
        h2_2m = hbar**2 / (2.0 * m_eff)

        V_c = V_coulomb(r_c)
        diag_c = 2.0 * h2_2m / dr_c**2 + V_c
        off_c = np.full(N_c - 1, -h2_2m / dr_c**2)

        evals_c, evecs_c = eigh_tridiagonal(
            diag_c, off_c, select='i', select_range=(0, 4)
        )

        E_exact = -m_eff * 1.27**2 / (2 * hbar**2)
        bound_c = evals_c[evals_c < 0]

        if len(bound_c) > 0:
            prob = evecs_c[:, 0]**2
            norm = np.sum(prob) * dr_c
            r_mean = np.sum(prob * r_c) * dr_c / norm

            print(f"  hbar={hbar:6.1f}: E_0 = {bound_c[0]:14.6e} (exact: {E_exact:14.6e}, "
                  f"err={abs(bound_c[0]-E_exact)/abs(E_exact)*100:.2f}%), "
                  f"<r> = {r_mean:.2f} (a_0 = {a0:.2f}, err={abs(r_mean-1.5*a0)/(1.5*a0)*100:.1f}%), "
                  f"N_bound={len(bound_c)}, N_r={N_c}, r_max={r_max_c:.0f}")
        else:
            print(f"  hbar={hbar:6.1f}: NO bound states found (N_r={N_c}, r_max={r_max_c:.0f})")

    # ── Fine scan for interesting regime ───────────────────────────────────

    print()
    print("=" * 100)
    print("FINE SCAN: transition region")
    print("=" * 100)
    print()

    # Find where N_bound transitions
    transition_hbar = []
    for i in range(len(results)-1):
        n0 = results[i]['n_bound']
        n1 = results[i+1]['n_bound']
        if n0 != n1 and n0 > 0 and n1 > 0:
            # Different number of bound states
            h0 = results[i]['hbar']
            h1 = results[i+1]['hbar']
            transition_hbar.extend(np.logspace(np.log10(h0), np.log10(h1), 8))
        if (n0 > 0 and n1 == 0) or (n0 == 0 and n1 > 0):
            h0 = results[i]['hbar']
            h1 = results[i+1]['hbar']
            transition_hbar.extend(np.logspace(np.log10(h0), np.log10(h1), 15))

    # Also scan regime where <r>/r_braid might be ~50000 (hydrogen-like)
    # a_0 = hbar^2/(m*alpha) = 50000*5 => hbar = sqrt(250000*1.27) ~ 563
    target_hbar = np.sqrt(50000 * r_braid * m_eff * 1.27)
    transition_hbar.extend([target_hbar * x for x in [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]])

    transition_hbar = sorted(set([h for h in transition_hbar if 0.001 <= h <= 50000]))

    fine_results = []
    if transition_hbar:
        print(f"{'hbar':>10s}  {'N_bnd':>5s}  {'E_0':>14s}  {'<r>_0':>12s}  {'<r>/r_b':>8s}  "
              f"{'E_0(H)':>14s}  {'a_0(H)':>12s}")
        print("-" * 85)

        for hbar in transition_hbar:
            sol, r_max_u, N_r_u = solve_adaptive(hbar, m_eff, l_ang)
            eigenvalues, eigenvectors, r, dr = sol

            bound_mask = eigenvalues < 0
            E_bound = eigenvalues[bound_mask]
            n_bound = len(E_bound)

            # Hydrogen reference
            E_H = -m_eff * 1.27**2 / (2 * hbar**2)
            a0_H = hbar**2 / (m_eff * 1.27)

            if n_bound > 0:
                psi0 = eigenvectors[:, np.where(bound_mask)[0][0]]
                prob0 = psi0**2
                norm0 = np.sum(prob0) * dr
                r0_mean = np.sum(prob0 * r) * dr / norm0 if norm0 > 0 else 0

                print(f"{hbar:10.3f}  {n_bound:5d}  {E_bound[0]:14.6e}  {r0_mean:12.2f}  "
                      f"{r0_mean/r_braid:8.1f}  {E_H:14.6e}  {a0_H:12.2f}")
                fine_results.append({
                    'hbar': hbar, 'n_bound': n_bound,
                    'E0': E_bound[0], 'r0_mean': r0_mean, 'r0_rms': 0
                })
            else:
                print(f"{hbar:10.3f}  {0:5d}  {'---':>14s}  {'---':>12s}  {'---':>8s}  "
                      f"{E_H:14.6e}  {a0_H:12.2f}")
                fine_results.append({
                    'hbar': hbar, 'n_bound': 0, 'E0': None, 'r0_mean': None, 'r0_rms': None
                })

    # ── l > 0 modes ────────────────────────────────────────────────────────

    print()
    print("=" * 100)
    print("ANGULAR MOMENTUM: l = 0, 1, 2, 3 at key hbar values")
    print("=" * 100)
    print()

    # Pick hbar values where l=0 has bound states
    l_scan_hbar = []
    for res in results + fine_results:
        if res.get('n_bound', 0) > 0 and res.get('n_bound', 0) <= 30:
            l_scan_hbar.append(res['hbar'])
    # Also add key values
    l_scan_hbar.extend([1.0, 10.0, 100.0, 563.0])
    l_scan_hbar = sorted(set(l_scan_hbar))
    # Keep only a few representative values
    if len(l_scan_hbar) > 8:
        indices = np.linspace(0, len(l_scan_hbar)-1, 8, dtype=int)
        l_scan_hbar = [l_scan_hbar[i] for i in indices]

    print(f"{'hbar':>10s}  {'l':>3s}  {'N_bnd':>5s}  {'E_0':>14s}  {'<r>_0':>12s}")
    print("-" * 55)

    for hbar in l_scan_hbar:
        for l in [0, 1, 2, 3]:
            sol, _, _ = solve_adaptive(hbar, m_eff, l)
            eigenvalues, eigenvectors, r, dr = sol
            bound_mask = eigenvalues < 0
            E_bound = eigenvalues[bound_mask]
            n_bound = len(E_bound)

            if n_bound > 0:
                psi0 = eigenvectors[:, np.where(bound_mask)[0][0]]
                prob0 = psi0**2
                norm0 = np.sum(prob0) * dr
                r0_mean = np.sum(prob0 * r) * dr / norm0 if norm0 > 0 else 0
                print(f"{hbar:10.3f}  {l:3d}  {n_bound:5d}  {E_bound[0]:14.6e}  {r0_mean:12.2f}")
            else:
                print(f"{hbar:10.3f}  {l:3d}  {0:5d}  {'---':>14s}  {'---':>12s}")

    # ── Summary ────────────────────────────────────────────────────────────

    print()
    print("=" * 100)
    print("SUMMARY: Hydrogen-like scaling")
    print("=" * 100)
    print()
    print("For V ~ -alpha/r^n with n=1.189 (between Coulomb n=1 and Newton n=2):")
    print("  - All n < 2 potentials support bound states")
    print("  - For Coulomb (n=1): E_0 = -m*alpha^2/(2*hbar^2), a_0 = hbar^2/(m*alpha)")
    print(f"  - alpha_eff = {total_factor:.2f} (well depth at r=5)")
    print()
    print("Key ratios:")
    for hbar in [1.0, 10.0, 100.0, 563.0, 1000.0]:
        a0 = hbar**2 / (m_eff * 1.27)
        E0_H = -m_eff * 1.27**2 / (2 * hbar**2)
        print(f"  hbar={hbar:8.1f}:  a_0={a0:12.1f}  E_0={E0_H:14.6e}  <r>/r_braid={a0*1.5/r_braid:10.0f}")

    # ── Save all results ───────────────────────────────────────────────────

    all_results = results + fine_results
    all_results.sort(key=lambda x: x['hbar'])
    # Deduplicate: keep the entry with more info (non-zero r_rms) for each hbar
    seen = {}
    for res in all_results:
        h = f"{res['hbar']:.10e}"
        if h not in seen:
            seen[h] = res
        else:
            # Keep the one with more complete data
            old = seen[h]
            if (old.get('r0_rms') is None or old.get('r0_rms') == 0) and \
               (res.get('r0_rms') is not None and res.get('r0_rms') != 0):
                seen[h] = res
    all_results = sorted(seen.values(), key=lambda x: x['hbar'])

    outfile = os.path.join(os.path.dirname(__file__), 'data', 'bound_states.tsv')
    with open(outfile, 'w') as f:
        f.write("# V35 Phase 2: Bound state eigenvalues\n")
        f.write("# m_eff=1.0, l=0, potential: V=-1.27/r^1.189\n")
        f.write("hbar\tn_bound\tE_0\tr_mean_0\tr_rms_0\tE_1\tr_mean_1\n")
        for res in all_results:
            h = res['hbar']
            nb = res['n_bound']
            E0 = f"{res['E0']:.8e}" if res.get('E0') is not None else "nan"
            r0 = f"{res['r0_mean']:.4f}" if res.get('r0_mean') is not None else "nan"
            rr0 = f"{res['r0_rms']:.4f}" if res.get('r0_rms') is not None else "nan"
            E1 = f"{res.get('E1', None):.8e}" if res.get('E1') is not None else "nan"
            r1 = f"{res.get('r1_mean', None):.4f}" if res.get('r1_mean') is not None else "nan"
            f.write(f"{h:.6e}\t{nb}\t{E0}\t{r0}\t{rr0}\t{E1}\t{r1}\n")

    print()
    print(f"Saved: {outfile}")

if __name__ == '__main__':
    main()
