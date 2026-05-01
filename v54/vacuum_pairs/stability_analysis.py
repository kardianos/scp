#!/usr/bin/env python3
"""
stability_analysis.py — Linear stability analysis for V54 vacuum pair extension.

Computes:
1. Critical amplitude A_crit for pair creation
2. Growth rate gamma(q) as function of wavenumber
3. Saturation amplitude theta_sat
4. Compares predictions to simulation diagnostics

Usage:
  python3 stability_analysis.py                          # just theory
  python3 stability_analysis.py test_stable_diag.tsv     # compare to sim
  python3 stability_analysis.py test_stable_diag.tsv test_unstable_diag.tsv
"""

import sys
import math
import numpy as np

# ============================================================
# Parameters (must match test configs)
# ============================================================
m_phi    = 1.5       # phi mass (m² = 2.25)
m_theta  = 1.0       # theta mass (m_θ² = 1.0)
sigma_c  = 0.5       # cross-potential coupling
lambda_s = 1.0       # theta self-potential
eta      = 0.5       # curl coupling
A_stable = 0.8       # amplitude for stable test
A_unstable = 1.5     # amplitude for unstable test


def stability_report(A, label):
    """Print stability predictions for a given amplitude."""
    m2_theta = m_theta ** 2
    phi2_avg = 1.5 * A ** 2  # <|phi|^2> for 3-field background
    A_crit = math.sqrt(2 * m2_theta / (3 * sigma_c)) if sigma_c > 0 else float('inf')

    print(f"\n{'='*60}")
    print(f"  {label}: A = {A}")
    print(f"{'='*60}")
    print(f"  m_theta² = {m2_theta}")
    print(f"  sigma_cross = {sigma_c}")
    print(f"  lambda_self = {lambda_s}")
    print(f"  <|phi|²> = {phi2_avg:.4f}")
    print(f"  A_crit = {A_crit:.4f}")
    print()

    omega2_q0 = m2_theta - sigma_c * phi2_avg
    print(f"  omega²(q=0) = m_theta² - sigma_cross*<|phi|²>")
    print(f"              = {m2_theta} - {sigma_c}*{phi2_avg:.4f}")
    print(f"              = {omega2_q0:.4f}")
    print()

    if omega2_q0 >= 0:
        omega_q0 = math.sqrt(omega2_q0)
        print(f"  STATUS: STABLE (omega² > 0)")
        print(f"  Theta oscillation frequency: omega = {omega_q0:.4f}")
        print(f"  Theta period: T = {2*math.pi/omega_q0:.2f}")
        print(f"  Prediction: theta_rms stays near zero (< 0.01)")
    else:
        gamma = math.sqrt(-omega2_q0)
        print(f"  STATUS: UNSTABLE (omega² < 0)")
        print(f"  Growth rate: gamma = {gamma:.4f}")
        print(f"  Growth time: tau = {1/gamma:.4f}")
        print(f"  e-folding: theta doubles every {math.log(2)/gamma:.2f} time units")
        print()

        # Saturation
        theta_sat = math.sqrt(-omega2_q0 / lambda_s) if lambda_s > 0 else float('inf')
        print(f"  Saturation: |theta|_sat = sqrt((-omega²)/lambda)")
        print(f"            = sqrt({-omega2_q0:.4f}/{lambda_s})")
        print(f"            = {theta_sat:.4f}")
        print(f"  theta_rms_sat ≈ {theta_sat/math.sqrt(3):.4f} (uniform)")
        print()

        # Time to saturation (from small perturbation)
        # theta(t) ~ theta_0 * exp(gamma*t), saturates when theta ~ theta_sat
        # log(theta_sat/theta_0) / gamma
        theta_0 = 1e-6  # numerical noise level
        t_sat = math.log(theta_sat / theta_0) / gamma if theta_sat > theta_0 else 0
        print(f"  Time to saturation (from noise ~{theta_0}):")
        print(f"    t_sat ~ ln({theta_sat:.4f}/{theta_0}) / {gamma:.4f} = {t_sat:.1f}")

    # Dispersion relation: omega²(q) = q² + m_theta² - sigma_cross * phi2_avg
    print(f"\n  Dispersion relation: omega²(q) = q² + {omega2_q0:.4f}")
    if omega2_q0 < 0:
        q_max = math.sqrt(-omega2_q0)
        print(f"  Unstable band: q < {q_max:.4f}")
        print(f"  Unstable wavelength: lambda > {2*math.pi/q_max:.2f}")
    else:
        print(f"  All modes stable (omega²(q) > 0 for all q)")


def load_diag(path):
    """Load diag TSV."""
    with open(path) as f:
        header = f.readline().strip().split('\t')
        data = []
        for line in f:
            vals = line.strip().split('\t')
            if len(vals) >= len(header):
                data.append([float(v) for v in vals])
    return header, np.array(data)


def compare_to_sim(path, label, A):
    """Compare predictions to simulation diagnostics."""
    header, data = load_diag(path)
    t_col = header.index('t')
    trms_col = header.index('theta_rms')
    et_col = header.index('E_total')
    ep_col = header.index('E_pot')

    t = data[:, t_col]
    trms = data[:, trms_col]
    E_total = data[:, et_col]
    E_pot = data[:, ep_col]

    # Energy drift
    E0 = E_total[0]
    drift = (E_total[-1] - E0) / abs(E0) * 100

    print(f"\n{'='*60}")
    print(f"  SIMULATION RESULTS: {label} (A={A})")
    print(f"{'='*60}")
    print(f"  Time range: {t[0]:.1f} - {t[-1]:.1f}")
    print(f"  Energy drift: {drift:+.3f}%")
    print(f"  theta_rms: {trms[0]:.6f} -> {trms[-1]:.6f}")
    print(f"  E_pot: {E_pot[0]:.1f} -> {E_pot[-1]:.1f}")
    print()

    m2_theta = m_theta ** 2
    phi2_avg = 1.5 * A ** 2
    omega2_q0 = m2_theta - sigma_c * phi2_avg

    if omega2_q0 >= 0:
        # Stable: theta should stay small
        max_trms = np.max(trms)
        print(f"  Prediction: theta stays small (stable vacuum)")
        print(f"  Measured max(theta_rms) = {max_trms:.6f}")
        if max_trms < 0.05:
            print(f"  MATCH: theta remained small ✓")
        else:
            print(f"  MISMATCH: theta grew beyond expected ✗")
    else:
        # Unstable: find growth rate from data
        gamma_pred = math.sqrt(-omega2_q0)
        theta_sat_pred = math.sqrt(-omega2_q0 / lambda_s) if lambda_s > 0 else float('inf')

        # Find growth phase: where theta is growing exponentially
        # Fit log(theta_rms) vs t in the early growth phase
        mask = (trms > 1e-6) & (trms < 0.5 * trms.max())
        if np.sum(mask) > 3:
            t_fit = t[mask]
            log_trms = np.log(trms[mask])
            # Linear fit: log(theta) = log(theta_0) + gamma * t
            coeffs = np.polyfit(t_fit, log_trms, 1)
            gamma_meas = coeffs[0]
            print(f"  Growth rate:")
            print(f"    Predicted: gamma = {gamma_pred:.4f}")
            print(f"    Measured:  gamma = {gamma_meas:.4f}")
            print(f"    Ratio: {gamma_meas/gamma_pred:.3f}")
        else:
            print(f"  [Not enough growth data points for exponential fit]")
            gamma_meas = None

        # Saturation level
        # Take mean of theta_rms in the last 20% of simulation
        late_mask = t > 0.8 * t[-1]
        theta_sat_meas = np.mean(trms[late_mask])
        print(f"  Saturation:")
        print(f"    Predicted: |theta|_sat = {theta_sat_pred:.4f}")
        print(f"    Measured:  theta_rms(late) = {theta_sat_meas:.4f}")


if __name__ == "__main__":
    print("V54 VACUUM PAIR EXTENSION — STABILITY ANALYSIS")
    print("="*60)

    # Theoretical predictions
    stability_report(A_stable, "TEST 1 (Stable)")
    stability_report(A_unstable, "TEST 2 (Unstable)")

    # Compare to simulations if data provided
    if len(sys.argv) >= 2:
        compare_to_sim(sys.argv[1], "Stable", A_stable)
    if len(sys.argv) >= 3:
        compare_to_sim(sys.argv[2], "Unstable", A_unstable)
