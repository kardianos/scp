#!/usr/bin/env python3
"""
fit_depletion.py — Fit radial depletion profile to Yukawa, power-law, and two-component models.

Usage: python3 fit_depletion.py data/depletion_t0200.tsv
"""

import sys
import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.stats import pearsonr

def yukawa_fixed_m(r, A, m=1.5):
    """Yukawa with fixed m=1.5"""
    return A * np.exp(-m * r) / r

def yukawa_free_m(r, A, m):
    """Yukawa with free m"""
    return A * np.exp(-m * r) / r

def power_law(r, B, n):
    """Power law: B / r^n"""
    return B / r**n

def two_component(r, A, B, n, m=1.5):
    """Two-component: Yukawa (m fixed) + power law"""
    return A * np.exp(-m * r) / r + B / r**n

def two_component_free(r, A, m, B, n):
    """Two-component with free m"""
    return A * np.exp(-m * r) / r + B / r**n

def r_squared(y_data, y_pred):
    """Coefficient of determination"""
    ss_res = np.sum((y_data - y_pred)**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    return 1.0 - ss_res / ss_tot

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 fit_depletion.py <depletion.tsv>")
        sys.exit(1)

    fname = sys.argv[1]

    # Read data
    r_all, rho_all, drho_all, counts_all = [], [], [], []
    with open(fname) as f:
        for line in f:
            if line.startswith('#') or line.startswith('r'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                r_all.append(float(parts[0]))
                rho_all.append(float(parts[1]))
                drho_all.append(float(parts[2]))
                counts_all.append(int(parts[3]))

    r_all = np.array(r_all)
    rho_all = np.array(rho_all)
    drho_all = np.array(drho_all)
    counts_all = np.array(counts_all)

    # Print full profile for RESULTS.md
    print("=" * 70)
    print("FULL RADIAL DEPLETION PROFILE")
    print("=" * 70)
    print(f"{'r':>8s} {'rho(r)':>14s} {'delta_rho':>14s} {'counts':>10s}")
    print("-" * 50)
    for i in range(len(r_all)):
        if i % 2 == 0 or r_all[i] <= 10:
            print(f"{r_all[i]:8.2f} {rho_all[i]:14.6e} {drho_all[i]:14.6e} {counts_all[i]:10d}")

    # Fitting range: r=5 to r=40 (avoid core and boundary)
    mask = (r_all >= 5.0) & (r_all <= 40.0) & (counts_all > 0)
    r = r_all[mask]
    drho = drho_all[mask]

    # Check sign of depletion
    print(f"\n\nFitting range: r = [{r[0]:.1f}, {r[-1]:.1f}], {len(r)} points")
    print(f"Mean delta_rho in fit range: {np.mean(drho):.6e}")
    print(f"Max |delta_rho| in fit range: {np.max(np.abs(drho)):.6e}")

    # Use absolute values for fitting (depletion is typically negative)
    # The sign tells us attraction vs repulsion
    sign = np.sign(np.mean(drho))
    y = drho  # keep sign for now
    abs_y = np.abs(drho)

    print(f"Sign of depletion: {'negative (attractive)' if sign < 0 else 'positive (excess)'}")

    # ============================================================
    # Model 1: Yukawa with fixed m=1.5
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL 1: Yukawa (fixed m=1.5)")
    print("  delta_rho = A * exp(-1.5*r) / r")
    print("=" * 70)

    try:
        popt, pcov = curve_fit(lambda r, A: yukawa_fixed_m(r, A, 1.5), r, y, p0=[sign * 1.0])
        y_pred = yukawa_fixed_m(r, popt[0], 1.5)
        R2 = r_squared(y, y_pred)
        print(f"  A = {popt[0]:.6e}")
        print(f"  R^2 = {R2:.6f}")
        print(f"  Residual range: [{np.min(y - y_pred):.4e}, {np.max(y - y_pred):.4e}]")
    except Exception as e:
        print(f"  FAILED: {e}")
        R2 = -999

    # ============================================================
    # Model 1b: Yukawa with free m
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL 1b: Yukawa (free m)")
    print("  delta_rho = A * exp(-m*r) / r")
    print("=" * 70)

    try:
        popt_yf, pcov_yf = curve_fit(yukawa_free_m, r, y, p0=[sign * 1.0, 0.5],
                                      bounds=([-np.inf, 0.001], [np.inf, 10.0]))
        y_pred_yf = yukawa_free_m(r, *popt_yf)
        R2_yf = r_squared(y, y_pred_yf)
        print(f"  A = {popt_yf[0]:.6e}")
        print(f"  m = {popt_yf[1]:.6f}  (Yukawa range = {1.0/popt_yf[1]:.3f})")
        print(f"  R^2 = {R2_yf:.6f}")
        print(f"  Residual range: [{np.min(y - y_pred_yf):.4e}, {np.max(y - y_pred_yf):.4e}]")

        if popt_yf[1] < 0.1:
            print(f"  >>> NOTE: m={popt_yf[1]:.4f} << 1.5 — this is effectively power-law, NOT Yukawa")
    except Exception as e:
        print(f"  FAILED: {e}")
        R2_yf = -999
        popt_yf = [0, 0]

    # ============================================================
    # Model 2: Power law
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL 2: Power law")
    print("  delta_rho = B / r^n")
    print("=" * 70)

    try:
        # For log-log fit: use |drho| > 0
        pos_mask = abs_y > 0
        if np.sum(pos_mask) > 2:
            # Initial guess from log-log linear fit
            log_r = np.log(r[pos_mask])
            log_y = np.log(abs_y[pos_mask])
            slope, intercept = np.polyfit(log_r, log_y, 1)
            n0 = -slope
            B0 = sign * np.exp(intercept)

            popt_pl, pcov_pl = curve_fit(power_law, r, y, p0=[B0, n0],
                                          bounds=([-np.inf, 0.1], [np.inf, 10.0]) if sign > 0
                                          else ([-np.inf, 0.1], [np.inf, 10.0]))
            y_pred_pl = power_law(r, *popt_pl)
            R2_pl = r_squared(y, y_pred_pl)

            print(f"  B = {popt_pl[0]:.6e}")
            print(f"  n = {popt_pl[1]:.4f}")
            print(f"  R^2 = {R2_pl:.6f}")
            print(f"  Residual range: [{np.min(y - y_pred_pl):.4e}, {np.max(y - y_pred_pl):.4e}]")

            # Physical interpretation
            if abs(popt_pl[1] - 2.0) < 0.3:
                print(f"  >>> CONSISTENT with 1/r^2 (gravitational-like)")
            elif abs(popt_pl[1] - 1.0) < 0.3:
                print(f"  >>> CONSISTENT with 1/r (Coulomb-like potential)")
            elif abs(popt_pl[1] - 4.0) < 0.5:
                print(f"  >>> CONSISTENT with 1/r^4 (dipole-like)")
            elif abs(popt_pl[1] - 6.0) < 0.5:
                print(f"  >>> CONSISTENT with 1/r^6 (van der Waals-like)")
        else:
            print("  SKIPPED: not enough positive data points")
            R2_pl = -999
            popt_pl = [0, 0]
    except Exception as e:
        print(f"  FAILED: {e}")
        R2_pl = -999
        popt_pl = [0, 0]

    # ============================================================
    # Model 3: Two-component (Yukawa + power law), m fixed
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL 3: Two-component (fixed m=1.5)")
    print("  delta_rho = A * exp(-1.5*r)/r + B / r^n")
    print("=" * 70)

    try:
        def two_comp_fixed(r, A, B, n):
            return two_component(r, A, B, n, m=1.5)

        popt_2c, pcov_2c = curve_fit(two_comp_fixed, r, y,
                                      p0=[sign * 0.1, sign * 1.0, 2.0],
                                      bounds=([-np.inf, -np.inf, 0.1], [np.inf, np.inf, 10.0]),
                                      maxfev=10000)
        y_pred_2c = two_comp_fixed(r, *popt_2c)
        R2_2c = r_squared(y, y_pred_2c)

        print(f"  A = {popt_2c[0]:.6e}  (Yukawa amplitude)")
        print(f"  B = {popt_2c[1]:.6e}  (power-law amplitude)")
        print(f"  n = {popt_2c[2]:.4f}  (power-law exponent)")
        print(f"  R^2 = {R2_2c:.6f}")

        # Evaluate contributions at r=10
        yuk_10 = yukawa_fixed_m(10, popt_2c[0], 1.5)
        pl_10  = power_law(10, popt_2c[1], popt_2c[2])
        print(f"  At r=10: Yukawa={yuk_10:.4e}, Power-law={pl_10:.4e}, ratio={abs(yuk_10/(pl_10+1e-50)):.2e}")

        yuk_30 = yukawa_fixed_m(30, popt_2c[0], 1.5)
        pl_30  = power_law(30, popt_2c[1], popt_2c[2])
        print(f"  At r=30: Yukawa={yuk_30:.4e}, Power-law={pl_30:.4e}, ratio={abs(yuk_30/(pl_30+1e-50)):.2e}")
    except Exception as e:
        print(f"  FAILED: {e}")
        R2_2c = -999

    # ============================================================
    # Model 3b: Two-component, m free
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL 3b: Two-component (free m)")
    print("  delta_rho = A * exp(-m*r)/r + B / r^n")
    print("=" * 70)

    try:
        popt_2cf, pcov_2cf = curve_fit(two_component_free, r, y,
                                        p0=[sign * 0.1, 1.0, sign * 1.0, 2.0],
                                        bounds=([-np.inf, 0.001, -np.inf, 0.1],
                                                [np.inf, 10.0, np.inf, 10.0]),
                                        maxfev=10000)
        y_pred_2cf = two_component_free(r, *popt_2cf)
        R2_2cf = r_squared(y, y_pred_2cf)

        print(f"  A = {popt_2cf[0]:.6e}  (Yukawa amplitude)")
        print(f"  m = {popt_2cf[1]:.6f}  (Yukawa mass)")
        print(f"  B = {popt_2cf[2]:.6e}  (power-law amplitude)")
        print(f"  n = {popt_2cf[3]:.4f}  (power-law exponent)")
        print(f"  R^2 = {R2_2cf:.6f}")
    except Exception as e:
        print(f"  FAILED: {e}")
        R2_2cf = -999

    # ============================================================
    # Log-log slope analysis (model-independent)
    # ============================================================
    print("\n" + "=" * 70)
    print("MODEL-INDEPENDENT: Local log-log slope")
    print("  d(log|delta_rho|)/d(log r)")
    print("=" * 70)

    nz = abs_y > 0
    if np.sum(nz) > 3:
        lr = np.log(r[nz])
        ly = np.log(abs_y[nz])

        # Sliding window slope (5-point)
        print(f"{'r_center':>10s} {'local_slope':>12s} {'|delta_rho|':>14s}")
        for i in range(2, len(lr) - 2):
            slope = (ly[i+2] - ly[i-2]) / (lr[i+2] - lr[i-2])
            print(f"{r[nz][i]:10.2f} {slope:12.4f} {abs_y[nz][i]:14.6e}")

    # ============================================================
    # Summary
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY OF FITS (r = 5 to 40)")
    print("=" * 70)
    fits = [
        ("Yukawa (m=1.5 fixed)", R2, "A*exp(-1.5r)/r"),
        ("Yukawa (m free)", R2_yf, f"A*exp(-{popt_yf[1]:.4f}r)/r" if R2_yf > -999 else "FAILED"),
        ("Power law", R2_pl, f"B/r^{popt_pl[1]:.4f}" if R2_pl > -999 else "FAILED"),
    ]
    if R2_2c > -999:
        fits.append(("Two-component (m=1.5)", R2_2c, f"Yukawa + B/r^{popt_2c[2]:.4f}"))
    if R2_2cf > -999:
        fits.append(("Two-component (m free)", R2_2cf, f"Yukawa(m={popt_2cf[1]:.4f}) + B/r^{popt_2cf[3]:.4f}"))

    print(f"{'Model':>30s} {'R^2':>10s} {'Form':>30s}")
    print("-" * 72)
    for name, r2, form in fits:
        print(f"{name:>30s} {r2:10.6f} {form:>30s}")

    # Determine winner
    best_r2 = max([r2 for _, r2, _ in fits])
    best_name = [name for name, r2, _ in fits if r2 == best_r2][0]
    print(f"\nBest fit: {best_name} (R^2 = {best_r2:.6f})")

    if "Power" in best_name or (R2_yf > -999 and popt_yf[1] < 0.1):
        print("\n>>> CONCLUSION: Depletion profile is POWER-LAW (not Yukawa)")
        if R2_pl > -999:
            n_val = popt_pl[1]
            print(f">>> Exponent n = {n_val:.3f}")
            if abs(n_val - 2.0) < 0.5:
                print(f">>> This is CONSISTENT with gravitational 1/r^2 scaling")
    elif "Yukawa" in best_name and "free" in best_name and popt_yf[1] < 0.3:
        print("\n>>> CONCLUSION: Depletion is effectively POWER-LAW (Yukawa m -> 0)")
    else:
        print("\n>>> CONCLUSION: Depletion is YUKAWA (massive, exponential decay)")

if __name__ == '__main__':
    main()
