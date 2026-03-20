#!/usr/bin/env python3
"""
fit_detailed.py — Detailed fitting with multiple ranges and diagnostics.

Usage: python3 fit_detailed.py data/depletion_t0100.tsv data/depletion_t0200.tsv
"""

import sys
import numpy as np
from scipy.optimize import curve_fit

def yukawa_fixed_m(r, A, m=1.5):
    return A * np.exp(-m * r) / r

def yukawa_free_m(r, A, m):
    return A * np.exp(-m * r) / r

def power_law(r, B, n):
    return B / r**n

def two_component_free(r, A, m, B, n):
    return A * np.exp(-m * r) / r + B / r**n

def r_squared(y_data, y_pred):
    ss_res = np.sum((y_data - y_pred)**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    return 1.0 - ss_res / ss_tot

def load_data(fname):
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
    return np.array(r_all), np.array(rho_all), np.array(drho_all), np.array(counts_all)

def fit_range(r_all, drho_all, counts_all, r_min, r_max, label):
    mask = (r_all >= r_min) & (r_all <= r_max) & (counts_all > 0)
    r = r_all[mask]
    y = drho_all[mask]

    # Skip if not enough points or all same sign issues
    if len(r) < 5:
        print(f"  SKIP: only {len(r)} points")
        return

    # Skip if all values are negative (can't fit power law to negative values)
    pos_frac = np.sum(y > 0) / len(y)
    print(f"  Range [{r_min:.0f}, {r_max:.0f}]: {len(r)} pts, "
          f"mean={np.mean(y):.3e}, positive fraction={pos_frac:.2f}")

    if pos_frac < 0.5:
        print(f"  SKIP: mostly negative (noise region)")
        return

    results = {}

    # Yukawa fixed m=1.5
    try:
        popt, _ = curve_fit(lambda r, A: yukawa_fixed_m(r, A, 1.5), r, y, p0=[1.0])
        R2 = r_squared(y, yukawa_fixed_m(r, popt[0], 1.5))
        results['Yukawa(m=1.5)'] = {'R2': R2, 'params': f'A={popt[0]:.4e}'}
    except:
        results['Yukawa(m=1.5)'] = {'R2': -999, 'params': 'FAILED'}

    # Yukawa free m
    try:
        popt, _ = curve_fit(yukawa_free_m, r, y, p0=[1.0, 0.3],
                            bounds=([-np.inf, 0.001], [np.inf, 10.0]))
        R2 = r_squared(y, yukawa_free_m(r, *popt))
        results['Yukawa(m free)'] = {'R2': R2, 'params': f'A={popt[0]:.4e}, m={popt[1]:.4f}'}
    except:
        results['Yukawa(m free)'] = {'R2': -999, 'params': 'FAILED'}

    # Power law
    try:
        # Log-log initial guess
        pos = y > 0
        if np.sum(pos) > 2:
            log_r = np.log(r[pos])
            log_y = np.log(y[pos])
            slope, intercept = np.polyfit(log_r, log_y, 1)
            n0 = -slope
            B0 = np.exp(intercept)
            popt, _ = curve_fit(power_law, r, y, p0=[B0, max(n0, 0.5)],
                                bounds=([0, 0.1], [np.inf, 10.0]))
            R2 = r_squared(y, power_law(r, *popt))
            results['Power law'] = {'R2': R2, 'params': f'B={popt[0]:.4e}, n={popt[1]:.4f}'}
        else:
            results['Power law'] = {'R2': -999, 'params': 'SKIP (no positive data)'}
    except Exception as e:
        results['Power law'] = {'R2': -999, 'params': f'FAILED: {e}'}

    # Two-component free
    try:
        popt, _ = curve_fit(two_component_free, r, y,
                            p0=[0.1, 1.0, 1.0, 2.0],
                            bounds=([-np.inf, 0.001, -np.inf, 0.1],
                                    [np.inf, 10.0, np.inf, 10.0]),
                            maxfev=10000)
        R2 = r_squared(y, two_component_free(r, *popt))
        results['Two-comp(free)'] = {
            'R2': R2,
            'params': f'A={popt[0]:.4e}, m={popt[1]:.4f}, B={popt[2]:.4e}, n={popt[3]:.4f}'
        }
    except Exception as e:
        results['Two-comp(free)'] = {'R2': -999, 'params': f'FAILED: {e}'}

    # Print results
    print(f"  {'Model':>25s} {'R^2':>10s} {'Parameters':>50s}")
    print(f"  {'-'*88}")
    for name, res in results.items():
        r2_str = f"{res['R2']:.6f}" if res['R2'] > -999 else "FAILED"
        print(f"  {name:>25s} {r2_str:>10s} {res['params']:>50s}")

    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 fit_detailed.py <file1.tsv> [file2.tsv]")
        sys.exit(1)

    for fname in sys.argv[1:]:
        r, rho, drho, counts = load_data(fname)

        print(f"\n{'='*80}")
        print(f"FILE: {fname}")
        print(f"{'='*80}")

        # Multiple fitting ranges
        ranges = [
            (5, 40, "Full range (r=5-40)"),
            (5, 15, "Near-field (r=5-15)"),
            (5, 25, "Inner half (r=5-25)"),
            (10, 30, "Mid-range (r=10-30)"),
            (15, 35, "Outer (r=15-35)"),
            (5, 12, "Core halo (r=5-12)"),
            (20, 40, "Far-field (r=20-40)"),
        ]

        for r_min, r_max, label in ranges:
            print(f"\n--- {label} ---")
            fit_range(r, drho, counts, r_min, r_max, label)

        # Log-log slope analysis (model-independent)
        print(f"\n--- LOG-LOG SLOPE ANALYSIS ---")
        mask = (r > 4) & (r < 45) & (counts > 0) & (drho > 0)
        r_pos = r[mask]
        y_pos = drho[mask]

        if len(r_pos) > 5:
            log_r = np.log10(r_pos)
            log_y = np.log10(y_pos)

            # Piecewise slope
            print(f"  {'r':>8s} {'d(logY)/d(logR)':>16s} {'delta_rho':>14s}")
            for i in range(2, len(r_pos)-2):
                slope = (log_y[i+2] - log_y[i-2]) / (log_r[i+2] - log_r[i-2])
                print(f"  {r_pos[i]:8.2f} {slope:16.4f} {y_pos[i]:14.6e}")

            # Overall log-log fit in different ranges
            for rlo, rhi in [(5,15), (10,25), (5,25)]:
                sub = (r_pos >= rlo) & (r_pos <= rhi)
                if np.sum(sub) > 3:
                    coeffs = np.polyfit(log_r[sub], log_y[sub], 1)
                    print(f"  Log-log linear fit r=[{rlo},{rhi}]: slope={coeffs[0]:.3f} "
                          f"(=> n={-coeffs[0]:.3f})")

    # Cross-time comparison
    if len(sys.argv) >= 3:
        print(f"\n{'='*80}")
        print(f"CROSS-TIME COMPARISON")
        print(f"{'='*80}")

        r1, _, drho1, c1 = load_data(sys.argv[1])
        r2, _, drho2, c2 = load_data(sys.argv[2])

        # Common radii
        common_r = sorted(set(r1) & set(r2))
        print(f"{'r':>8s} {'drho_t100':>14s} {'drho_t200':>14s} {'ratio':>10s}")
        for rv in common_r:
            i1 = np.argmin(np.abs(r1 - rv))
            i2 = np.argmin(np.abs(r2 - rv))
            if c1[i1] > 0 and c2[i2] > 0:
                d1 = drho1[i1]
                d2 = drho2[i2]
                ratio = d2 / d1 if abs(d1) > 1e-10 else float('nan')
                if rv % 2 < 0.3 or rv <= 10:
                    print(f"{rv:8.2f} {d1:14.6e} {d2:14.6e} {ratio:10.3f}")

if __name__ == '__main__':
    main()
