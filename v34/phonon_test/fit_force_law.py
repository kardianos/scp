"""
Phonon Hypothesis Test: Fit V33-C1 force data to three models:
  1. Pure Yukawa: F = A × e^{-mD} / D²
  2. Pure power law: F = B / D^n
  3. Two-component: F = A × e^{-mD} / D² + B / D^n

If the two-component model fits significantly better AND the power-law
component is dominant at large D, the phonon carries the long-range force.
"""

import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# V33-C1 data: D_init, dD/dt (velocity = proxy for force)
# From V33_RESULTS.md — using dD/dt as force proxy (F ∝ dD/dt for equal-mass braids)
# Only using attractive regime (dD/dt < 0, D ≥ 12)
data = [
    # (D, -dD/dt)  — positive = attractive force
    (12, 0.04217),
    (15, 0.04819),
    (18, 0.00793),
    (20, 0.00623),
    (25, 0.00451),
    (30, 0.01266),
    (40, 0.00673),
    (80, 0.00529),
]

D = np.array([d[0] for d in data])
F = np.array([d[1] for d in data])

m = 1.5  # mass parameter

print("=" * 70)
print("PHONON HYPOTHESIS TEST: Force Law Decomposition")
print("=" * 70)
print(f"\nData points: {len(data)} (attractive regime, D=12-80)")
print(f"Mass parameter: m = {m}")
print(f"Yukawa range: 1/m = {1/m:.3f}")
print(f"Yukawa at D=15: e^(-mD) = {np.exp(-m*15):.2e}")
print(f"Yukawa at D=80: e^(-mD) = {np.exp(-m*80):.2e}")
print()

# Model 1: Pure Yukawa
def yukawa(D, A):
    return A * np.exp(-m * D) / D**2

# Model 2: Pure power law
def power_law(D, B, n):
    return B / D**n

# Model 3: Two-component (Yukawa + power law)
def two_comp(D, A, B, n):
    return A * np.exp(-m * D) / D**2 + B / D**n

# Model 4: Pure power law (single parameter, n fixed)
def power_law_fixed(D, B):
    return B / D**2

print("--- Model 1: Pure Yukawa F = A × e^{-mD} / D² ---")
try:
    popt, pcov = curve_fit(yukawa, D, F, p0=[1e10])
    F_pred = yukawa(D, *popt)
    ss_res = np.sum((F - F_pred)**2)
    ss_tot = np.sum((F - F.mean())**2)
    r2 = 1 - ss_res / ss_tot
    print(f"  A = {popt[0]:.4e}")
    print(f"  R² = {r2:.6f}")
    print(f"  Residuals: {ss_res:.6e}")
    print(f"  At D=80: F_pred = {yukawa(80, *popt):.2e} (actual: {F[-1]:.2e})")
except Exception as e:
    print(f"  FAILED: {e}")
    r2_yuk = -1

print()
print("--- Model 2: Pure Power Law F = B / D^n ---")
try:
    popt2, pcov2 = curve_fit(power_law, D, F, p0=[10, 2])
    F_pred2 = power_law(D, *popt2)
    ss_res2 = np.sum((F - F_pred2)**2)
    r2_pl = 1 - ss_res2 / ss_tot
    print(f"  B = {popt2[0]:.4f}")
    print(f"  n = {popt2[1]:.4f}")
    print(f"  R² = {r2_pl:.6f}")
    print(f"  At D=80: F_pred = {power_law(80, *popt2):.2e} (actual: {F[-1]:.2e})")
except Exception as e:
    print(f"  FAILED: {e}")
    r2_pl = -1

print()
print("--- Model 3: Two-Component F = A×e^{-mD}/D² + B/D^n ---")
try:
    popt3, pcov3 = curve_fit(two_comp, D, F, p0=[1e5, 1, 1.5], maxfev=10000)
    F_pred3 = two_comp(D, *popt3)
    ss_res3 = np.sum((F - F_pred3)**2)
    r2_tc = 1 - ss_res3 / ss_tot
    A, B, n = popt3
    print(f"  Yukawa:    A = {A:.4e}")
    print(f"  Power law: B = {B:.4f}, n = {n:.4f}")
    print(f"  R² = {r2_tc:.6f}")
    
    # Decompose at each D
    print(f"\n  Decomposition:")
    print(f"  {'D':>5s}  {'F_total':>10s}  {'F_yukawa':>10s}  {'F_power':>10s}  {'%_power':>8s}")
    for d in [12, 15, 20, 30, 50, 80]:
        fy = A * np.exp(-m * d) / d**2
        fp = B / d**n
        ft = fy + fp
        pct = 100 * fp / ft if ft > 0 else 0
        print(f"  {d:5d}  {ft:10.4e}  {fy:10.4e}  {fp:10.4e}  {pct:7.1f}%")
except Exception as e:
    print(f"  FAILED: {e}")
    r2_tc = -1

print()
print("--- Model 4: Pure 1/D² (Newton) ---")
try:
    popt4, _ = curve_fit(power_law_fixed, D, F, p0=[1])
    F_pred4 = power_law_fixed(D, *popt4)
    ss_res4 = np.sum((F - F_pred4)**2)
    r2_n = 1 - ss_res4 / ss_tot
    print(f"  B = {popt4[0]:.4f}")
    print(f"  R² = {r2_n:.6f}")
except Exception as e:
    print(f"  FAILED: {e}")

print()
print("=" * 70)
print("COMPARISON")
print("=" * 70)
print(f"  Pure Yukawa:    R² = {r2:.6f}" if r2 > -1 else "  Pure Yukawa: FAILED")
print(f"  Pure power law: R² = {r2_pl:.6f} (n={popt2[1]:.2f})" if r2_pl > -1 else "  Pure power law: FAILED")
print(f"  Two-component:  R² = {r2_tc:.6f}" if r2_tc > -1 else "  Two-component: FAILED")
print(f"  Pure 1/D²:      R² = {r2_n:.6f}" if r2_n > -1 else "  Pure Newton: FAILED")
print()
if r2_tc > max(r2, r2_pl, r2_n):
    print("  TWO-COMPONENT WINS → phonon hypothesis supported")
elif r2_pl > max(r2, r2_tc, r2_n):
    print("  PURE POWER LAW WINS → force is NOT Yukawa at all")
elif r2 > max(r2_pl, r2_tc, r2_n):
    print("  PURE YUKAWA WINS → no phonon component")
