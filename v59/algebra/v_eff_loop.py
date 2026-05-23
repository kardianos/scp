#!/usr/bin/env python3
"""
v59/algebra/v_eff_loop.py

ALGEBRAIC 1-loop V_eff(φ_X) for the Brannen-Z₃ Yukawa structure.

In the v59 framework, the "loop integral" is a FINITE SUM over the
sector subspace of Cl(7)_even (NOT a continuum integral on a hand-imposed
spacetime lattice).  This file builds the infrastructure for that.

For each sector X, the 1-loop V_eff(φ_X) is:
  V_eff^{(X)}(φ) = -(N_X/64π²) · Tr[M_X(φ)² · log(M_X(φ)²/μ_X²)]
where:
  - M_X(φ) is the Brannen-Z₃ 3×3 mass matrix at phase φ
  - Tr is the sum over the flavor space (3 generations)
  - N_X is a sector-specific counting factor (= dim of contributing
    G_2-equivariant subspace in the sector ambient)
  - μ_X is the sector's renormalization scale (dim-density-dependent)

The COEFFICIENT of cos(3φ) in V_eff determines the phase pinning.  At the
minimum: ∂V_eff/∂φ = 0 → cos(3φ) at specific value.

For consistency with empirical φ_X ≈ N_X · α(μ_X):
  N_X = dim G_2 = 14 for d-quark
  N_X = dim Spin(5) = 10 for u-quark
  α(μ_X) at sector dim-density scale

This file:
  1. Computes V_eff(φ) symbolically using sympy.
  2. Extracts the coefficient of cos(3φ).
  3. Identifies the sector counting factor structurally.
  4. Compares to the empirical Brannen phase prediction.
"""
import sympy as sp
import math
import json
import os

# Load Cl(7)_even and embedding data
algebra_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(algebra_dir, "cl7_even_structure.json")) as f:
    struct = json.load(f)
with open(os.path.join(algebra_dir, "brannen_embeddings.json")) as f:
    embed = json.load(f)


# =============================================================================
# Part 1: Symbolic Brannen kernel and its spectrum
# =============================================================================
print("=" * 78)
print("Part 1: Symbolic Brannen kernel spectrum")
print("=" * 78)

a, t, phi = sp.symbols('a t phi', positive=True, real=True)

# Brannen amplitudes s_k = a(1 + 2t cos(2πk/3 + φ))
def s(k):
    return a * (1 + 2*t*sp.cos(2*sp.pi*k/3 + phi))

# Masses m_k = s_k²
m_list = [sp.simplify(s(k)**2) for k in range(3)]

# Power sums
P1 = sp.simplify(sum(s(k) for k in range(3)))  # = 3a
P2 = sp.simplify(sum(m_list))                    # = sum m_k
P3 = sp.simplify(sum(s(k)**3 for k in range(3)))  # contains cos(3φ)
P4 = sp.simplify(sum(s(k)**4 for k in range(3)))  # = sum m_k²

print(f"  P_1(s) = Σ s_k = {P1}")
print(f"  P_2(s) = Σ m_k = {sp.trigsimp(P2)}")
print(f"  P_3(s) = Σ s_k³ = {sp.trigsimp(P3)}")
print(f"  P_4(s) = Σ m_k² = {sp.trigsimp(P4)}")


# =============================================================================
# Part 2: 1-loop V_eff(φ) expansion
# =============================================================================
print()
print("=" * 78)
print("Part 2: V_eff(φ) expansion in cos(3φ)")
print("=" * 78)
print("""
The Coleman-Weinberg 1-loop V_eff for a fermion of mass m(Φ) is:
  V_eff^{1-loop}(Φ) ~ -(N_c/64π²) · m(Φ)^4 · [log(m(Φ)²/μ²) - 3/2]

For Brannen-Z₃ kernel M with eigenvalues m_k(φ):
  Tr V_eff = sum_k m_k² · [log(m_k²) - 3/2]   (per renormalization scheme)

We compute Tr V_eff(φ) and extract its φ-dependence.
""")

# Compute Tr[M² log M²] = Σ m_k log m_k
# Expand around a "reference" φ_0 (taken to be 0 for symmetry-breaking).
# To leading order in t·cos(3φ), the φ-dependent piece is:

# Power-series expansion of the masses around t = 0:
# m_k = a²·(1 + 2t cos θ_k)² = a²·(1 + 4t cos θ_k + 4t² cos² θ_k)
# log m_k = log(a²) + log(1 + 4t cos θ_k + 4t² cos² θ_k)
#        ≈ log(a²) + 4t cos θ_k + 4t² cos² θ_k - 8t² cos² θ_k + ...
#        = log(a²) + 4t cos θ_k - 4t² cos² θ_k + ...

# Σ m_k log m_k ≈ a²·(1 + 2t² + 2t³ cos(3φ) terms higher up) × log

# Numerical extraction: for t² = 1/2 (lepton), compute V_eff numerically as a function of φ
# and extract the cos(3φ) Fourier coefficient.

a_num = 1.0  # normalize a = 1 for symbolic clarity
mu_sq = 1.0  # renormalization scale

import numpy as np

def V_eff_lepton(phi_val, t_val):
    """1-loop V_eff for lepton sector as function of φ (numerical), with t fixed.
    Sums m_k log m_k (up to overall normalization)."""
    s_vals = [a_num * (1 + 2*t_val*math.cos(phi_val + 2*math.pi*k/3)) for k in range(3)]
    m_vals = [s_val**2 for s_val in s_vals]
    return sum(m * (math.log(abs(m)/mu_sq) - 3/2) for m in m_vals)

# Compute V_eff as a function of φ for lepton sector
print("\nNumerical V_eff(φ) for lepton sector (t² = 1/2):")
phi_grid = np.linspace(0, 2*np.pi, 1000)
V_grid = np.array([V_eff_lepton(p, math.sqrt(1/2)) for p in phi_grid])

# Fourier decomposition: V(φ) = c_0 + c_3·cos(3φ) + c_6·cos(6φ) + ...
def fourier_cos_coeff(V_vals, n, phi_vals):
    """Extract coefficient of cos(n·φ) in V(φ)."""
    dphi = phi_vals[1] - phi_vals[0]
    return (1/np.pi) * np.sum(V_vals * np.cos(n * phi_vals)) * dphi

c_0 = (1/(2*np.pi)) * np.sum(V_grid) * (phi_grid[1] - phi_grid[0])
c_3 = fourier_cos_coeff(V_grid, 3, phi_grid)
c_6 = fourier_cos_coeff(V_grid, 6, phi_grid)
c_9 = fourier_cos_coeff(V_grid, 9, phi_grid)
c_12 = fourier_cos_coeff(V_grid, 12, phi_grid)

print(f"  V(φ) ≈ {c_0:.6f} + {c_3:.6f}·cos(3φ) + {c_6:.6f}·cos(6φ) + ...")
print(f"  cos(3φ) coefficient: c_3 = {c_3:.6f}")
print(f"  cos(6φ) coefficient: c_6 = {c_6:.6f}")
print(f"  cos(9φ) coefficient: c_9 = {c_9:.6f}")

# Critical points: ∂V/∂φ = 0
# For V = c_3 cos(3φ) + c_6 cos(6φ) + ...
# ∂V/∂φ = -3 c_3 sin(3φ) - 6 c_6 sin(6φ) - ...
# = sin(3φ)·(-3 c_3 - 12 c_6 cos(3φ) - ...)

# So critical points at sin(3φ) = 0 (φ = 0, π/3, 2π/3, π, 4π/3, 5π/3)
# OR at cos(3φ) = -3 c_3 / (12 c_6) = -c_3/(4 c_6)

cos_3phi_critical = -c_3 / (4 * c_6)
print(f"\n  Critical cos(3φ) from V'(φ)=0:  cos(3φ) = -c_3/(4c_6) = {cos_3phi_critical:.5f}")
print(f"    → φ_critical = arccos(cos)/3")
if abs(cos_3phi_critical) <= 1:
    phi_crit = math.acos(cos_3phi_critical) / 3
    print(f"    φ_critical ≈ {phi_crit:.5f} rad")
    print(f"    vs empirical φ_l = 2/9 = {2/9:.5f}")

# =============================================================================
# Part 3: Sector-specific counting via dim-density
# =============================================================================
print()
print("=" * 78)
print("Part 3: Sector counting factor N_X (algebraic interpretation)")
print("=" * 78)
print("""
The empirical observation:
  φ_d ≈ +14·α(M_Z)  →  N_d = 14 = dim G_2
  φ_u ≈ -10·α(0)   →  N_u = 10 = dim Spin(5) (or 2·Killing-index)

In the algebraic 1-loop calculation, N_X arises from SUMMING over the
sector-specific contributing states in Cl(7)_even.

For d-quark (F = Λ⁴ ambient):
  Loop integrand sums over Λ⁴ states.  The G_2-EQUIVARIANT projection
  picks out 14 states (the dim of the G_2 adjoint rep acting in F).

For u-quark (L⊕F ambient):
  Loop integrand sums over BOTH L and F.  The "10" = dim Spin(5) might
  correspond to the BICOMPLEX structure of L⊕F under the relevant
  Spin(7) sub-action.
""")

# Numerical comparison
alpha_0 = 1.0/137.036
alpha_MZ = 1.0/127.952
phi_d_pred = 14 * alpha_MZ
phi_u_pred = -10 * alpha_0
phi_d_emp = 0.10859
phi_u_emp = -0.07251
print(f"Numerical predictions:")
print(f"  φ_d = 14·α(M_Z) = 14/127.952 = {phi_d_pred:.6f}  (empirical {phi_d_emp:+.5f}, gap {abs(phi_d_pred-phi_d_emp)/phi_d_emp*100:.3f}%)")
print(f"  φ_u = -10·α(0) = -10/137.036 = {phi_u_pred:.6f}  (empirical {phi_u_emp:+.5f}, gap {abs(phi_u_pred-phi_u_emp)/abs(phi_u_emp)*100:.3f}%)")


# =============================================================================
# Part 4: Symbolic cos(3φ) coefficient at order t³
# =============================================================================
print()
print("=" * 78)
print("Part 4: Symbolic structure — cos(3φ) coefficient at leading order in t")
print("=" * 78)

# Expand m_k log m_k around t = 0 to extract cos(3φ) coefficient
# We need higher-order term where cos(3φ) appears.

# log m_k = log(s_k²) = 2 log|s_k|
# For s_k = a(1 + 2t cos θ_k), and small t:
# log s_k ≈ log a + 2t cos θ_k - 2t² cos² θ_k + (8/3)t³ cos³ θ_k - ...

# m_k log m_k = s_k² · 2 log|s_k|
# Σ over k: Σ s_k² · 2 log s_k = 2 a² · Σ (1 + 2t cos θ_k)² · log[a(1 + 2t cos θ_k)]
#         = 2 a² log a · Σ (1 + 2t cos θ_k)² + 2 a² Σ (1+2t cos θ_k)² log(1 + 2t cos θ_k)

# The first term gives the φ-INDEPENDENT P_2(s)/a · 2 log a
# The second term contains the φ-dependence.

# Expand log(1 + x) = x - x²/2 + x³/3 - x⁴/4 + ...
# with x = 2t cos θ_k.

# (1 + x)² log(1 + x) = (1 + 2x + x²)(x - x²/2 + x³/3 - x⁴/4 + ...)
# = x - x²/2 + x³/3 - x⁴/4 + ...
# + 2x² - x³ + 2x⁴/3 - x⁵/2 + ...
# + x³ - x⁴/2 + x⁵/3 - x⁶/4 + ...
# = x + 3x²/2 + (1/3 - 1 + 1)x³ + (...) x⁴ + ...
# = x + 3x²/2 + x³/3 + ...

# So Σ (1+x_k)² log(1+x_k) = Σ x_k + (3/2) Σ x_k² + (1/3) Σ x_k³ + ...
# with x_k = 2t cos θ_k.

# Σ x_k = 2t Σ cos θ_k = 0 (φ-indep)
# Σ x_k² = 4t² Σ cos² θ_k = 4t² · 3/2 = 6t² (φ-indep!)
# Σ x_k³ = 8t³ Σ cos³ θ_k = 8t³ · (3/4) cos(3φ) = 6t³ cos(3φ)  ← cos(3φ)!
# Σ x_k⁴ = 16t⁴ · (9/8) = 18t⁴  (φ-indep)
# Σ x_k⁵ = 32t⁵ · (5/8) cos(3φ) = 20 t⁵ cos(3φ)

# So the leading cos(3φ) piece comes at order t³:
# Σ (1+x_k)² log(1+x_k) [order t³] = (1/3) · 6 t³ cos(3φ) = 2 t³ cos(3φ)

# Therefore: 1-loop V_eff ∝ 2 a² · 2 t³ cos(3φ) = 4 a² t³ cos(3φ) + higher orders

# Combined with the (N_c/64π²) factor and the sector counting N_X:
# V_eff(φ) = -(N_X/64π²) · 4 a² t³ cos(3φ) + ...

# At the minimum, ∂V/∂φ = (N_X/64π²) · 12 a² t³ sin(3φ) = 0
# → sin(3φ) = 0 → φ = 0, π/3, 2π/3 (commensurate)

# But we observe φ_X = small non-zero values.  This means the LEADING t³ contribution
# alone doesn't pin φ_X.  Higher orders (t⁵, t⁷, ...) bring in additional cos(6φ),
# cos(9φ) terms that break the commensurate degeneracy.

print("Expansion of Σ m_k log m_k in powers of t:")
print()
print("  Σ m_k log m_k = (φ-indep) + 4·a²·t³·cos(3φ) + (higher orders in t)")
print()
print("This is the LEADING cos(3φ) coefficient at order t³.")
print()
print("For sector X with t² = 1 - 14/D_X:")
for sec, t_sq, D in [("lepton",  1/2, 28),
                     ("d-quark", 3/5, 35),
                     ("u-quark", 7/9, 63)]:
    t_val = math.sqrt(t_sq)
    leading_coeff = 4 * t_val**3
    print(f"  {sec:8s}: 4·t³ = 4·({t_val:.4f})³ = {leading_coeff:.5f}")

print()
print("=" * 78)
print("Part 5: Putting it together — predicted φ_X from V_eff")
print("=" * 78)
print("""
At 1-loop, with V_eff ∝ -N_X · α · cos(3φ) · (4 t³) + higher:

If we ASSUME a small additional symmetry-breaking term ε·cos(6φ) (from
2-loop or higher), then:
  V_eff(φ) ≈ -α·N_X · 4t³·cos(3φ) + ε·cos(6φ)
  ∂V/∂φ = +α·N_X · 12 t³ · sin(3φ) - 6 ε · sin(6φ)
  At minimum: sin(3φ)·[12·α·N_X·t³ - 12·ε·cos(3φ)] = 0
  → cos(3φ) = α·N_X·t³/ε

  φ_X = (1/3)·arccos(α·N_X·t³/ε)

For φ_X SMALL (which we observe), we need cos(3φ_X) ≈ 1 → arccos ≈ 0.
This means α·N_X·t³/ε ≈ 1.

Expansion around cos(3φ) = 1:
  3φ ≈ √(2·(1 - α·N_X·t³/ε))
  φ ≈ (1/3)·√(2 - 2α·N_X·t³/ε)

For SMALL φ, this requires ε just slightly above α·N_X·t³.

This gives φ_X as a small NON-ZERO value, with magnitude controlled by
the ratio of N_X·α·t³ to ε.  The empirical φ_X = N_X·α (linear, not sqrt)
suggests a DIFFERENT loop structure — perhaps:
  V_eff(φ) ∝ -α·N_X · 4t³ · cos(3φ) + (sector-specific tilting)
where the tilting comes from cross-sector mixing terms or from L-vs-F
G_2-equivariant projection.

CONCLUSION OF THIS SETUP:
  - Cl(7)_even algebra constructed (cl7_even.py)
  - Sector embeddings defined (brannen_kernel.py)
  - V_eff symbolic + numerical extraction (this file)
  - Leading cos(3φ) coefficient at order t³: 4·a²·t³ (universal)
  - Sector counting N_X enters as overall factor — needs explicit G_2/Spin(5)
    representation theory calculation to derive N_X structurally.

NEXT STEP: explicit calculation of the G_2-equivariant projection of the F
ambient (for d-quark) and the Spin(5)-equivariant projection for u-quark.
This requires building the action of G_2 ⊂ Spin(7) on Λ⁴(R⁷) and the
Spin(5) action on the L⊕F structure.
""")


# Save results
out_path = os.path.join(algebra_dir, "v_eff_loop_results.json")
results = {
    "leading_cos_3phi_coefficient_at_t3": 4.0,
    "sector_t_cubed": {
        "lepton":   math.sqrt(0.5)**3,
        "d-quark":  math.sqrt(0.6)**3,
        "u-quark":  math.sqrt(7/9)**3,
    },
    "predicted_phi_X_from_alpha_N": {
        "phi_d_14_alpha_MZ": 14*alpha_MZ,
        "phi_u_neg_10_alpha_0": -10*alpha_0,
    },
    "empirical_phi_X": {
        "phi_l": 2/9,
        "phi_d": 0.10859,
        "phi_u": -0.07251,
    },
    "fourier_lepton_V_eff": {
        "c_0": c_0,
        "c_3": c_3,
        "c_6": c_6,
        "c_9": c_9,
        "cos_3phi_critical": cos_3phi_critical,
    },
    "notes": "1-loop V_eff infrastructure set up.  Sector counting N_X (= dim G_2 "
             "for d-quark, dim Spin(5) for u-quark) requires explicit "
             "G_2-equivariant projection of Λ⁴(R⁷) — next computational step.",
}
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"Results saved to {out_path}")
print()
print("INFRASTRUCTURE READY.  Next:")
print("  1. Build G_2 action on Λ⁴(R⁷) → identify 14-dim invariant subspace")
print("  2. Compute the 1-loop V_eff restricted to G_2-equivariant content")
print("  3. Verify the cos(3φ_X) coefficient = N_X · α(μ_X) structurally")
