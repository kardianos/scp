#!/usr/bin/env python3
"""
04_lagrangian_modes.py

Full v59 strain Lagrangian — mode decomposition and natural-scale analysis.

The user's caution (this session): the strain Lagrangian should be written
*in concert* with the rest of the v59 structure, not in isolation.  Here we
decompose the natural quaternion-valued field

    ξ(x) ∈ ℍ ≅ ℝ⁴,         constraint  |ξ|² = ρ₀²  (the v59 S³)

into its four independent modes at every spacetime point, compute the kinetic
term in each, and identify the candidate physical sectors and their natural
coupling scales.

Modes (at the Brannen equilibrium point ξ = ρ₀(cos δ + sin δ · n̂), where
n̂ ∈ S² is a unit pure-imaginary quaternion):

  (R) RADIAL          δρ := |ξ| − ρ₀  →  off-S³, gives δQ from the constraint
                          (changes Koide ratio away from 2/3 if excited).
                      Candidate for G_e (cross-sector ratio depends on |ξ|).

  (A) ACTIVE phase    ψ := arctan(|Im ξ|/Re ξ) − δ_B
                          →  on S³, in the (1,i)-plane, changes Brannen phase
                          → flavor-dependent lepton mass shifts (ruled out at
                          gravity strength: experiment 01).
                      Candidate for α via U(1) phase / Hopf fibration.

  (S) SILENT direction n̂ ∈ S² (2-dim sphere of imaginary-direction rotations).
                      EXACTLY invariant for the Brannen spectrum (Lean theorem
                      `SilentDirection.silent_pair`; numerical 4×10⁻¹⁵).
                      Candidate for SU(2)_L (Furey identification, now derived
                      kinematically).

The kinetic term decomposes as:

    L_kin = (1/2) (∂μ ξ)·(∂μ ξ)
          = (1/2) (∂ρ)²                       (R)
            + (ρ²/2) (∂ψ)²                    (A)
            + (ρ²/2) sin²(ψ) (∂n̂)²             (S, sigma-model on S²)

At ρ = ρ₀ = 1/√2 and ψ = δ_B = 2/9 (the Brannen equilibrium):

    f_S² := ρ₀ · sin(δ_B)        ← silent-direction decay constant
    f_A  := ρ₀                   ← active-direction "decay constant" (× modulator)
    f_R  := 1                    ← radial direction (unconstrained by ρ₀)

The relations to empirical couplings are read off from the gauged sigma model
on each sector.  We compute the natural scales here and compare to:

    α     = 1/137.036              (U(1)_em fine-structure)
    α_W   = g₂²/(4π) ≈ 1/29.78     (SU(2)_L weak coupling)
    G_e   ≈ 1.75 × 10⁻⁴⁵           (dimensionless gravitational coupling)
"""

import numpy as np

# ----------------------------------------------------------------------
# v59 dimensionless parameters
# ----------------------------------------------------------------------
rho_0 = 1 / np.sqrt(2)         # the S³ constraint radius |ξ| = 1/√2
delta_B = 2.0 / 9.0            # the Brannen phase, 2/9 rad

# Equilibrium quaternion components at the physical Brannen point
xi_re = rho_0 * np.cos(delta_B)              # Re ξ
xi_im_mag = rho_0 * np.sin(delta_B)          # |Im ξ|
# Sanity
assert abs(xi_re**2 + xi_im_mag**2 - rho_0**2) < 1e-14

print("=" * 72)
print("v59 Brannen-equilibrium values (dimensionless, in units of ρ₀ = 1/√2)")
print("=" * 72)
print(f"  ρ₀         = {rho_0:.6f}  (S³ constraint radius)")
print(f"  δ_B        = {delta_B:.6f}  (Brannen phase, in rad, = 2/9)")
print(f"  Re ξ       = {xi_re:.6f}")
print(f"  |Im ξ|     = {xi_im_mag:.6f}")
print(f"  |ξ|²       = {xi_re**2 + xi_im_mag**2:.6f}  (= 1/2 exactly)")

# ----------------------------------------------------------------------
# Mode-by-mode "decay constants" (kinetic-term coefficients)
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Mode kinetic-term coefficients (decay constants f_X²)")
print("=" * 72)

f_R_sq = 1.0                                # radial: just (∂ρ)², coefficient 1
f_A_sq = rho_0**2                           # active: ρ₀² (∂ψ)²
f_S_sq = rho_0**2 * np.sin(delta_B)**2     # silent: ρ₀² sin²(δ_B) (∂n̂)²

print(f"  f_R² (radial)  = {f_R_sq:.8f}")
print(f"  f_A² (active)  = {f_A_sq:.8f}  = ρ₀²")
print(f"  f_S² (silent)  = {f_S_sq:.8f}  = ρ₀² sin²(δ_B)")
print()
print(f"  Ratio f_S²/f_A² = sin²(δ_B) = {np.sin(delta_B)**2:.6f}")
print(f"  Ratio f_S²/f_R² = ρ₀² sin²(δ_B) = {f_S_sq:.6f}")

# ----------------------------------------------------------------------
# Geometric quantities: areas/volumes of the mode submanifolds
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Geometric quantities of the mode submanifolds")
print("=" * 72)

# Active direction: 1D circle of phase ψ on S³ at fixed |ξ|.
# Its circumference depends on the parameterization; for the standard
# embedding, the circle is at radius ρ₀, so circumference = 2π ρ₀.
active_circumference = 2 * np.pi * rho_0
print(f"  Active (ψ) circumference         = 2π ρ₀ = {active_circumference:.6f}")

# Silent direction: S² of radius |Im ξ| = ρ₀ sin(δ_B)
# Area = 4π × (|Im ξ|)² = 4π ρ₀² sin²(δ_B) = 4π f_S²
silent_area = 4 * np.pi * f_S_sq
print(f"  Silent (n̂) S² area               = 4π f_S² = {silent_area:.6f}")
print(f"  Silent (n̂) S² radius             = ρ₀ sin(δ_B) = {np.sqrt(f_S_sq):.6f}")

# Full S³: volume = 2π² ρ₀³
S3_volume = 2 * np.pi**2 * rho_0**3
print(f"  Full S³ volume                  = 2π² ρ₀³ = {S3_volume:.6f}")

# ----------------------------------------------------------------------
# Candidate identifications with α, α_W, G
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Candidate identifications with physical couplings")
print("=" * 72)

alpha_emp = 1.0 / 137.035999
alpha_W_emp = 1.0 / 29.78        # ≈ g₂² / (4π) with g₂ ≈ 0.654
g2_emp = np.sqrt(4 * np.pi * alpha_W_emp)
G_e_emp = 1.75e-45                # dimensionless gravity coupling (Planck units)

print(f"  Empirical values:")
print(f"    α     = {alpha_emp:.6e}  ( = 1/137.036 )")
print(f"    α_W   = {alpha_W_emp:.6e}  ( ≈ g₂²/(4π), g₂ ≈ 0.654 )")
print(f"    g₂    = {g2_emp:.6f}")
print(f"    G_e   = {G_e_emp:.6e}")

print()
print("-" * 72)
print("(S) Silent SU(2)/U(1) sector — candidate for α_W")
print("-" * 72)
print()
print("Sigma-model decay constant: f_S² = ρ₀² sin²(δ_B)")
print()
print("Natural candidate identifications:")
print()
print(f"  α_W candidate (1):  f_S²                = {f_S_sq:.6f}")
print(f"  α_W candidate (2):  4π f_S² = Area      = {silent_area:.6f}")
print(f"  α_W candidate (3):  f_S² / (2π)         = {f_S_sq / (2*np.pi):.6f}")
print(f"  α_W candidate (4):  sin²(δ_B)            = {np.sin(delta_B)**2:.6f}")
print(f"  Empirical α_W                            = {alpha_W_emp:.6f}")
print()
print(f"  Ratios (empirical/candidate):")
print(f"    α_W / f_S²               = {alpha_W_emp / f_S_sq:.4f}")
print(f"    α_W / (4π f_S²)           = {alpha_W_emp / silent_area:.4f}")
print(f"    α_W / sin²(δ_B)           = {alpha_W_emp / np.sin(delta_B)**2:.4f}")

print()
print("-" * 72)
print("(A) Active U(1) sector — candidate for α")
print("-" * 72)
print()
print("Sigma-model decay constant: f_A² = ρ₀² = 1/2")
print()
print(f"  α candidate (1):  1/(4π f_A²) = 1/(2π)  = {1/(2*np.pi):.6f}")
print(f"  α candidate (2):  f_A² / V(S³)          = {f_A_sq / S3_volume:.6f}")
print(f"  α candidate (3):  f_S² / f_A² = sin²    = {np.sin(delta_B)**2:.6f}")
print(f"  Empirical α                              = {alpha_emp:.6f}")
print()
print(f"  Ratios (empirical/candidate):")
print(f"    α / (1/2π)                = {alpha_emp * 2 * np.pi:.6f}")
print(f"    α / (f_A²/V(S³))           = {alpha_emp / (f_A_sq/S3_volume):.4f}")
print(f"    α^(-1) / V(S³)             = {(1/alpha_emp) / S3_volume:.4f}")

print()
print("-" * 72)
print("Cross-sector ratio (v59 result): α_W / α")
print("-" * 72)
print()
print(f"  α_W / α    (empirical)  = {alpha_W_emp / alpha_emp:.6f}")
print(f"  f_S² / f_A² = sin²(δ_B)  = {np.sin(delta_B)**2:.6f}")
print(f"  Ratio (empirical / sin²)  = {(alpha_W_emp/alpha_emp) / np.sin(delta_B)**2:.4f}")
print()
print(f"  α_W / α  reformulated  : 137.036 / 29.78 = {137.036 / 29.78:.4f}")
print(f"  Compare to (1/sin² δ_B) / dim Spin(7)/14")
print(f"             = {1/np.sin(delta_B)**2 / (21/14):.4f}")

# ----------------------------------------------------------------------
# Specific Higgs-VEV style identification
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Higgs-mechanism style: m_W² = g₂² f² (with f = sigma decay constant)")
print("=" * 72)
print()
print("Standard SM:  m_W = g₂ v / 2  where v = Higgs VEV = 246 GeV.")
print("Predicted m_W from g₂ and Higgs VEV:")
m_W_emp = 80.4    # GeV
v_higgs = 246.0   # GeV
print(f"  m_W (emp) = {m_W_emp} GeV,  v (emp) = {v_higgs} GeV")
print(f"  g₂ = 2 m_W / v = {2 * m_W_emp / v_higgs:.6f}")
print()
print("In our v59 framework: f_S = ρ₀ sin(δ_B) is dimensionless. To make a")
print("dimensional W-boson mass, we'd need an external mass scale Λ:")
print(f"  m_W² ~ g₂² (Λ × f_S)² ?")
print(f"  Λ × f_S = v/2 = {v_higgs/2} GeV  →  Λ = {v_higgs/2 / np.sqrt(f_S_sq):.2f} GeV")
print(f"  Compare to Brannen mass scale a²: a² ≈ 314 MeV² → a ≈ 17.7 √MeV.")
print(f"  These scales differ by ~3 orders of magnitude — so the lepton-mass")
print(f"  Brannen scale is NOT the Higgs scale.  The two sectors decouple.")

# ----------------------------------------------------------------------
# What we can SAY about the prediction
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Honest assessment")
print("=" * 72)
print("""
What's now clean (numerical / Lean-verified):

  • Brannen amplitudes sₖ = a (1 + 2t cos(2πk/3 + φ)) — analytic.
  • Koide Q = (1 + 2t²)/3, hence Q = 2/3 ⟺ t² = 1/2.
  • dim Spin(7) = 21 = (7 choose 2); dim G₂ = 14 = 49 − 35.
  • Silent SU(2)/U(1) direction: spectrum invariant under SO(3) rotations of
    Im ξ (Lean theorem `SilentDirection.silent_pair`).

What's natural-scale (this script):

  • The silent direction is a 2-sphere of radius |Im ξ| = ρ₀ sin(δ_B) ≈ 0.156.
  • Its area is 4π f_S² ≈ 0.306.
  • The naive identification α_W = f_S² gives 0.0243 vs empirical 0.0336
    (factor 1.38 off — IN THE BALLPARK).
  • The α / α_W ratio comes out at sin²(δ_B) = 0.0486 in our framework,
    vs empirical 4.60 (factor of 90+ off — wrong sign? Or our identifications
    α↔A, α_W↔S need to be reversed).

What's needed to predict g₂ precisely:

  • A specific action functional with both YM kinetic term and matter kinetic
    term on the silent S² — only the RATIO of their coefficients gives g₂.
  • A natural mechanism setting that ratio (e.g. anomaly matching, integrating
    out heavy modes, gauge bundle topology).
  • A density coupling that gives a *specific* mode of the silent SU(2) — i.e.
    which combination of the 3 SU(2) generators is excited by ∇ρ.

What's worth investigating next:

  (a) Swap candidate identifications: maybe α (not α_W) is on the silent S²,
      and α_W is on the active circle.  Test: α_emp/α_W_emp ≈ 0.22.
      f_S²/f_A² = sin²(δ_B) ≈ 0.0486.  Ratio 0.22/0.0486 ≈ 4.5 — closer.
      But not exact.
  (b) Treat α_W as a Hopf-fibration holonomy on the silent S² (the S²
      has a natural Hopf bundle structure).  Predicts α_W from quantised
      Chern numbers.
  (c) Reverse engineer: what FORM of action gives α_W exactly from v59
      data?  This is dangerous numerology but can reveal the right shape.
""")

# Save the natural-scale data
np.savez('/home/d/code/scp/v59/cosserat_experiment/04_modes.npz',
         rho_0=rho_0, delta_B=delta_B,
         f_R_sq=f_R_sq, f_A_sq=f_A_sq, f_S_sq=f_S_sq,
         silent_area=silent_area, S3_volume=S3_volume,
         alpha_emp=alpha_emp, alpha_W_emp=alpha_W_emp, G_e_emp=G_e_emp,
         g2_emp=g2_emp)
print("Saved natural-scale data to 04_modes.npz")
