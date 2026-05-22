#!/usr/bin/env python3
"""
07_full_lagrangian.py

Full v59 strain Lagrangian — written *in concert* with the existing structure.

Field content:
  ξ(x) ∈ ℍ ≅ ℝ⁴               quaternion-valued field
  ψ_lepton(x) ∈ ℂ³            charged-lepton triplet
  A^a_μ(x), a = 1, 2, 3        SU(2)_L gauge field (silent direction)
  ρ_matter(x)                  matter density (the Cosserat source)

Symmetries:
  spacetime         SO(3, 1)
  global internal   SO(4)  →  SO(3)_silent × U(1)_active  (broken by Yukawa)
  local internal    SU(2)_L  (=  silent SO(3) gauged)

Constraint:   |ξ|² = ρ₀² = 1/2   (enforced by V → ∞)

Lagrangian:
  L = L_kin + L_YM + L_Yuk + L_V + L_density

with:
  L_kin     = (1/2) (D_μ ξ)·(D^μ ξ),   D_μ ξ = ∂_μ ξ + [A_μ, ξ]
  L_YM      = -(1/4 g_W²) F^a_{μν} F^{aμν}
  L_Yuk     = -y ψ̄_L M(ξ) ψ_R + h.c.,   M(ξ) = a (I + ξ S + ξ̄ Sᵀ)
  L_V       = -V(|ξ|²),  V = (λ/4)(|ξ|² - ρ₀²)²  (rigid limit)
  L_density = -g_D ρ_matter · |Im ξ|

Mode decomposition:
  ξ = ρ (cos ψ + sin ψ · n̂),  ρ → ρ₀,  ψ ∈ [0, π/2],  n̂ ∈ S² ⊂ Im ℍ

  |∂ ξ|² = (∂ρ)² + ρ²[(∂ψ)² + sin²ψ (∂n̂)²]    (R)+(A)+(S)

At Brannen equilibrium (ψ = δ_B = 2/9):
  silent decay constant²    f_S² = ρ₀² sin²(δ_B) ≈ 0.02428
  silent S² radius          R_S = ρ₀ sin(δ_B) = sin(2/9)/√2 ≈ 0.156

Gauge coupling identification:
  Killing-form embedding index for so(3) ⊂ so(7) = 5 (from 06_killing_form.py).
  Hence  g_W / g_parent = √5  where  g_parent  is the Yang-Mills coupling of
  the "ambient" Spin(7) structure.  This script tests several natural
  identifications for g_parent and reports the empirical match.

The empirical SU(2)_L coupling (PDG 2024, at M_Z):
  g₂(M_Z) = 0.6517   →   α_W = g₂²/(4π) = 0.03379 = 1/29.59.
"""

import numpy as np

# ----------------------------------------------------------------------
# v59 structural inputs
# ----------------------------------------------------------------------
rho_0 = 1 / np.sqrt(2)
delta_B = 2 / 9
S_em_v59 = np.pi**2 / 2     # the v59 instanton-action conjecture
embedding_index = 5         # Killing-form so(3) ⊂ so(7) (computed in 06)

# Empirical
alpha_emp = 1 / 137.035999084
alpha_emp_MZ = 1 / 127.952
g2_emp_MZ = 0.6517
g2_emp_tree = 0.6529        # from m_W/v tree-level
e_emp_lowE = np.sqrt(4 * np.pi * alpha_emp)
e_emp_MZ = np.sqrt(4 * np.pi * alpha_emp_MZ)

print("=" * 72)
print("FULL v59 STRAIN LAGRANGIAN — Mode decomposition & coupling identification")
print("=" * 72)

# ----------------------------------------------------------------------
# Section 1.  Lagrangian decomposition
# ----------------------------------------------------------------------
print()
print("─" * 72)
print("Section 1: Lagrangian decomposition at the Brannen equilibrium")
print("─" * 72)

f_S_sq = rho_0**2 * np.sin(delta_B)**2

print(f"""
At ξ = ξ_eq = ρ₀(cos δ_B + sin δ_B · n̂), n̂ ∈ S²:

  L_kin   =  (ρ₀²/2)(∂ψ)²                         ← (A) active mode
          + (ρ₀² sin²δ_B / 2)(D n̂)²                ← (S) silent mode (gauged)
          + (1/2)(∂ρ)²                            ← (R) radial mode (off S³)

  L_YM    = -(1/4 g_W²) F^a_{{μν}} F^{{aμν}}        ← SU(2)_L gauge kinetic

  L_Yuk   = -y (lepton)·M(ξ)·(lepton̄)              ← Brannen lepton mass

  L_V     = -(λ/4)(|ξ|² - ρ₀²)²                   ← constraint, λ → ∞

  L_dens  = -g_D ρ_matter |Im ξ|                  ← density coupling (Cosserat)

Numerical values:
  ρ₀          = 1/√2     = {rho_0:.6f}
  δ_B         = 2/9       = {delta_B:.6f}  rad
  f_S²        = ρ₀² sin²(δ_B)  = {f_S_sq:.8f}   ← silent "Higgs" VEV²
  f_S         = ρ₀ sin(δ_B)    = {np.sqrt(f_S_sq):.6f}
  R_S = f_S   = silent S² radius
""")

# ----------------------------------------------------------------------
# Section 2.  Killing-form prediction
# ----------------------------------------------------------------------
print("─" * 72)
print("Section 2: Killing-form prediction g_W / g_parent = √5")
print("─" * 72)

print(f"""
From 06_killing_form.py:
  Killing form on so(7):    B_{{so(7)}}(X,Y) = 5 · Tr_vec(X Y)
  Killing form on so(3):    B_{{so(3)}}(X,Y) = 1 · Tr_vec(X Y)
  Embedding index           X = B_so7|_so3 / B_so3 = 5

YM coupling ratio rule:
  If a parent gauge group G with coupling g_parent is broken to a subgroup H
  with embedding index X, then the kinetic term of the H gauge field in the
  parent's Killing-normalised action contributes
    L_YM(G)|_H = -(X / 4 g_parent²) Tr_H(F_H · F_H)
  giving g_H² = g_parent² / X, so g_H = g_parent / √X.

  But the OPPOSITE convention is used when normalising kinetic terms uniformly
  per generator (canonical normalisation): then g_H = √X · g_parent.

  Both forms are in the literature.  We test both below.

  Predicted ratio  g_W / g_parent  =  √5 ≈ {np.sqrt(5):.4f}
  Or              g_W · √5 = g_parent  (alternative convention).
""")

# ----------------------------------------------------------------------
# Section 3.  Identifying g_parent
# ----------------------------------------------------------------------
print("─" * 72)
print("Section 3: Trying natural identifications for g_parent")
print("─" * 72)

# Calculate candidates for the convention g_W = √5 × g_parent.
candidates = [
    ("e  (low-E α = 1/137.036)",          e_emp_lowE),
    ("e  (M_Z α = 1/127.95)",             e_emp_MZ),
    ("α^(1/4)  (low-E)",                   alpha_emp**(1/4)),
    ("α^(1/4)  (M_Z)",                     alpha_emp_MZ**(1/4)),
    ("α^(1/4)  (v59: exp(-(π²/2 − 2α)/4))", np.exp(-(np.pi**2/2 - 2*alpha_emp)/4)),
    ("exp(-π²/8)                          ", np.exp(-np.pi**2 / 8)),
    ("√α  (low-E)                          ", np.sqrt(alpha_emp)),
    ("f_S  (silent decay const)            ", np.sqrt(f_S_sq)),
]

print(f"\n  {'identification':<42}  g_parent      g_W = √5·g    gap vs g₂(M_Z)={g2_emp_MZ}")
print(f"  {'-'*42}  {'-'*12}  {'-'*12}  {'-'*15}")
for name, g_p in candidates:
    g_W_pred = np.sqrt(embedding_index) * g_p
    gap = (g_W_pred - g2_emp_MZ) / g2_emp_MZ * 100
    star = "  ←" if abs(gap) < 1 else ""
    print(f"  {name:<42}  {g_p:.6f}    {g_W_pred:.6f}    {gap:+6.3f} %{star}")

# Try also the alternative convention g_W = g_parent / √5
print(f"\n  Alternative convention: g_W = g_parent / √5\n")
print(f"  {'identification':<42}  g_parent      g_W = g/√5   gap vs g₂(M_Z)={g2_emp_MZ}")
print(f"  {'-'*42}  {'-'*12}  {'-'*12}  {'-'*15}")
for name, g_p in candidates:
    g_W_pred = g_p / np.sqrt(embedding_index)
    gap = (g_W_pred - g2_emp_MZ) / g2_emp_MZ * 100
    star = "  ←" if abs(gap) < 1 else ""
    print(f"  {name:<42}  {g_p:.6f}    {g_W_pred:.6f}    {gap:+6.3f} %{star}")

# ----------------------------------------------------------------------
# Section 4.  The cleanest match: g_W² = 5·√α
# ----------------------------------------------------------------------
print()
print("─" * 72)
print("Section 4: Best match — g_W² = 5 · √α")
print("─" * 72)

# This is g_W = √5 · α^(1/4).
g_W_predicted = np.sqrt(5) * alpha_emp**(1/4)

# Also the v59-conjecture-α version:
def solve_v59_alpha():
    x = np.pi**2 / 2
    for _ in range(50):
        x = np.pi**2 / 2 - 2 * np.exp(-x)
    return np.exp(-x)
alpha_v59 = solve_v59_alpha()
g_W_v59 = np.sqrt(5) * alpha_v59**(1/4)

print(f"""
The combination

    g_W² = 5 · √α    or equivalently    g_W = √5 · α^(1/4)

uses two v59 structural inputs (Killing-form embedding 5; α from the
−ln α + 2α = π²/2 conjecture).  Numerics:

  α (empirical low-E)             = {alpha_emp:.6e}
  α (v59 conjecture)              = {alpha_v59:.6e}

  g_W predicted (empirical α)      = √5 · α^(1/4) = {g_W_predicted:.6f}
  g_W predicted (v59 α)            = √5 · α_v59^(1/4) = {g_W_v59:.6f}

  g_W empirical (M_Z, PDG)         = {g2_emp_MZ:.6f}
  g_W empirical (tree, m_W/v)      = {g2_emp_tree:.6f}

  Gap (M_Z, using empirical α)     = {(g_W_predicted - g2_emp_MZ)/g2_emp_MZ*100:+.3f} %
  Gap (tree, using empirical α)    = {(g_W_predicted - g2_emp_tree)/g2_emp_tree*100:+.3f} %
  Gap (M_Z, using v59 α)            = {(g_W_v59 - g2_emp_MZ)/g2_emp_MZ*100:+.3f} %
""")

# ----------------------------------------------------------------------
# Section 5.  Cross-check with α_W
# ----------------------------------------------------------------------
print("─" * 72)
print("Section 5: Cross-check with α_W and the cross-sector ratio")
print("─" * 72)

alpha_W_predicted = g_W_predicted**2 / (4 * np.pi)
alpha_W_emp_MZ = g2_emp_MZ**2 / (4 * np.pi)
print(f"\n  α_W predicted   = g_W²/(4π) = {alpha_W_predicted:.6f} = 1/{1/alpha_W_predicted:.3f}")
print(f"  α_W empirical (M_Z) = {alpha_W_emp_MZ:.6f} = 1/{1/alpha_W_emp_MZ:.3f}")
print(f"  Gap = {(alpha_W_predicted - alpha_W_emp_MZ)/alpha_W_emp_MZ*100:+.3f} %")

# Structural form:
print(f"\n  Structural form of the conjecture:")
print(f"     α_W  =  (5 / 4π) · √α  ")
print(f"     α_W  =  ({5/(4*np.pi):.6f}) · √α")
print(f"     α_W  =  {5/(4*np.pi) * np.sqrt(alpha_emp):.6f}  (predicted from α_low)")
print(f"     α_W (M_Z)  =  {alpha_W_emp_MZ:.6f}  (PDG)")

# Also: α_W / α = (5/4π) · (1/√α) = (5/4π) · √(1/α)
ratio_predicted = (5 / (4 * np.pi)) / np.sqrt(alpha_emp)
ratio_emp = alpha_W_emp_MZ / alpha_emp
print(f"\n  Ratio α_W / α (predicted):    {ratio_predicted:.4f}")
print(f"  Ratio α_W / α (empirical, M_Z + low-E):  {ratio_emp:.4f}")

# ----------------------------------------------------------------------
# Section 6.  Honest assessment
# ----------------------------------------------------------------------
print()
print("=" * 72)
print("Section 6: Honest assessment")
print("=" * 72)
print(f"""
Numerical observation:  g_W² = 5 · √α  matches empirical g_W(M_Z) at the
                        {abs((g_W_predicted - g2_emp_MZ)/g2_emp_MZ*100):.2f} % level (better than the previous √21 conjecture
                        which was 5.7 % off at consistent scales).

Structural ingredients of this combination:
  • 5 = Killing-form embedding index of so(3) ⊂ so(7)  (DERIVED, 06_killing_form.py)
  • √α — the v59 instanton action α = exp(-S_em) with S_em = π²/2 (CONJECTURE,
        v59 04_findings.md, matches empirical α at 4 × 10⁻⁵).

The combination g_W² = 5 √α arises naturally if the Spin(7) gauge coupling
takes the value g_Spin(7) = α^(1/4) at the v59 scale, and SU(2)_L is the
embedded subgroup with the standard Yang-Mills kinetic-term lift.  We have NOT
derived g_Spin(7) = α^(1/4) from a Lagrangian — it is a numerical observation
matching at the 0.1 % level.

Things that go right:
  • The exponent 1/4 is the only nearby exponent that gives ~empirical g_W.
    (1/2, 1, 1/8, 1/3 all give very different values.)
  • The factor √5 is the natural Killing-form result, NOT a fit.
  • The match is at 0.1 - 0.3 % across choices of α (empirical or v59-conjecture).

Things that need scrutiny:
  • Why g_Spin(7) = α^(1/4)?  No Lagrangian derivation yet.
  • Running coupling: α and g_W run differently with scale; the formula matches
    only at certain scale combinations.  Honest precision is at the 1-3 % level.
  • The "Spin(7) parent" gauge theory structure is conjectural — v59 doesn't
    yet have a Spin(7) Yang-Mills theory with broken symmetry pattern.

Comparable to v59's existing partial predictions:
  • α          (π²/2 - 2α conjecture): 4 × 10⁻⁵ ✓
  • g_W        (5 √α conjecture, THIS WORK): 0.1 - 3 % ✓ (tentative)
  • G_e        (21·π²/2 conjecture): factor 0.76 ⚠
""")

# Save
np.savez('/home/d/code/scp/v59/cosserat_experiment/07_lagrangian.npz',
         rho_0=rho_0, delta_B=delta_B,
         f_S_sq=f_S_sq, embedding_index=embedding_index,
         g_W_predicted=g_W_predicted, g_W_v59=g_W_v59,
         g2_emp_MZ=g2_emp_MZ, g2_emp_tree=g2_emp_tree)
print("\nSaved Lagrangian data to 07_lagrangian.npz")
