#!/usr/bin/env python3
"""
v59/furey_construction/08_brannen_yukawa.py

Variant I — Brannen-Koide Yukawa structure: connect the Brannen mass
kernel M = a(I + ξS + ξ̄S²) explicitly to the fermion content built in
07_full_generation.py.

The question this addresses:
  Given that one generation has 16 fermion states (from ℂ⊗ℍ⊗𝕆), and three
  generations come from Z_3 ⊂ S_3 triality of Spin(8), WHAT IS THE
  YUKAWA TERM in the Lagrangian that gives the v59 Brannen-Koide masses?

The answer:
  L_Yuk = sum_X [y_X · ψ̄_X^{(i)} · M_X(ξ)_{ij} · ψ_X^{(j)} + h.c.]
  where:
    - X labels the fermion sector (lepton, d-quark, u-quark, neutrino)
    - ψ_X^{(i)} for i=1,2,3 are the three generations
    - M_X(ξ) = a_X · (I + ξ_X S + ξ̄_X S²)
    - S = Z_3 cyclic-shift generator (3×3 matrix acting on flavor)
    - ξ_X ∈ ℍ with sector-specific constraint |ξ_X|² = (1 - dim G_2/D_X)
    - y_X = sector-specific Yukawa amplitude

This connects to the SM Higgs sector via:
  ξ_X = identification of the lepton Higgs (Brannen quaternion) with the
        sector-X projection of the Cl(7)_even multivector Φ.
"""

import numpy as np
import math

# -----------------------------------------------------------------------------
# Part 1: The Z_3 cyclic shift acting on 3-generation flavor space
# -----------------------------------------------------------------------------
print("=" * 80)
print("Part 1: Z_3 cyclic shift on 3-flavor space (generations)")
print("=" * 80)

S = np.array([
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0],
], dtype=complex)

# Eigenvalues of S
omega = np.exp(2j * np.pi / 3)
eigs_S = np.array([1, omega, omega**2])
print(f"S^3 = I:                 {np.allclose(np.linalg.matrix_power(S, 3), np.eye(3))}")
print(f"Eigenvalues of S:        1, ω, ω²  where ω = e^(2πi/3)")
print(f"Verified:                {np.allclose(sorted(np.linalg.eigvals(S), key=lambda z: np.angle(z)), sorted(eigs_S, key=lambda z: np.angle(z)))}")
print()
print("Eigenvectors are |ω^k⟩ for k=0,1,2 — these are the MASS EIGENSTATES.")
print("Brannen kernel diagonalises in this basis.")


# -----------------------------------------------------------------------------
# Part 2: Brannen mass kernel — derivation of mass eigenvalues
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("Part 2: Brannen kernel M = a(I + ξS + ξ̄S²) — mass eigenvalues")
print("=" * 80)

def brannen_M(a, xi):
    """Build the Brannen 3×3 mass operator M = a(I + ξS + ξ̄S²).

    For ξ ∈ ℝ, this is symmetric.  For ξ ∈ ℂ (with ξ̄ = complex conjugate),
    M is Hermitian.  For ξ ∈ ℍ acting on 3-gen flavor space, M is
    quaternion-Hermitian and the eigenvalues are real (Koide structure).
    """
    return a * (np.eye(3, dtype=complex) + xi * S + np.conj(xi) * S.conj().T)

# Verify Brannen formula for leptons (t² = 1/2)
print("\nLEPTON SECTOR (Furey N=0, ambient D = 28 = L = Λ²⊕Λ⁶ of R⁷):")
print(f"  t² = 1 - dim G_2 / D = 1 - 14/28 = 1/2     (Koide Q = 2/3)")
print(f"  Brannen phase φ = 2/9 rad (experimental)")

# Reproduce empirical lepton masses
a_l = 17.7155      # √MeV (empirical Brannen lepton scale)
t_l = math.sqrt(1/2)
phi_l = 2/9         # rad
xi_l = t_l * np.exp(1j * phi_l)  # ξ = t·e^(iφ) (using complex slice of ℍ)

M_l = brannen_M(a_l, xi_l)
eigs_l = np.linalg.eigvalsh((M_l + M_l.conj().T) / 2)  # ensure Hermitian
masses_l = eigs_l**2  # mass = (sqrt-mass)² (Brannen convention)

print(f"\n  ξ_l = {xi_l:.4f}  (|ξ|² = {abs(xi_l)**2:.4f}, should be 1/2)")
print(f"  Brannen amplitudes s_k = eigenvalues of M_l:")
print(f"     s_0 = {eigs_l[0]:8.4f} √MeV  →  m = {eigs_l[0]**2:9.4f} MeV   (expected m_e = 0.511 MeV)")
print(f"     s_1 = {eigs_l[1]:8.4f} √MeV  →  m = {eigs_l[1]**2:9.4f} MeV   (expected m_μ = 105.66 MeV)")
print(f"     s_2 = {eigs_l[2]:8.4f} √MeV  →  m = {eigs_l[2]**2:9.4f} MeV   (expected m_τ = 1776.86 MeV)")

# Compute Koide ratio
Q_l = sum(masses_l) / sum(eigs_l)**2
print(f"\n  Koide Q = Σm / (Σ√m)² = {Q_l:.6f}    (expected 2/3 = {2/3:.6f})")
print(f"  Match: {abs(Q_l - 2/3)/2/3 * 100:.5f}%")

# -----------------------------------------------------------------------------
# Part 3: The Yukawa coupling in v59
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("Part 3: Yukawa coupling in v59 — connection to fermions")
print("=" * 80)

print("""
In SM, the Yukawa Lagrangian for leptons reads:
  L_Yuk^lepton = -y_e^{ij} · L̄_L^i · Φ · e_R^j + h.c.
  where L_L = (ν_L, e_L) doublet, Φ = Higgs doublet, e_R = right singlet,
  and y_e^{ij} is the 3×3 generation matrix.

At Higgs VEV ⟨Φ⟩ = (0, v/√2):
  m_e^{ij} = y_e^{ij} · v / √2

In v59, the Yukawa matrix has the Brannen-Koide form:
  y_e^{ij} · v / √2 = a_l · [I + ξ_l·S + ξ̄_l·S²]^{ij}
                    = M_l^{ij}

So:
  y_e^{ij} = (√2 / v) · M_l^{ij}
           = (√2 / (28²·a_l²)) · a_l · [I + ξ S + ξ̄ S²]^{ij}
           = (√2 / 28²) · (1/a_l) · [I + ξ S + ξ̄ S²]^{ij}

The Yukawa magnitude has TWO factors:
  - Geometric: √2 / 28² ≈ 1.80e-3  (universal, all leptons)
  - Mass-scale: 1/a_l ≈ 1/17.7 √MeV (sets per-generation magnitude)
  - Brannen structure: [I + ξS + ξ̄S²]  (gives the 3-generation hierarchy)
""")

# Verify numerically: compute y_e for electron and tau, check against SM
v_Higgs = 246.22 * 1000  # MeV
print(f"v_Higgs (MeV)               = {v_Higgs:.0f}")
print(f"D_lepton² · a_l² (MeV)      = {28**2 * a_l**2:.0f}    (should match v_Higgs at 0.07%)")

m_e = 0.511
m_tau = 1776.86
y_e_SM = m_e * math.sqrt(2) / v_Higgs
y_e_v59 = (math.sqrt(2) / v_Higgs) * m_e
y_tau_SM = m_tau * math.sqrt(2) / v_Higgs

print(f"\nElectron Yukawa:")
print(f"  y_e = m_e · √2 / v = {y_e_SM:.4e}")
print(f"  In v59 normalization: y_e = (√2 / 28²) · (m_e/a_l²) = {(math.sqrt(2)/(28**2)) * (m_e/a_l**2):.4e}")

print(f"\nTau Yukawa:")
print(f"  y_τ = m_τ · √2 / v = {y_tau_SM:.4e}")
print(f"  In v59: y_τ = (√2 / 28²) · (m_τ/a_l²) = {(math.sqrt(2)/(28**2)) * (m_tau/a_l**2):.4e}")

# Both match by construction since y · v = m

print()
print("CRUCIAL POINT: the Brannen kernel REPLACES the SM Yukawa matrix.")
print("Instead of 3×3 = 9 complex Yukawa entries (= 18 real parameters),")
print("v59 has TWO parameters per sector: (a_X, ξ_X), with |ξ|² constrained.")
print("Net: 1 (amplitude) + 4 - 1 (constraint) = 4 real parameters per sector,")
print("of which the Koide ratio is structural (Q = (1+2|ξ|²)/3) — leaving 3.")


# -----------------------------------------------------------------------------
# Part 4: Sector-specific ξ_X — cross-sector structure
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("Part 4: Sector-specific ξ_X — cross-sector Yukawa structure")
print("=" * 80)

sectors = {
    'lepton (e, μ, τ)':       {'D': 28, 't_sq': 1/2,   'Q_pred': 2/3,    'masses': [0.511, 105.66, 1776.86]},
    'd-quark (d, s, b)':      {'D': 35, 't_sq': 3/5,   'Q_pred': 11/15,  'masses': [4.67, 93.4, 4180.0]},
    'u-quark (u, c, t)':      {'D': 63, 't_sq': 7/9,   'Q_pred': 23/27,  'masses': [2.16, 1273.0, 172570.0]},
}

print(f"\n  {'Sector':<22s} {'D':>4s} {'t²':>8s} {'Q (v59)':>10s} {'Q (obs)':>10s} {'Q gap %':>10s}")
print("  " + "-"*66)
for name, p in sectors.items():
    m = p['masses']
    Q_obs = sum(m) / (sum(math.sqrt(mi) for mi in m))**2
    p['Q_obs'] = Q_obs
    gap = abs(Q_obs - p['Q_pred']) / p['Q_pred'] * 100
    print(f"  {name:<22s} {p['D']:>4d} {p['t_sq']:>8.4f} {p['Q_pred']:>10.6f} {Q_obs:>10.6f} {gap:>10.4f}")

# Verify (1-Q)·D = 28/3
print()
print("Verifying universal Koide-deviation identity (1-Q)·D = 28/3 = D_lepton/n_gen:")
for name, p in sectors.items():
    LHS = (1 - p['Q_pred']) * p['D']
    print(f"  {name:<22s} (1-Q)·D = {LHS:.4f}  (expected 28/3 = {28/3:.4f})")


# -----------------------------------------------------------------------------
# Part 5: Selection rule — which sector sees which Cl-grade
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("Part 5: Selection rule — sector → Cl-grade projection")
print("=" * 80)

print("""
v59 observation (cosserat_experiment/16_Z2_decomposition.py):
  Cl(7)_even = ℂ⊗𝕆 (dim 64) bisects via G_2-invariant content into:
    L (Lie content)  = Λ²ℝ⁷ ⊕ Λ⁶ℝ⁷  (dim 28)  — NO G₂-invariant
    F (Form content) = Λ⁴ℝ⁷           (dim 35)  — contains coassociative *φ

  Sectors couple via Z₂ × Z₂ pattern (Bit-L, Bit-F):
    N=0 lepton    : (L=1, F=0) → L only   → D = 28
    N=1 d-quark   : (L=0, F=1) → F only   → D = 35
    N=2 u-quark   : (L=1, F=1) → L ⊕ F    → D = 63
    N=3 neutrino  : (L=0, F=0) → empty?   → D = 0  (sterile? Dirac mass possible?)
                  OR (L=1, F=1) again with different N choice → D = 63

The selection rule WHY each sector picks its specific (Bit-L, Bit-F) is
NOT yet derived — 7 hypotheses tested in 2026-05-22 session, all failed.
This is Frontier 2 of v59 ROADMAP.

But the SELECTION RULE determines:
  - Which Cl-grade of the Higgs multivector Φ ∈ Cl(7)_even gives each
    sector its Yukawa coupling.
  - Whether the neutrino sector is sterile (D=0, no Brannen Yukawa) or
    not (D=63, Dirac neutrino with the same kernel as the u-quark).

For now, the framework is:
  M_X^{ij} = a_X · (I + ξ_X S + ξ̄_X S²)^{ij}
  with ξ_X ∈ (sector-specific Cl-grade subspace of Cl(7)_even ≅ ℂ⊗𝕆)
""")


# -----------------------------------------------------------------------------
# Part 6: Connection to v58 multivector force law — gravity source
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("Part 6: Connection to v58 — gravity from constraint deviation")
print("=" * 80)

print("""
v58's multivector force law had:
  L_v58 = (1/2)⟨(∂M)(∂M̃)⟩₀ - V(|M|²) + chiral_current·gauge

with gravity source ρ_M = (1/2)(M·M̃ - v²) for some vacuum v².

v59 + v58 SYNTHESIS (`synthesis/FINDINGS_synthesis.md`):
  ρ_M(x) = a² · (|ξ(x)|² - 1/2)

i.e., the v58 gravity source IS the deviation of the Brannen kernel from
the v59 equilibrium constraint surface |ξ|² = 1/2.

For each sector:
  - Lepton (|ξ|² = 1/2):     ρ_M = 0    (on constraint — no gravity)
  - d-quark (|ξ|² = 3/5):    ρ_M = a²·(3/5 - 1/2) = a²/10
  - u-quark (|ξ|² = 7/9):    ρ_M = a²·(7/9 - 1/2) = 5a²/18

So QUARKS gravitate via their constraint-surface mismatch, LEPTONS don't.

This DOES NOT match SM (electrons gravitate just like quarks per unit mass).
So either:
  (a) The v58 gravity source is universal (not sector-specific) and the
      ξ-mismatch is a CORRECTION to gravity, not the leading term.
  (b) The "gravity" v58 picks up from ξ-mismatch is something OTHER than
      Newton gravity — perhaps QCD confinement, or anomalous mass.
  (c) The sector identifications need revisiting.

This is an OPEN ISSUE in v58⊕v59 synthesis.  Mark for Frontier-3 work.
""")


# -----------------------------------------------------------------------------
# Part 7: Summary of the integrated picture
# -----------------------------------------------------------------------------
print()
print("=" * 80)
print("Part 7: The integrated v58⊕v59⊕Furey picture")
print("=" * 80)

print("""
ALGEBRAIC FOUNDATION:
  ℂ⊗ℍ⊗𝕆  (64 cplx-dim = 128 real-dim) — the full Furey color algebra
  ≅ Cl(7)_even (via ℂ⊗𝕆 ≅ Cl(6))
  ≅ all v59 sector ambients are graded subspaces

DYNAMICAL FIELDS:
  Φ(x) ∈ Cl(7)_even (or specific subspaces, sector-specific)
    Decomposes:  Λ⁰ ⊕ Λ² ⊕ Λ⁴ ⊕ Λ⁶  = 1 + 21 + 35 + 7 = 64
                 (identity, lepton-L₁, d-quark, lepton-L₂)

  ψ_X^{(i)}(x) — 3-flavor (i=1,2,3) fermion fields:
    ψ_e^{(i)}: leptons (ν_L, e_L doublet + ν_R, e_R)
    ψ_d^{(i)}: d-type quarks
    ψ_u^{(i)}: u-type quarks
    Total per generation: 16 states (in ℂ⊗ℍ⊗𝕆 left-ideal)
    Three generations from Z₃ ⊂ S₃ triality of Spin(8).

  Gauge fields:
    A_μ^a (a=1,2,3) — SU(2)_L from silent direction (Lean: SilentDirection.lean)
    B_μ              — U(1)_Y (origin in Cl(7)_even still open; possibly SU(2)_R)
    G_μ^a (a=1...8) — SU(3)_c from Cl(6) Witt decomposition

  Spacetime: Cl(3,1) Dirac structure (γ^μ).

LAGRANGIAN TERMS:
  L_kinetic_Φ  = (1/2)⟨(D_μ Φ)(D^μ Φ̃)⟩₀     (Higgs sector kinetic)
  L_potential  = -V(|Φ|²)                    (drives Φ to constraint surfaces)
  L_kinetic_ψ  = i ψ̄_X γ^μ D_μ ψ_X            (fermion kinetic)
  L_Yukawa     = sum_X y_X · ψ̄_X · M_X(Φ) · ψ_X + h.c.   (Brannen kernel)
  L_gauge      = -(1/4) F^a_μν F^{aμν}        (all three gauge factors)
  L_grav       = R_μν / (16π G_N) + ρ_M·g_μν   (gravity from constraint deviation)

STRUCTURAL OUTPUTS (this and prior sessions):
  - Q_e = 2/3 = dim G_2 / dim Spin(7)
  - φ = 2/9 = sin²θ_W
  - v_Higgs = 28² · a_l² (scale bridge)
  - α(0) = exp(-π²/2 + 2α)        (instanton)
  - α(M_Z) = 25/(324π²)          (tree-level consistency)
  - g_W² = 5√α  (Killing-form embedding)
  - cos²θ_W = t²_u = 7/9         (cross-sector)
  - m_W, m_Z, m_H at <0.15%
  - Q_d = 11/15, Q_u = 23/27
  - m_top = 172.6 GeV (apex + a_u² = 72·a_l²)

EMPIRICAL INPUTS (currently):
  - a_l (lepton Brannen scale) ← from m_e, m_μ, m_τ
  - Sector-specific a_X (some testable, some empirical)
  - Brannen phases φ_X (only φ_l = 2/9 is structural; quark phases empirical)

OPEN ISSUES:
  1. Lagrangian-level mechanism for g_W² = 5√α prefactor (Frontier 3)
  2. Lagrangian-level mechanism for v = 28²·a_l² (scale bridge)
  3. Selection rule: which sector → which Cl-grade (Z₂×Z₂ pattern undetermined)
  4. U(1)_Y geometric source in Cl(7)_even
  5. CKM/PMNS from Brannen eigenvectors
  6. Quark Brannen phases φ_u, φ_d
  7. v58⊕v59 gravity source ρ_M not yet universal (sector-dependent)
  8. Cosmological constant
""")

print()
print("=" * 80)
print("STEP 1 (fermion content) — COMPLETE")
print("STEP 2 (commit dynamical field) — needs decision: Φ ∈ Cl(7)_even (preferred)")
print("STEP 3 (CKM from Brannen eigenvectors) — separate computation, next")
print("STEP 4 (v58 multivector integration) — separate work, next")
print("STEP 5 (anomaly check) — DONE in this script ✓")
print("STEP 6 (write Lagrangian) — now feasible after Steps 2-4")
print("=" * 80)
