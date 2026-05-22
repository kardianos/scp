#!/usr/bin/env python3
"""
v59 Cosserat-coupling experiment.

The hypothesis (user's framing, 2026-05-22):
    The Brannen phase φ is not a global constant — it is a *local field* φ(x).
    A local-density rotation of the lepton-kernel frame shifts the spectral
    decomposition (eigenvalues of the cyclic mass operator), inducing a Cosserat-
    style strain.  Simple/scalar test masses don't notice the rotation; the Z₃
    cyclic kernel does (its eigenvalues are flavor-dependent functions of φ).

This script runs three steps:
    1. Fit (a, δ) of the Brannen form to PDG lepton masses; verify Koide = 2/3.
    2. Compute response coefficients ∂(ln m_k)/∂φ for τ, μ, e.
    3. Solve a static, spherically-symmetric Cosserat field equation for φ(r)
       around a proton density bump, then convert to predicted lepton mass
       shifts and compare to experimental bounds (atomic clocks, EP).

Output: terminal report + saves data arrays to .npz for further analysis.
"""

import numpy as np
from scipy.integrate import solve_bvp

# ----------------------------------------------------------------------
# PDG 2024 lepton masses (MeV)
# ----------------------------------------------------------------------
m_e = 0.5109989461
m_mu = 105.6583755
m_tau = 1776.86

# Brannen ordering: heaviest → k=0, then 2π/3 and 4π/3 phase shifts.
# Empirical assignment (from v59 fit): k=0 → τ, k=1 → e, k=2 → μ.
lepton_names = ('τ', 'e', 'μ')
lepton_masses = np.array([m_tau, m_e, m_mu])

# ----------------------------------------------------------------------
# Step 1.  Brannen fit.
# ----------------------------------------------------------------------
print("=" * 70)
print("Step 1: Brannen fit to PDG lepton masses")
print("=" * 70)

sqrt_m = np.sqrt(lepton_masses)
a = sqrt_m.sum() / 3.0  # the overall scale
# The Brannen phase δ is determined by α_τ = δ:
# cos(δ) = (√m_τ - a) / (√2 · a)
cos_delta = (sqrt_m[0] - a) / (np.sqrt(2.0) * a)
delta = np.arccos(cos_delta)
print(f"  Brannen scale  a = {a:.6f}  (√MeV)")
print(f"  Brannen phase  δ = {delta:.8f} rad")
print(f"                   = {delta * 9 / 2:.10f} (in units of 2/9)")
print(f"  Empirical 2/9     = {2/9:.10f}")
print(f"  δ - 2/9          = {delta - 2/9:+.3e} rad")

# Construct α_k for k=0,1,2
alpha = np.array([delta, 2*np.pi/3 + delta, 4*np.pi/3 + delta])
sqrt_m_pred = a * (1.0 + np.sqrt(2.0) * np.cos(alpha))
m_pred = sqrt_m_pred ** 2

print(f"\n  Reconstructed Brannen amplitudes √m_k = a(1 + √2 cos(α_k)):")
for i, name in enumerate(lepton_names):
    rel_err = (m_pred[i] - lepton_masses[i]) / lepton_masses[i]
    print(f"    {name}: m_pred = {m_pred[i]:12.6f} MeV, "
          f"m_obs = {lepton_masses[i]:12.6f} MeV, rel-err = {rel_err:+.3e}")

# Koide sanity check
Q = lepton_masses.sum() / sqrt_m.sum()**2
print(f"\n  Koide Q = Σm / (Σ√m)² = {Q:.12f}")
print(f"  2/3                     = {2/3:.12f}")
print(f"  Q - 2/3                 = {Q - 2/3:+.3e}")

# ----------------------------------------------------------------------
# Step 2.  Response coefficients ∂(ln m_k) / ∂φ.
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 2: Response coefficients ∂(ln m_k) / ∂φ")
print("=" * 70)
print("\n  m_k = (a (1 + √2 cos α_k))² = (√m_k)², α_k = 2πk/3 + δ.")
print("  ∂(√m_k)/∂φ = -√2 a sin α_k.")
print("  ∂(ln m_k)/∂φ = -2 √2 a sin α_k / √m_k.\n")

response = -2.0 * np.sqrt(2.0) * a * np.sin(alpha) / sqrt_m
print(f"  Response coefficients (dimensionless, per rad of δφ):")
for i, name in enumerate(lepton_names):
    print(f"    ∂(ln m_{name})/∂φ = {response[i]:+10.4f}")

# Sum of responses (should reflect Σ m_k invariance modulo phase invariance)
print(f"\n  Sum of m_k × response_k: {(lepton_masses * response).sum():+.3e}")
print(f"  (Compare Σ ∂ m_k / ∂φ from direct: ", end='')
direct_sum = (-2 * np.sqrt(2) * a * np.sin(alpha) * sqrt_m).sum()
print(f"{direct_sum:+.3e})")

# Save coefficients
np.savez('/home/d/code/scp/v59/cosserat_experiment/01_fit.npz',
         a=a, delta=delta, alpha=alpha,
         m_obs=lepton_masses, m_pred=m_pred,
         sqrt_m=sqrt_m, sqrt_m_pred=sqrt_m_pred,
         response=response, lepton_names=lepton_names, Q=Q)

# ----------------------------------------------------------------------
# Step 3.  Cosserat field around a proton-like density bump.
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 3: Static Cosserat field φ(r) around a proton")
print("=" * 70)
print("""
The φ-field equation (radial, spherical):
    -∇²φ + m_φ² (φ - δ) = -λ · (ρ - ρ_bg) / ρ_bg
    -1/r² d/dr(r² dφ/dr) + m_φ² (φ - δ) = -λ · (ρ(r) - ρ_bg) / ρ_bg

ρ(r) = ρ_p exp(-(r/R_p)²) + ρ_bg (Gaussian proton plus background)
R_p = 0.84 fm (proton charge radius)
Boundary: φ'(0) = 0 (regularity), φ(L) = δ (far-field equilibrium).
""")

# Numerical setup
R_p = 0.84  # fm (proton charge radius)
L = 30.0    # fm (radial extent)
rho_p = 1.0  # arbitrary density scale (will cancel via δρ/ρ_bg)
rho_bg = 1e-6  # very small background — δρ/ρ_bg ~ 1/ρ_bg at the center
                # → large source magnitude at proton center.
                # We choose this so that δρ/ρ_bg is order unity in the
                # surrounding "vacuum" near the bump's edge, representing
                # the typical scale of density variation in a soliton context.

# Actually a cleaner formulation: just use the dimensionless source
# S(r) = exp(-(r/R_p)²), peak amplitude 1 at the center, decays away.
# Then the equation is:
#   -∇²φ + m_φ² (φ - δ) = -λ · S(r)
# with λ now setting the overall coupling strength.

def source(r):
    return np.exp(-(r/R_p)**2)

# We test several m_φ values (the mediator mass).
# m_φ = 0 (massless, Newtonian-like), 1/fm (~200 MeV, pion-like),
# 5/fm (~1 GeV, hadronic).
m_phi_values = [0.0, 0.5, 1.0, 2.0]  # in 1/fm

# λ = 1 for now (we'll scale results to compare to gravity).
lam = 1.0

# Use a radial BVP solver.  Define y = (φ, dφ/dr).
def make_ode(m_phi):
    def ode(r, y):
        phi, dphi_dr = y[0], y[1]
        # d²φ/dr² = m²(φ - δ) - λS(r) - (2/r) dφ/dr
        # at r ≈ 0 use limiting form: d²φ/dr² = m²(φ-δ)/3 - λS(0)/3 (factor from regularity)
        # Avoid divide-by-zero at r=0 by clipping.
        r_safe = np.maximum(r, 1e-8)
        d2phi_dr2 = m_phi**2 * (phi - delta) - lam * source(r) - 2.0 * dphi_dr / r_safe
        return np.vstack([dphi_dr, d2phi_dr2])
    return ode

def bc(ya, yb):
    # φ'(0) = 0, φ(L) = δ
    return np.array([ya[1], yb[0] - delta])

r_grid = np.linspace(0.01, L, 400)
phi_solutions = {}

for m_phi in m_phi_values:
    # Initial guess: φ = δ + a small bump at the center
    y_init = np.zeros((2, len(r_grid)))
    y_init[0, :] = delta - 0.1 * source(r_grid)

    sol = solve_bvp(make_ode(m_phi), bc, r_grid, y_init, tol=1e-6, max_nodes=10000)
    if not sol.success:
        print(f"  WARNING: BVP solve failed for m_φ = {m_phi}")
        continue
    phi_solutions[m_phi] = sol

    # Sample at a few characteristic radii
    r_samples = np.array([0.01, R_p, 2*R_p, 5*R_p, 10*R_p, L])
    phi_samples = sol.sol(r_samples)[0]
    delta_phi_samples = phi_samples - delta

    print(f"\n  m_φ = {m_phi:.2f} /fm ({m_phi * 197.327:.1f} MeV equivalent)")
    print(f"    r [fm]  |  δφ(r)  =  φ(r) - δ_eq")
    for rs, dps in zip(r_samples, delta_phi_samples):
        print(f"    {rs:6.2f}  | {dps:+.6e}")

# ----------------------------------------------------------------------
# Step 4.  Predict lepton mass shifts at the proton center.
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 4: Lepton mass shifts at the proton center (r → 0)")
print("=" * 70)

for m_phi, sol in phi_solutions.items():
    phi_at_zero = sol.sol(0.01)[0]  # near origin
    delta_phi = phi_at_zero - delta
    print(f"\n  m_φ = {m_phi:.2f} /fm:  δφ(0) = {delta_phi:+.6e} rad  (with λ=1)")
    for i, name in enumerate(lepton_names):
        d_ln_m = response[i] * delta_phi
        print(f"    δ(ln m_{name})/Δλ = {d_ln_m:+.3e}  →  m_{name} shift = {d_ln_m*100:+.5f}%")

# ----------------------------------------------------------------------
# Step 5.  Compare to experimental bounds.
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 5: Compare predicted shifts to experimental bounds")
print("=" * 70)
print("""
Experimental sensitivity scales:
  • Atomic clock comparison of optical transitions (Sr, Yb, Al+, ...):
    Δν/ν sensitivity ~ 10⁻¹⁸ over a year, can probe Δ(m_e/m_p) and Δα
    at the same level over Earth's Φ_grav/c² ~ 6.96 × 10⁻¹⁰ (Sun) or
    ~3 × 10⁻¹⁰ (annual variation in Earth-Sun distance).
  • Equivalence-principle (EP) violation bounds (MICROSCOPE, lunar laser):
    differential acceleration η ~ 10⁻¹⁵ — limits flavor-dependent forces.
  • Direct laboratory bounds on Δα in gravitational potential: 10⁻⁷ /(Φ/c²).
""")

# The v59 prediction was δα/α ~ 10⁻¹⁰ at Earth's surface (Φ/c² ~ 10⁻⁹).
# If φ has the same coupling as α to the gravitational potential, then
# δφ in Earth's well is also ~10⁻¹⁰ rad.

# Earth's surface gravitational potential (relative to infinity):
# Φ_Earth/c² ≈ -G M_Earth / (R_Earth · c²) ≈ -6.96 × 10⁻¹⁰
Phi_Earth_over_c2 = 6.96e-10

# If the Cosserat coupling sends Φ/c² → δφ at order unity:
delta_phi_earth = Phi_Earth_over_c2  # assumes coupling constant ~1

print(f"  At Earth's surface  Φ/c² = {Phi_Earth_over_c2:.3e}")
print(f"  If δφ ~ Φ/c² (order-unity coupling):")
print(f"    δφ ≈ {delta_phi_earth:.3e} rad")
print()
print(f"  → Predicted lepton mass shifts at Earth's surface:")
print(f"    flavor | δ(ln m) | δ(m_lepton/m_proton)")
m_proton = 938.272  # MeV
for i, name in enumerate(lepton_names):
    d_ln_m = response[i] * delta_phi_earth
    # δ(m_lepton/m_proton) ~ δ(ln m_lepton) - δ(ln m_proton)
    # assume δ(ln m_proton) = 0 (proton is hadronic, not from Brannen kernel)
    d_ratio = d_ln_m
    print(f"    {name:>4} | {d_ln_m:+.3e} | {d_ratio:+.3e}")

# Atomic clock bounds:
# Best published bound on Δ(m_e/m_p) in Earth's gravity: ~10⁻⁷ per Φ/c²
# (e.g., Lange et al 2021 Yb+ clock comparison).
# So δ(m_e/m_p) bound at Earth's surface ≈ 10⁻⁷ · 10⁻⁹ ~ 10⁻¹⁶.
print()
print(f"  Atomic-clock bound on δ(m_e/m_p) at Earth's surface ≈ 10⁻¹⁶")
print(f"  Our prediction (order-unity Cosserat):                  {response[1] * delta_phi_earth:+.3e}")
print()
print(f"  Ratio (prediction / bound):     {abs(response[1] * delta_phi_earth) / 1e-16:.2e}")

# ----------------------------------------------------------------------
# Step 6.  Summary.
# ----------------------------------------------------------------------
print()
print("=" * 70)
print("Step 6: Assessment — is the simple Cosserat coupling viable?")
print("=" * 70)

ratio = abs(response[1] * delta_phi_earth) / 1e-16
print(f"""
  • The Brannen response coefficient for the electron is large
    ({response[1]:.2f} per rad of δφ), because √m_e is much smaller than the
    Brannen scale a, so a small phase shift moves m_e by a large relative
    amount.

  • If the Cosserat coupling has order-unity strength to the gravitational
    potential (same as for α at the 10⁻¹⁰ level), the predicted shift
    Δ(m_e/m_p) ~ {response[1]*delta_phi_earth:+.3e} at Earth's surface.

  • This is {ratio:.1e}× the atomic-clock sensitivity bound of ~10⁻¹⁶ on
    Δ(m_e/m_p) per unit Φ_grav/c².  → Simple Cosserat coupling at this
    magnitude is RULED OUT by current data.

  • Two routes to revive the hypothesis:
     (a) The coupling strength to the lepton-mass sector is much WEAKER
         than to α: the Cosserat strain affects α via the EM-bivector
         grade but not the spectral-kernel grade.  Then δφ_kernel would
         be much smaller than δα/α.
     (b) The Brannen identification m_k = s_k² is wrong; perhaps
         instead m_k = a (1 + 2t cos α_k) directly (eigenvalues are
         masses, not √masses), in which case the response coefficients
         are radically different (the electron is no longer hypersensitive).

  Next experiment: redo Step 2 with the alternative identification
  m_k = s_k (linear, not quadratic) and see if the responses look
  more compatible with bounds.
""")

# Verify ordering and save final summary.
np.savez('/home/d/code/scp/v59/cosserat_experiment/01_predictions.npz',
         m_phi_values=np.array(m_phi_values),
         delta_phi_earth=delta_phi_earth,
         response=response,
         predicted_dln_m=response * delta_phi_earth,
         ratio_to_bound=ratio,
         Phi_Earth_over_c2=Phi_Earth_over_c2)

print("Saved fit and predictions to .npz files.")
