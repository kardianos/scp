#!/usr/bin/env python3
"""Step 3: Compute virial profile integrals from V45 SFA data.

Uses the D=80 baseline (non-interacting baryons) to extract single-baryon
profile integrals needed for the η₁ self-consistency equation.

Reads the diag.tsv for global energetics, then computes the virial
residual R and the topology-weighted integral J₁.

η₁ = -R / (2·η₀·J₁)

Output: step3_integrals_output.txt
"""
import csv
import numpy as np
import sys

# Parameters
m2 = 2.25
mu = -41.345
kappa = 50.0
eta0 = 0.5
A_bg = 0.1

def load_diag(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}

# Load D=80 baseline (two non-interacting baryons)
d80 = load_diag('/home/d/code/scp/v45/v45_D80_diag.tsv')
mask = (d80['t'] >= 200) & (d80['t'] <= 400)

# Time-averaged energy components
E_kin = d80['E_phi_kin'][mask].mean()      # ½∫(∂φ/∂t)² = kinetic
E_grad = d80['E_grad'][mask].mean()         # ½∫|∇φ|²
E_mass = d80['E_mass'][mask].mean()         # ½∫m²φ²
E_pot = d80['E_pot'][mask].mean()           # ∫V(P)
E_tkin = d80['E_theta_kin'][mask].mean()    # ½∫(∂θ/∂t)²
E_tgrad = d80['E_tgrad'][mask].mean()       # ½∫|∇θ|²
E_coupling = d80['E_coupling'][mask].mean() # η₀∫φ·(∇×θ)
E_total = d80['E_total'][mask].mean()

print("="*70)
print("STEP 3: NUMERICAL PROFILE INTEGRALS FROM V45 D=80")
print("="*70)
print()
print("Time-averaged energetics (t=200-400):")
print(f"  E_phi_kin   = {E_kin:.1f}")
print(f"  E_grad      = {E_grad:.1f}")
print(f"  E_mass      = {E_mass:.1f}")
print(f"  E_pot       = {E_pot:.1f}")
print(f"  E_theta_kin = {E_tkin:.1f}")
print(f"  E_tgrad     = {E_tgrad:.1f}")
print(f"  E_coupling  = {E_coupling:.1f}")
print(f"  E_total     = {E_total:.1f}")
print()

# Background energy estimates (to subtract)
N = 512; L = 100.0
dx = 2*L/(N-1)
dV = dx**3
N3 = N**3
k_bg = np.pi/L
omega_bg = np.sqrt(k_bg**2 + m2)

# Background contributes to E_kin, E_mass, E_grad but NOT to E_pot, E_coupling
E_kin_bg = 0.5 * 3 * (omega_bg**2 * A_bg**2 * 0.5) * N3 * dV
E_mass_bg = 0.5 * m2 * 3 * (A_bg**2 * 0.5) * N3 * dV
E_grad_bg = 0.5 * 3 * (k_bg**2 * A_bg**2 * 0.5) * N3 * dV

print("Background energy estimates:")
print(f"  E_kin_bg    = {E_kin_bg:.1f}")
print(f"  E_mass_bg   = {E_mass_bg:.1f}")
print(f"  E_grad_bg   = {E_grad_bg:.1f}")
print()

# Soliton-only energies (perturbation above background)
# Note: these can be negative (depletion zones)
E_kin_sol = E_kin - E_kin_bg
E_mass_sol = E_mass - E_mass_bg
E_grad_sol = E_grad - E_grad_bg

print("Soliton-only energies (total minus background):")
print(f"  E_kin_sol   = {E_kin_sol:.1f}")
print(f"  E_mass_sol  = {E_mass_sol:.1f}")
print(f"  E_grad_sol  = {E_grad_sol:.1f}")
print(f"  E_pot_sol   = {E_pot:.1f}  (all from soliton)")
print(f"  E_tgrad_sol = {E_tgrad:.1f}  (all from soliton)")
print(f"  E_coup_sol  = {E_coupling:.1f}  (all from soliton)")
print()

# ================================================================
# Virial residual R
# ================================================================
# The virial identity: E₁ + 2·E₂ + 3·E₃ = 0
# E₁ = ½I_grad + ½I_tgrad (λ¹ terms)
# E₂ = E_coupling (λ² terms)
# E₃ = ½I_mass + I_pot (λ³ terms)
#
# Using the TOTAL energies (not soliton-only) because the virial
# applies to the whole field configuration:

E1 = E_grad + E_tgrad  # Already multiplied by ½ in the diagnostics
E2 = E_coupling
E3 = E_mass + E_pot    # Already ½I_mass in diagnostics

R_total = E1 + 2*E2 + 3*E3

print("Virial analysis (TOTAL field, including background):")
print(f"  E₁ (λ¹ terms) = E_grad + E_tgrad = {E1:.1f}")
print(f"  E₂ (λ² terms) = E_coupling = {E2:.1f}")
print(f"  E₃ (λ³ terms) = E_mass + E_pot = {E3:.1f}")
print(f"  R = E₁ + 2E₂ + 3E₃ = {R_total:.1f}")
print()

# The background itself satisfies its own virial (it's a plane wave solution).
# Background virial:
# E1_bg = E_grad_bg (no θ background)
# E2_bg = 0 (no coupling in background)
# E3_bg = E_mass_bg (no V(P) in background since P_bg≈0)
R_bg = E_grad_bg + 0 + 3*E_mass_bg

print("Background virial (should be ~0 for plane wave):")
print(f"  R_bg = E_grad_bg + 3·E_mass_bg = {R_bg:.1f}")
print()

# Soliton virial residual (subtract background virial)
R_sol = R_total - R_bg

print("Soliton virial residual:")
print(f"  R_sol = R_total - R_bg = {R_sol:.1f}")
print()

# ================================================================
# Estimate J₁ = ∫|P|·φ·∇×[G*(∇×φ)] d³x
# ================================================================
# At leading order (θ = η₀ × G*(∇×φ)):
# I_C1 = ∫|P|·φ·(∇×θ) d³x ≈ η₀ × J₁
# We don't have I_C1 directly from diagnostics.
# But we can estimate it from the coupling energy and P_int:
#
# I_C0 = ∫φ·(∇×θ) d³x = E_coupling / η₀
# I_C1 = ∫|P|·φ·(∇×θ) d³x ≈ ⟨|P|⟩_coupling × I_C0
#
# where ⟨|P|⟩_coupling is the coupling-energy-weighted average of |P|.
# From V45 data: P_int ≈ 642, and the core has |P|≈0.08-0.1.
# The coupling energy is concentrated at the core, so ⟨|P|⟩ ≈ 0.05-0.1.

I_C0 = E_coupling / eta0
print("Coupling integrals:")
print(f"  I_C0 = E_coupling/η₀ = {I_C0:.1f}")
print()

# Estimate I_C1 for different assumed ⟨|P|⟩ values
print("Estimated I_C1 = ⟨|P|⟩ × I_C0:")
for P_avg in [0.01, 0.03, 0.05, 0.08, 0.1]:
    I_C1_est = P_avg * I_C0
    J1_est = I_C1_est / eta0
    eta1_est = -R_sol / (2 * eta0 * J1_est) if J1_est != 0 else float('inf')
    print(f"  ⟨|P|⟩ = {P_avg:.2f}: I_C1 = {I_C1_est:.2f}, J₁ = {J1_est:.2f}, "
          f"η₁ = {eta1_est:.1f}")

print()
print("="*70)
print("CROSS-CHECK: η₁ = |μ|/η₀")
print("="*70)
print()
eta1_conj = abs(mu) / eta0
print(f"  |μ|/η₀ = {eta1_conj:.1f}")
print()

# What ⟨|P|⟩ would make η₁ = |μ|/η₀?
# η₁ = -R_sol / (2η₀²·⟨|P|⟩·I_C0/η₀) = -R_sol / (2η₀·⟨|P|⟩·I_C0)
# |μ|/η₀ = -R_sol / (2η₀·⟨|P|⟩·I_C0)
# ⟨|P|⟩ = -R_sol / (2η₀·I_C0·|μ|/η₀) = -R_sol / (2|μ|·I_C0)
if I_C0 != 0:
    P_avg_needed = -R_sol / (2 * abs(mu) * I_C0)
    print(f"  For η₁ = |μ|/η₀, need ⟨|P|⟩_coupling = {P_avg_needed:.4f}")
    print(f"  This is {'REASONABLE' if 0.01 < P_avg_needed < 0.2 else 'UNREASONABLE'}")
    print(f"  (Expected range: 0.01-0.1 based on braid core P values)")

print()
print("="*70)
print("ALTERNATIVE: Direct force balance estimate")
print("="*70)
print()
# From step1: η₁ ≈ (1/η₀ - η₀)/A_core³
for A_core in [0.15, 0.2, 0.25, 0.3, 0.35]:
    eta1_fb = (1/eta0 - eta0) / A_core**3
    print(f"  A_core = {A_core}: η₁ = {eta1_fb:.1f}, "
          f"η_eff(|P|=0.1) = {eta0 + eta1_fb*0.1:.1f}")

print()
print("="*70)
print("SUMMARY")
print("="*70)
print()
print(f"Virial residual R_sol = {R_sol:.1f}")
print(f"Coupling integral I_C0 = {I_C0:.1f}")
print(f"Conjectured η₁ = |μ|/η₀ = {eta1_conj:.1f}")
if I_C0 != 0:
    print(f"Required ⟨|P|⟩ for conjecture = {P_avg_needed:.4f}")
print(f"Force balance η₁ range: 35-450 (A_core = 0.35-0.15)")
print(f"|μ|/η₀ = {eta1_conj:.1f} falls in the force balance range: YES")
