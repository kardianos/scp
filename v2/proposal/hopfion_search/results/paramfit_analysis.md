# Parameter Fitting: Skyrme Model to Proton

## Summary

Two physical inputs (proton mass M_p = 938.272 MeV, charge radius r_p = 0.841 fm)
constrain the Skyrme model parameters. The constraint leaves one free parameter (e).

## Numerical profile constants (B=1 hedgehog, sigma model)

From the radial solver (verified at multiple e, ρ₀ values):

| Quantity | Value | Definition |
|----------|-------|------------|
| K | 103.13 | E_sol × e / ρ₀³ (universal) |
| C_R | 1.4961 | R_rms × e / ρ₀ (universal) |
| C̃ | 1.058 | R_rms / √c₄ (rescaled RMS radius) |
| E/E_FB | 1.2315 | Ratio to Faddeev-Bogomolny bound |

Scaling laws: E_sol = 103.13 × ρ₀³/e,  R_rms = 1.496 × ρ₀/e

## Physical constraint

E_sol × R_rms = M_p × r_p / ℏc = 4.00 (dimensionless)

This fixes: **ρ₀⁴/e² = 0.02591**

Remaining freedom parameterized by e (the Skyrme coupling).

## Standard Skyrme convention (E₄ ∝ 1/e², no F_π)

In the standard convention where the Skyrme term coefficient is 1/(2e²) (independent of F_π):

| Quantity | Value |
|----------|-------|
| **e** | **3.69** |
| **F_π** | **95.0 MeV** |
| F_π (experimental, reduced) | 92.1 MeV |
| Agreement | 3.1% |

This gives excellent agreement with the experimental pion decay constant.

## Our code convention (E₄ ∝ ρ₀⁴/e²)

Our code includes ρ₀⁴ in the Skyrme term (from expanding |[∂q, ∂q]|² = ρ₀⁴|[∂q̂, ∂q̂]|²).
F_π = ρ₀ × Λ_E is constant at 56.5 MeV regardless of e — a factor ~3.3× below experimental.

This discrepancy reflects the non-standard normalization, NOT a physical disagreement. The
conversion factors (Λ_E, Λ_L) are convention-independent.

## Conversion factors (e=1, ρ₀=1 — scattering simulation parameters)

| Code unit | Physical value |
|-----------|---------------|
| 1 energy | 9.098 MeV |
| 1 length | 0.5624 fm |
| 1 time | 1.875 × 10⁻²⁴ s |
| Soliton mass | 103.13 code = 938.3 MeV |
| Soliton radius | 1.495 code = 0.841 fm |

## Particle masses in code units (at e=1, ρ₀=1)

| Particle | Mass (MeV) | Code units | Accessible at v=0.5c? |
|----------|-----------|------------|----------------------|
| pion | 139.6 | 15.3 | YES |
| rho(770) | 775.3 | 85.2 | YES |
| proton | 938.3 | 103.1 | — (rest mass) |
| Delta(1232) | 1232 | 135.4 | YES |
| W± | 80,369 | 8,834 | NO (need γ≈86) |
| Z⁰ | 91,188 | 10,023 | NO (need γ≈97) |
| Higgs | 125,200 | 13,761 | NO (need γ≈133) |

At v=0.5c, CM collision energy ≈ 237 code units ≈ 2160 MeV.

## Key finding

**W/Z/Higgs are inaccessible** in our current simulations — they require ultra-relativistic
collisions (γ ≈ 86–133), far beyond what the lattice can handle.

**Pion, rho, and Delta features ARE accessible** at current energies and should be the
primary targets for spectral analysis of the scattering simulations.

## Scattering results in physical units

### B+B repulsive (v=0.5c)
- Duration: 2.35 code = 4.4 × 10⁻²⁴ s
- Closest approach: 1.23 code = 0.69 fm
- Deceleration: 18× → collision time ~ 2 × 10⁻²⁴ s

### B+B̄ annihilation (v=0.5c)
- Duration: 5.0 code = 9.4 × 10⁻²⁴ s
- Energy radiated: 104 code = 946 MeV = 1.01 M_p
- Post-collision equipartition: E_pot ≈ E_kin ≈ 113 code = 1028 MeV

### Timescale comparison
- W/Z lifetime: ~3 × 10⁻²⁵ s = 0.16 code time units
- rho lifetime: ~4.5 × 10⁻²⁴ s = 2.4 code time units
- Delta lifetime: ~5.6 × 10⁻²⁴ s = 3.0 code time units
- Collision time: ~4.4 × 10⁻²⁴ s = 2.35 code time units

The collision timescale (~2.4 code) is comparable to the rho and Delta lifetimes,
suggesting these resonances should leave observable signatures in the scattering field.
