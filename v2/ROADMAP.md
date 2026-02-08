# CHPT Numerical Program — Roadmap

**Status**: B=1–4 Skyrmions solved. Finite-λ effects mapped. Degenerate sector decoupling proven.

## Next Steps (Priority Order)

### 1. ~~Higher-B Solitons via Rational Map Ansatz~~ ✓ DONE
B=1–4 solved. See results below and `src/rational_map.c`.

### 2. ~~Finite-λ Effects~~ ✓ DONE
Self-consistent solutions from λ=10⁸ to λ=8000 via coupled f-shooting + ρ-BVP Newton with under-relaxation. Key findings: 28% mass reduction at λ=8000, FB bound violated at λ≈9000, collapse instability below λ≈7000–8000. σ-model soliton is a LOCAL energy minimum; global minimum = collapse (ρ→0). Mass formula: Mc²=2E₄-2E_V. See `src/finite_lambda.c` and `spec/math/05_mass_mechanism.md` §7.

### 3. ~~Degenerate Sector Analysis~~ ✓ DONE (Theoretical Result)
**Critical finding**: The degenerate sector (P, J) is **completely decoupled** from the soliton in the current Lagrangian. The scalar extraction <...>₀ in Cl(3,0,1) kills all e₀ terms, so E₂, E₄, and E_V depend only on the bulk quaternion q. Moreover, <(∂Ψ)(∂Ψ̃)>₀ = |∂q|² — the degenerate sector has **no kinetic energy at all** (neither gradient nor temporal). The EOM is algebraic: μ²J=0, μ²P=0 → forced to zero. The modes are non-dynamical. See `spec/15_open_problems.md` B3/B5.

**Implication**: To make the degenerate sector physical, the theory needs:
- A **degenerate kinetic term**: L₂,D = (1/2c²)|∂_t p|² - (1/2)|∇p|²
- A **coupling term** between bulk and degenerate sectors (e.g., through a modified Skyrme term or explicit interaction)

### 4. 3D Field Reconstruction
Initialize a 3D grid with the known radial profile (properly resolved). Run full 8-component gradient flow to verify 3D energy matches 1D prediction. Validates the pipeline for scattering simulations.

### 5. Soliton-Soliton Scattering
Time-dependent 3D simulation of two solitons colliding. Observe repulsion, attraction, elastic scattering, or soliton production.

### 6. Sourced Maxwell Equations
Extract the electromagnetic field produced by a soliton. Connect topological charge Q to the Coulomb field. Derive the source current J from the nonlinear knot solution.

### 7. Extended Lagrangian (Degenerate Coupling)
If the degenerate sector is to have physical content (weak force connection), the Lagrangian must be extended. Options:
- Replace <∂Ψ(∂Ψ̃)>₀ with component-wise norm Σ_α(∂Ψ_α)² (gives propagation, no coupling)
- Add cross-coupling like g²|q|²|∇p|² (gives soliton-dependent effective mass for P, J)
- Modified Skyrme term using full component norm instead of <...>₀
Study requires theoretical analysis followed by numerical verification.

## Completed

- [x] Lagrangian and EOM established (5 parameters: ρ₀, λ, e, μ, c)
- [x] 3D energy functional and gradient implemented (field.c)
- [x] Skyrme force derived analytically and verified to ~1e-8 (verify.c)
- [x] Static decoupling theorem confirmed (bulk = Skyrme, degenerate → 0)
- [x] B=1 hedgehog Skyrmion found via shooting method (radial.c)
- [x] E/E_FB = 1.232 confirmed, virial E₂=E₄ verified, Q=1.000000
- [x] FB bound derived: E_FB = 6√2 π² ρ₀³/e
- [x] Profile visualization (viz_profile.py)
- [x] B=1–4 solitons via rational map ansatz (rational_map.c)
  - Angular integrals: I(1)=1.00, I(2)=5.81, I(3)=13.58, I(4)=20.65
  - Near-origin behavior: f(r) ~ π - ar^α, α=(-1+√(1+8B))/2
  - E/E_FB: B=1: 1.231, B=2: 1.208, B=3: 1.184, B=4: 1.137
  - Binding: E(B)/(B·E₁) = 0.981 (B=2), 0.962 (B=3), 0.923 (B=4)
  - Virial E₂=E₄ verified for all B, Q=B exact
- [x] Finite-λ effects mapped (finite_lambda.c)
  - Self-consistent solutions: λ=10⁸ to λ=8000
  - ρ(0) drops from 1.0 to 0.803; E/E_FB from 1.232 to 0.891
  - Virial E₂-E₄+3E_V=0 verified; Mass Mc²=2E₄-2E_V
  - Collapse instability below λ≈7000–8000 (global minimum = ρ→0)
  - ρ BVP exponentially stiff: must use Thomas algorithm, not shooting
- [x] Degenerate sector decoupling proven (algebraic + numerical)
  - <(∂Ψ)(∂Ψ̃)>₀ = |∂q|² — no degenerate kinetic energy
  - E₄ scalar extraction kills e₀ terms — no Skyrme coupling
  - Degenerate modes non-dynamical in current Lagrangian

## Code Location

`proposal/hopfion_search/` — C (gcc + OpenMP) + Python visualization
