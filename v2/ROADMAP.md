# CHPT Numerical Program — Roadmap

**Status**: B=1–4 Skyrmions solved via rational map ansatz. Binding energies computed.

## Next Steps (Priority Order)

### 1. ~~Higher-B Solitons via Rational Map Ansatz~~ ✓ DONE
B=1–4 solved. See results below and `src/rational_map.c`.

### 2. Finite-λ Effects
Extend the radial solver beyond the sigma model limit (λ→∞). Include the bulk potential E_V in the hedgehog ODE and shooting method. Study how finite λ modifies the soliton profile, mass formula (Mc² = 2E₄ - 2E_V), and size. Map out the (λ, e) parameter space.

### 3. Degenerate Sector Perturbations
Linearize the degenerate sector (P, J) around the B=1 soliton background. Solve for bound states and scattering modes. Test whether parity-violating interactions emerge (weak force connection). This is where CHPT diverges from the standard Skyrme model.

### 4. 3D Field Reconstruction
Initialize a 3D grid with the known radial profile (properly resolved). Run full 8-component gradient flow to verify 3D energy matches 1D prediction. Validates the pipeline for scattering simulations.

### 5. Soliton-Soliton Scattering
Time-dependent 3D simulation of two solitons colliding. Observe repulsion, attraction, elastic scattering, or soliton production.

### 6. Sourced Maxwell Equations
Extract the electromagnetic field produced by a soliton. Connect topological charge Q to the Coulomb field. Derive the source current J from the nonlinear knot solution.

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

## Code Location

`proposal/hopfion_search/` — C (gcc + OpenMP) + Python visualization
