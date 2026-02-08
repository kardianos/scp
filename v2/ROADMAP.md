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

### 4. 3D Field Reconstruction — IN PROGRESS
Initialize a 3D grid from 1D radial profiles. Verify 3D energy matches 1D. Run gradient flow to confirm stationarity and topology preservation.

**Status (Phase 8)**:
- `src/verify3d.c` written: loads profile_B{1,2,3,4}.dat, initializes 3D hedgehog or rational map grid, computes energy/Q, optionally runs sigma-model gradient flow
- Batch script: `run_verify3d.sh`

**Initialization accuracy** (VERIFIED — no relaxation needed):

| N | L | h | E error | Q error | E₂/E₄ | Core pts |
|---|---|---|---------|---------|--------|----------|
| 64 | 6 | 0.188 | -10.1% | -0.103 | 1.154 | 1.9 |
| 96 | 6 | 0.125 | -2.3% | -0.025 | 1.040 | 2.8 |
| 128 | 6 | 0.094 | -0.5% | -0.008 | 1.020 | 3.8 |
| 160 | 6 | 0.075 | +0.1% | -0.004 | 1.016 | 4.7 |
| 192 | 6 | 0.063 | +0.4% | -0.002 | 1.015 | 5.7 |
| 256 | 8 | 0.063 | +0.1% | -0.002 | 1.009 | 5.7 |

Soliton core radius: √c₄ ≈ 0.354 (e=4). "Core pts" = core radius / h.

**Topology loss during gradient flow** (CONCLUDED — topology lost at all N):

| N | h | Core pts | Steps to Q<0.5 | Final Q |
|---|---|----------|----------------|---------|
| 128 | 0.094 | 3.8 | ~8 | 6e-8 |
| 160 | 0.075 | 4.7 | ~22 | 3e-8 |
| 192 | 0.063 | 5.7 | ~35 | 5e-6 |

Root cause: Skyrmion is a saddle point of the unconstrained lattice energy. Continuous topology is only approximately conserved on a discrete grid. Higher N delays but cannot prevent unwinding.

**Higher-B initialization** (ALL B=1–4 VERIFIED):

| B | N | L | E error | Q error | E₂/E₄ |
|---|---|---|---------|---------|-------|
| 1 | 256 | 6 | +0.74% | -0.001 | 1.018 |
| 2 | 256 | 6 | -0.02% | -0.001 | 1.002 |
| 3 | 256 | 6 | -0.04% | -0.002 | 1.001 |
| 4 | 256 | 6 | -0.08% | -0.001 | 1.000 |

B=3 bug found and fixed: verify3d.c had wrong rational map denominator (z³ instead of z²). After fix, B=3 converges normally (0.005% at N=512 L=8). See `results/README.md` for full details.

**Archived results**: `results/` directory with run logs and `README.md`.

**Status**: 3D field reconstruction COMPLETE. All B=1–4 initialized and verified. Gradient flow loses topology at all N (lattice saddle point). 1D radial solver is authoritative.

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
