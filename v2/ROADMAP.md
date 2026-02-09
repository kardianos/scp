# CHPT Numerical Program — Roadmap

**Status**: Phase 8 — Soliton dynamics, scattering, and extended Lagrangian analysis. All prior phases complete.

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

### 4. ~~3D Field Reconstruction~~ ✓ DONE
Initialize a 3D grid from 1D radial profiles. Verify 3D energy matches 1D.

All B=1–4 initialized and verified to <0.1% accuracy (N=256). Gradient flow loses topology at all N — Skyrmion is a saddle point on the lattice. 1D radial solver is authoritative. Full details in `results/README.md`.

### 5. ~~Soliton-Soliton Scattering~~ ✓ COLLISION OBSERVED

**Code**: `src/scatter.c` — leapfrog (Störmer-Verlet) integrator with product ansatz initialization, Lorentz-boosted initial velocity, Dirichlet boundary clamping, 3D centroid tracking.

**Successful collision parameters**: e=1, N=192, L=10, σ-model profile, v=0.5c, z₀=1.5 (sep=3), λ=5000, dt=0.001.

**Repulsive channel (B+B, π-isorotation around ê₁)**:
- Q=1.9999 maintained for 2.35 time units (entire approach + deep interpenetration)
- 3D centroid tracking: solitons approach along z, decelerate 18× (0.74/t → 0.04/t)
- Point of closest approach: r≈1.23 (0.87 core radii apart)
- No transverse (90°) scattering — purely axial dynamics
- Bounce extrapolated at t≈3.2, but lattice topology loss starts at t≈2.4

**Key technical insights**:
- Core resolution is the critical parameter: ≥13 grid pts across core radius for ~2.4t stability
- σ-model profile eliminates product-ansatz breathing (|q₁·q₂/ρ₀|=ρ₀ exactly)
- 3D charge-weighted centroid essential (z-axis minimum tracker gives false "merger" picture)
- Smaller e = wider core = better lattice stability (e=1 vs e=2 is 3× improvement)
- All damping strategies (velocity, radial, settling) harm topology — conservative evolution only

**Dead ends**: finite-λ profile (breathing → crash t≈0.5), radial damping (shrinks core → crash), e=2 at any N≤256 (insufficient core resolution).

### 6. ~~Sourced Maxwell Equations~~ ✓ DONE (B=1)

**Code**: `src/maxwell.c` — computes EM field from soliton charge current.

B=1 hedgehog: Q=1.000 (exact), radial Coulomb field E_r = Q/(4πr²), E_EM = 0.28% of soliton mass. B>1 has known bug (hedgehog formula without rational map angular factors, gives Q=1 for all B).

### 7. ~~Extended Lagrangian (Degenerate Coupling)~~ ✓ ANALYSIS COMPLETE

Three coupling options analyzed theoretically and numerically:

| | Option 1 (component-wise L₂) | Option 2 (g²\|q\|²\|∇p\|²) | Option 3 (full Skyrme norm) |
|---|---|---|---|
| Propagation | ✓ | ✓ | ✓ |
| Soliton coupling | None | Attractive well | Repulsive barrier |
| Bound states | No | Yes (at finite λ) | No |

**Key result**: Modified Skyrme term (Option 3) gives **repulsive** (positive-definite) effective potential V_eff(r) for degenerate modes — no bound states. Option 2 gives **attractive** well at finite λ where ρ(r)<ρ₀.

**Recommendation**: Hybrid approach (Options 1+2+3) gives Lennard-Jones-like potential: short-range repulsion (Skyrme) + medium-range attraction (bulk coupling). Full Lagrangian: L₂ + L₂,D + L₄,full + g²|q|²|∇p|²/2 - V - V_D.

See `results/extended_lagrangian_analysis.md` for derivation and numerical V_eff data.

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
- [x] 3D field reconstruction verified (verify3d.c)
  - All B=1–4 to <0.1% accuracy at N=256
  - 4th-order convergence in h (matches stencil order)
  - Gradient flow loses topology at all N (lattice saddle point)
  - B=3 rational map bug found and fixed (z³→z² denominator)
- [x] Sourced Maxwell equations (maxwell.c)
  - B=1: Q=1.000, Coulomb field, E_EM = 0.28% of soliton mass
- [x] Extended Lagrangian analysis (veff.c + analysis document)
  - Three coupling options analyzed (component-wise, bulk coupling, full Skyrme)
  - Modified Skyrme term: repulsive V_eff (positive-definite)
  - Bulk coupling g²|q|²|∇p|²: attractive well at finite λ
  - Hybrid approach recommended (Lennard-Jones-like potential)
- [x] Soliton scattering (scatter.c) — COLLISION OBSERVED
  - Leapfrog integrator, product ansatz, Lorentz boost, 3D centroid tracking
  - Repulsive B+B at v=0.5c: Q=1.9999 for 2.35t, solitons decelerate 18×
  - Point of closest approach r≈1.23 (0.87 core radii)
  - Key: e=1 σ-model profile, ≥13 grid pts across core, no damping
  - Attractive channel: needs more resolution (crashes at t≈1.7)

## Phase 9 — Proposed Next Steps (Priority Order)

### 9.1 ~~B+B̄ Annihilation Simulation~~ ✓ DONE — INELASTIC PASS-THROUGH
At v=0.5c, B+B̄ does NOT annihilate — the pair passes through each other. Attractive interaction accelerates approach (E_kin: 9.5→153, 16× increase), 66% of E_pot converted to E_kin at closest approach (r≈1.45), then pair separates. 48% of rest mass energy permanently radiated (E_pot: 217→113). Late-time equipartition E_pot≈E_kin≈113. Q=0.0000 maintained for full 5.0 time units. Complete annihilation would require lower velocity. See `results/README.md` and `data/scatter_anti_v0.50.dat`.

### 9.2 ~~Observe the Bounce~~ ✓ DONE — BOUNCE CONFIRMED
N=256, L=10, e=1 (h=0.078, 18.1 pts across core). Bounce turning point at **r=1.222** (t≈2.70) when Q=1.92 (96% preserved). Consistent with N=192 extrapolation. Post-bounce separation visible (r: 1.222→1.248) but polluted by lattice topology crash (Q→1.2 by t=2.9). The bounce is physical — the repulsive interaction genuinely reverses the approach at r≈0.86 core radii. See `results/README.md`.

### 9.3 Implement Degenerate Sector Dynamics (HIGH IMPACT, HIGH EFFORT)
Implement the hybrid extended Lagrangian (L₂,D + g²|q|²|∇p|² + full Skyrme norm) in the 3D time evolution code. Requires adding 4 more field components (j1,j2,j3,p) and computing bound states numerically. Addresses the biggest theoretical gap (B6/B3).

### 9.4 Visualization (MEDIUM IMPACT, LOW EFFORT)
Publication-quality plots: scattering trajectories, energy density isosurfaces during collision, radial profile comparisons, V_eff potential well for degenerate modes.

### 9.5 Spec Chapter Drafting (for proposal)
Turn numerical results into proposal narrative, particularly the soliton dynamics/scattering chapter.

## Code Inventory

`proposal/hopfion_search/` — C (gcc + OpenMP) + Python visualization. `make` builds all into `bin/`.

| Directory | Contents |
|-----------|----------|
| `src/` | C source files and headers |
| `bin/` | Compiled executables |
| `data/` | Generated profiles, Maxwell fields, scatter output |
| `scripts/` | Shell batch scripts, Python visualization |
| `results/` | Analysis documents, archived run logs, figures |

| Source | Binary | Purpose |
|--------|--------|---------|
| `src/radial.c` | `bin/radial_solver` | B=1 shooting method solver |
| `src/rational_map.c` | `bin/rational_map_solver` | B=1–4 rational map ansatz solver |
| `src/finite_lambda.c` | `bin/finite_lambda_solver` | Coupled f-shooting + ρ-BVP Newton |
| `src/scatter.c` | `bin/scatter` | Soliton-soliton scattering simulator |
| `src/verify3d.c` | `bin/verify3d` | 3D initialization from 1D profiles |
| `src/maxwell.c` | `bin/maxwell` | Sourced Maxwell equations |
| `src/veff.c` | `bin/veff` | Effective potential for degenerate modes |
| `src/field.c` | (library) | 3D energy functional & gradient (4th-order) |
| `src/verify.c` | `bin/verify_gradient` | Gradient verification |
| `src/main.c` | `bin/soliton_search` | 3D gradient flow (historical) |
