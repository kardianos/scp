# CHPT Numerical Program — Roadmap

**Status**: Phase 14 complete — K=1 angular modes analyzed, no bound pion vibrational mode. Pion = quantized isorotational zero mode. Near-binding in mixed channel (4% margin). All bound-state searches (Phases 10–14) definitively negative.

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
- [x] Degenerate sector coupling (coupling.c, verify_coupling.c)
  - E_{2,D} (free propagation) + E_{4,C} (Skyrme cross-coupling) + E_int (bulk-degenerate)
  - Gradient verified to <1e-7; energy conservation ΔE/E < 2e-4
  - E_{4,C} repulsion confirmed: degenerate field expelled 94–99% from core
  - 8-component evolution integrated into scatter.c
- [x] Bound state eigenvalue search (bound.c) — NO BOUND STATES (Phase 10+11)
  - Shooting + bisection on ω², square-well verification to <0.1%
  - Skyrme coupling: V_eff ≥ μ²/2 everywhere, β²/(4K) correction negligible
  - Geometric coupling: cross-term ≡ 0 (Lie algebra orthogonality), no well
  - Momentum formulation (g, p=r²(2Kg'+βg)) avoids β' noise, grid-converged
  - Fixed: GL weight formula, profile format detection, β' numerical instability
- [x] Normal mode analysis (normal_modes.c) — NO BREATHING RESONANCE (Phase 12 Step 3)
  - K=0 breathing mode: no bound states in continuum, spectral peak is box mode
  - Moment of inertia: Λ=141.6 (massless), Λ=79.18 (massive)
  - Massive pion profile: no instabilities, no bound breathing modes below m_π² threshold
- [x] Pion mass in 3D (field.c, scatter.c) — CFL FIXED (Phase 13)
  - V_pi = m_pi² ρ₀² (1 - q·s/|q|), force tangent to S³
  - CFL for Skyrme term scales as h² not h; dt≤0.015 at h=0.083
  - Massive scattering qualitatively same as massless
- [x] K=1 angular modes (angular_modes.c) — NO BOUND STATES (Phase 14)
  - Channel A: zero mode (pion = quantized isorotation)
  - Channel B: no bound states (λ=0.484 > threshold 0.243)
  - Channel C: near-binding (W/m min=0.252, only 4% above threshold 0.243)
  - Bug fixed: pion mass term sin²f·cosf → cosf

## Phase 9 — Proposed Next Steps (Priority Order)

### 9.1 ~~B+B̄ Annihilation Simulation~~ ✓ DONE — INELASTIC PASS-THROUGH
At v=0.5c, B+B̄ does NOT annihilate — the pair passes through each other. Attractive interaction accelerates approach (E_kin: 9.5→153, 16× increase), 66% of E_pot converted to E_kin at closest approach (r≈1.45), then pair separates. 48% of rest mass energy permanently radiated (E_pot: 217→113). Late-time equipartition E_pot≈E_kin≈113. Q=0.0000 maintained for full 5.0 time units. Complete annihilation would require lower velocity. See `results/README.md` and `data/scatter_anti_v0.50.dat`.

### 9.2 ~~Observe the Bounce~~ ✓ DONE — BOUNCE CONFIRMED
N=256, L=10, e=1 (h=0.078, 18.1 pts across core). Bounce turning point at **r=1.222** (t≈2.70) when Q=1.92 (96% preserved). Consistent with N=192 extrapolation. Post-bounce separation visible (r: 1.222→1.248) but polluted by lattice topology crash (Q→1.2 by t=2.9). The bounce is physical — the repulsive interaction genuinely reverses the approach at r≈0.86 core radii. See `results/README.md`.

### 9.3 ~~Implement Degenerate Sector Dynamics~~ ✓ DONE — COUPLING IMPLEMENTED & TESTED
Three coupling terms implemented in `src/coupling.c` + `src/coupling.h`:
- **E_{2,D}**: degenerate gradient energy (free propagation)
- **E_{4,C}**: Skyrme cross-coupling (repulsive near soliton core, from full 8-component norm)
- **E_int**: bulk-degenerate coupling g²|q|²|∇w|² (attractive well at finite-λ)

Gradient verified to <1e-7 relative error (`src/verify_coupling.c`). Integration tests (e=1, λ=10000, N=128, L=8): energy conservation ΔE/E < 2e-4, E_{4,C} repulsion confirmed (degenerate field expelled 94–99% from core), E_int retains 14% more coupling energy (weak trapping). `scatter.c` updated with `-degenerate -g <val>` flags for 8-component evolution.

### 9.4 Visualization (MEDIUM IMPACT, LOW EFFORT)
Publication-quality plots: scattering trajectories, energy density isosurfaces during collision, radial profile comparisons, V_eff potential well for degenerate modes.

### 9.5 Spec Chapter Drafting (for proposal)
Turn numerical results into proposal narrative, particularly the soliton dynamics/scattering chapter.

---

## Phase 10 — Degenerate Bound State Search ✓ COMPLETE: NO BOUND STATES

**Code**: `src/bound.c` — standalone eigenvalue solver (shooting + bisection on ω²).

**Method**: For each radial point, compute angular-averaged Skyrme cross-coupling coefficients (K, W, β) via Gauss-Legendre quadrature over the hedgehog background, including E_{2,D}, E_int, and E_{4,C} contributions. Solve the radial Schrödinger-like eigenvalue equation g'' = f(K, W, β, ω²) g via RK4 shooting from r=δ to R_max, then scan ω² ∈ [0, c²μ²] for sign changes of g(R_max) and refine by bisection.

**Two channels tested**:
- **P-mode** (ℓ=0, pseudoscalar): g(0)=1, g'(0)=0
- **J-mode** (ℓ=1, hedgehog-like vector): g~r near origin

**Verification**: Square-well test matches analytical eigenvalues to <0.1% (ground state 0.013%, 3 eigenvalues found).

**Parameter scan**: 50 combinations of g_c ∈ {0.5, 1, 2, 5, 10} × μ ∈ {0.1, 0.5, 1, 2, 5} on two profiles:
- e=2, λ=5000: ρ(0)=0.969, 23% density deficit — **zero bound states everywhere**
- e=1, λ=8000: ρ(0)=0.996 — **zero bound states everywhere**

**Root cause**: E_{4,C} is positive-definite (squared norm |F^w|²), creating a massive repulsive kinetic barrier K(r) ~ 200 at the soliton core. Even where E_int creates an attractive well W_eff ~ -4, the effective Schrödinger potential W_eff/K ~ -0.02 is far too shallow to support bound states. The J-mode is even worse (centrifugal + Skyrme barrier). The problem is fundamental: any squared-norm coupling term is always repulsive.

**Implication**: The ad hoc E_{4,C} = (1/4e²)Σ|F^w_{ij}|² cannot produce weak boson candidates. A geometrically motivated coupling is needed → Phase 11.

---

## Phase 11 — Geometric Covariant Derivative Coupling ✓ COMPLETE: NO BOUND STATES

### Motivation

The null result of Phase 10 has a clear root cause: E_{4,C} is a squared norm, hence always repulsive. The geometrically natural alternative is the **covariant derivative** D_i p = ∂_i p + g_cov[A_i, p], where A_i = q̃₀ ∂_i q₀ is the SU(2) connection from the hedgehog.

### Implementation

Added `-geom` / `-skyrme` mode to `src/bound.c` with `compute_cov_coupling()` function. Computes angular-averaged |D_i p|² for both P-mode (singlet) and J-mode (adjoint).

### Bug fixes during implementation

1. **Gauss-Legendre weight formula**: w_i = 2(1-x²)/(n²P_{n-1}²), NOT 2/((1-x²)n²P_{n-1}²). Off by factor (1-x²)². Fixed in bound.c and veff.c. (rational_map.c was correct.)
2. **Profile format**: 6-column sigma-model profiles have baryon_dens as col 4, not ρ. read_profile now detects nc≥5 and uses ρ=1.
3. **β' numerical noise**: Old W_eff = W - β'/2 - β/r formula amplified noise from β'. Reformulated shooting as momentum system (g, p = r²(2Kg'+βg)) which avoids β' entirely. V_eff = W - β²/(4K) is now grid-converged.

### Results: DEFINITIVE NULL

**Geometric coupling (covariant derivative)**:
- **P-mode (singlet)**: [A_i, scalar] = 0, no coupling at all. Propagates freely.
- **J-mode (adjoint)**: Cross-term ⟨∂_i p, [A_i, p]⟩ = 0 identically due to **Lie algebra orthogonality**: [A, P₀] ⊥ P₀ for any pure imaginary quaternion P₀.
- β_cov ≡ 0 for hedgehog J-mode. The W correction is 4sin²f/r² × g_cov(g_cov-1):
  - g_cov = 1: W correction = 0 (neutral, connection absorbed into covariant derivative)
  - g_cov = 1/2: centrifugal barrier saturated but W_eff = cos²f/r² + μ²/2 ≥ μ²/2
- **No bound states** for any g_cov ∈ {0, 0.3, 0.5, 1, 2, 5, 10}.

**Skyrme coupling (E_{4,C}) with corrected formulation**:
- Momentum formulation V_eff = W - β²/(4K) ≥ μ²/2 for all parameters
- β²/(4K) correction is tiny compared to the large W from Skyrme coupling
- **No bound states** for any (e, gc, μ, λ) tested, grid-converged

**Physical conclusion**: The B=1 hedgehog soliton background does not support trapped degenerate modes under either coupling mechanism tested. The Skyrme coupling is always repulsive; the covariant derivative coupling has zero cross-term due to Lie algebra orthogonality. This null result is reinterpreted in Phase 12: bound states are the wrong framework for short-lived particles.

---

## Phase 12 — Bosons as Transient Geometric Deformations ✓ STEPS 1–4 COMPLETE

### Conceptual reframing

The Phase 10–11 null result (no bound states) assumed bosons should be **geometrically stable** — persistent trapped modes in a soliton background. But this contradicts a fundamental distinction in CHPT:

- **Geometrically stable particles** (proton, neutron, electron): topologically protected solitons. Long-lived because the topology cannot unwind. Correspond to bound state eigenvalues.
- **Geometrically unstable particles** (W, Z, mesons): transient field configurations with no topological protection. Short-lived because the geometry is not a local energy minimum. Do NOT correspond to bound states.

Searching for bound states of bosons was looking for the wrong mathematical object. The absence of bound states is **consistent** with bosons being short-lived — it confirms they are not geometrically stable.

### Alternate view: bosons as expanded knots

In a stable soliton, the field is in a compact topological "knot" — energy concentrated in steep field gradients at the core. The knot's compactness minimizes total energy (virial balance E₂ = E₄).

During soliton-soliton interactions, portions of this compact geometry can temporarily **deform or expand**:
- The field is no longer in its energy-minimizing configuration
- The expanded region carries energy proportional to the deformation
- The configuration is geometrically unstable — it relaxes back or radiates away
- The "mass" detected for a boson is the **energy cost of the geometric deformation**, not the mass of a separate persistent object

This explains why bosons "appear" to have high mass in detector readings: the field has briefly left its compact equilibrium, and the energy of the expanded configuration registers as invariant mass in the decay products.

### Key distinction from QFT

In QFT, bosons and fermions have the same ontological status — both are excitations of quantum fields. In CHPT, they have fundamentally different geometries:
- **Fermions**: components of topologically stable soliton configurations
- **Bosons**: transient geometric features of soliton interactions — not objects in themselves, but descriptions of how energy and quantum numbers flow between solitons during deformation

A boson is to a soliton what a wave is to the ocean — a pattern of energy transfer, not a separate thing.

### Evidence from existing simulations

The B+B̄ annihilation simulation (Phase 9) already contains boson-like phenomena:
- At closest approach: E_pot drops from 217 to 74, E_kin surges to 153
- The field in the interaction region is a high-energy, non-equilibrium, non-topological configuration
- 48% of rest mass permanently radiated — energy carried away by field waves
- This transient expanded state is exactly the geometry described above

### Approach: back-calculate boson geometry from experiment

Rather than derive boson properties from the Lagrangian (forward problem), invert the problem:

1. **Fit soliton parameters** (ρ₀, e) to known stable particle properties (proton mass, charge radius) using PDG data
2. **Spectrally analyze scattering simulations** — Fourier transform the field in the interaction region during B+B and B+B̄ collisions. Peaks in the spectrum correspond to characteristic deformation energies.
3. **Compare energy ratios** — check whether spectral peaks fall at the right ratios relative to soliton mass:
   - W/proton ≈ 86 (W mass 80.4 GeV / proton mass 938 MeV)
   - Z/proton ≈ 97 (Z mass 91.2 GeV / proton mass 938 MeV)
4. **Characterize decay** — measure how quickly each spectral component radiates away. The decay rate gives the boson width (lifetime).
5. **Reconstruct geometry** — from the spectral analysis, identify which specific field deformations carry the boson quantum numbers.

### Required experimental data (from PDG)

For parameter fitting:
- Proton mass: 938.272 MeV/c²
- Proton charge radius: 0.8414 fm
- Proton magnetic moment: 2.793 μ_N

Target boson properties:
- W± mass: 80,377 MeV/c², width: 2,085 MeV (lifetime ~3×10⁻²⁵ s)
- Z⁰ mass: 91,188 MeV/c², width: 2,495 MeV (lifetime ~3×10⁻²⁵ s)
- Pion mass: 139.6 MeV/c², width: very narrow (lifetime 2.6×10⁻⁸ s — long for a meson)

### Implementation plan

1. **Parameter fitting**: Map (ρ₀, e) → (proton mass, charge radius) using radial solver + maxwell outputs ✓ DONE
2. **Spectral analysis of existing data**: FFT of post-collision energy time series ✓ DONE
3. **Normal mode analysis**: Breathing mode eigenfrequencies of the hedgehog ✓ DONE
4. ~~High-resolution scattering: Re-run at finer time output~~ (NOT NEEDED — see below)
5. **Comparison with PDG**: Check mass ratios, lifetime ratios, production thresholds

### Step 1 results (parameter fitting)

Conversion factors at e=1, ρ₀=1 (scattering simulation parameters):
- 1 code energy = 9.098 MeV, 1 code length = 0.5624 fm, 1 code time = 1.875×10⁻²⁴ s
- Soliton mass = 103.13 code = 938.3 MeV, radius = 1.495 code = 0.841 fm
- Physical constraint: ρ₀⁴/e² = 0.0259, standard Skyrme: e_std=3.69, F_π=95.0 MeV (3% of exp.)
- Accessible particles at v=0.5c: pion (15.3 code), rho (85.2), Delta (135.4)
- W/Z/Higgs inaccessible (need γ≈86–133)

See `results/paramfit_analysis.md` for full details including PDG data.

### Step 2 results (spectral analysis)

FFT of B+B̄ scattering (v=0.5c, e=1, σ-model, post-collision t>2.0):
- Dominant peak at ω = 2.09 code = 19.1 MeV (breathing mode)
- Second peak at ω = 4.19 code = 38 MeV (2nd harmonic)
- No pion-mass feature (expected: σ-model has massless pions)
- Nyquist limit 572 MeV (rho at 775 MeV unresolvable)

### Step 3 results (normal modes) — DEFINITIVE

**Code**: `src/normal_modes.c` — breathing mode eigenvalue solver (Sturm-Liouville, RK4 shooting)

Computes second variation of E = E₂ + E₄ around hedgehog, solves -(Pg')'+Wg = λmg.

**Key finding: NO breathing resonance in the sigma model.**

The 19 MeV spectral peak (box mode n=6, ω=2.09) matches the soliton light-crossing time 2R/c = 3.0 code time units. However, this mode is only 2% above the free-particle box eigenvalue — it is NOT a soliton-specific resonance. The soliton's effect on the K=0 continuum is a weak perturbation:
- n=0: 86% above free (strong soliton interaction)
- n=3: within 0.3% of free (negligible soliton effect)
- n=6 (19 MeV peak): 2% above free (standing wave, not resonance)

Effective potential: repulsive core (W/m~36 at r=0), shallow well (W/m~-8 at r≈0.5), returns to W/m→0 at large r. Well too shallow for bound states in the massless continuum.

**Moment of inertia** Λ = 141.6 (e=1) → Delta-N splitting = 0.1 MeV. At e=4: Λ=2.32, ΔE=23.5 MeV. Both far from experiment (293.7 MeV).

**Implications**:
1. Meson candidates require angular modes (K=1 for rho, K=2 for f₂), not K=0 breathing
2. Pion mass essential for genuine resonances (creates mass gap)
3. The 19 MeV peak is the soliton's natural ringing frequency, not a meson

See `results/normal_mode_analysis.md` for full analysis.

### Step 4 results (massive pion profile) — DEFINITIVE

**Code**: `src/radial.c` — extended with `-mpi` parameter for pion mass term V=m_pi^2*rho_0^2*int(1-cosf)d^3x.

Modified ODE adds +m_pi^2*r^2*sin(f) to numerator. Virial: E_2 - E_4 + 3*E_V = 0. Mass: Mc^2 = 2E_4 - 2E_V. Asymptotic: f ~ exp(-m_pi*r)/r.

Physical pion mass in code units: m_pi = 0.398 code^{-1} (1/length, NOT energy).

**Key findings at e=1, rho_0=1, m_pi=0.398 (physical pion mass)**:
- Soliton mass: 110.19 code (vs 103.14 massless) — +6.8% increase
- E_2 = 44.73, E_4 = 60.27, E_V = 5.18 (E_2 no longer equals E_4)
- R_rms = 1.264 code = 0.711 fm (vs 1.496 = 0.841 fm) — 15.5% more compact
- Virial: E_2 - E_4 + 3*E_V = 1.4e-8 (machine precision)

**Breathing modes with self-consistent massive profile**:
- **NO bound states below continuum threshold** m_pi^2 = 0.158 (fine scan confirms)
- **NO instabilities** (all eigenvalues positive — the sigma-model artifact is gone)
- Indicial exponent nu = 4.31 (vs 3.79 massless)
- Box modes essentially unchanged (continuum modes >> m_pi^2)

**Previous instability RESOLVED**: Using the sigma-model profile with pion mass gave lambda_0 < 0 (false instability). Using the correct self-consistent massive profile gives all lambda > 0. Always use profiles solved at the same m_pi.

**Conclusion**: K=0 breathing modes are definitively ruled out as meson candidates in both massless and massive cases. Meson identification requires K >= 1 angular modes.

See `results/massive_pion_analysis.md` for full analysis.

---

## Phase 13 — Pion Mass in 3D Scattering ✓ COMPLETE

**Code**: `src/field.c` (pion mass energy + force), `src/scatter.c` (added `-mpi` CLI flag).

**3D pion mass term**: $V_\pi = m_\pi^2 \rho_0^2 (1 - q \cdot s/|q|)$. Force is tangent to $S^3$ ($F \cdot q = 0$, no radial projection needed). Gradient verified to ~10 digits.

**CFL discovery**: The Skyrme $E_4$ term imposes a CFL condition scaling as $h^2$, not $h$. At $N=192$, $L=8$ ($h=0.0833$), $dt$ must be $\leq 0.015$ (not 0.025). At $N=128$, $L=8$ ($h=0.125$), $dt=0.025$ is stable. Finer grids require *smaller* timesteps — the $E_4$ fourth-derivative stencil dominates.

**Massive scattering results** ($v=0.5c$, $e=1$, $m_\pi=0.493$, $N=192$, $L=8$, $dt=0.015$):
- **B+B̄** ($Q=0$): Stable to $t \approx 2.7$. $E_{\text{kin}}: 7 \to 96$, $Q = 0.000$ preserved throughout. Solitons approach from $r=7.7$ to $r=5.1$. Energy blowup at $t=3.0$ (close approach exceeds grid resolution).
- **B+B** ($Q=2$): $Q$ drops from 2.0 to 0.05 at $t \approx 2.1$ (topology loss at close approach, same as massless case). Energy blows up at $t=3.0$.
- Both channels qualitatively similar to massless scattering. Pion mass slightly increases interaction strength but does not change the qualitative picture.

**Profile**: `data/profiles/profile_massive_e1_mpi0.493.dat` (solved at physical pion mass).

---

## Phase 14 — K=1 Angular Modes ✓ COMPLETE: NO BOUND STATES

**Code**: `src/angular_modes.c` — K=1 perturbation eigenvalue solver (Sturm-Liouville, shooting with momentum formulation).

Grand spin $K = L + I$ decomposition of perturbations around the B=1 hedgehog. Three channels at $K=1$:

### Channel A: Isorotational ($L=0$, $I=1$)
- Uniform isorotation $\delta q = \epsilon_a \sigma_a q$ — no spatial dependence.
- **Zero mode at $\omega = 0$**: This is the isorotational Goldstone mode.
- **The physical pion IS this zero mode**, promoted to a collective coordinate via quantization.
- Pion mass arises from collective quantization: $m_\pi = 3/(2\Lambda)$ where $\Lambda$ is the moment of inertia.
- $\Lambda = 79.18$ (massive profile), $\Lambda = 141.55$ (massless). Delta-N splitting: $\Delta E = 3/(2\Lambda) = 0.17$ MeV at $e=1$ (vs 293.7 MeV experiment — wrong scale at $e=1$).
- No vibrational bound states.

### Channel B: Translational ($L=1$, $I=0$)
- Spatial dipole perturbation $g(r) \hat{r}$ — translation-like.
- **No bound states**: Lowest eigenvalue $\lambda = 0.484$, well above continuum threshold $m_\pi^2 = 0.243$.
- Strong centrifugal + Skyrme barrier prevents binding.

### Channel C: Mixed ($L=1$, $I=1$, $K=1$)
- Pion-like channel with both spatial and isospin structure.
- **No bound states**, but **near-binding**: $W/m$ minimum = 0.252, only **4% above** the continuum threshold $m_\pi^2 = 0.243$.
- **Bug found and fixed**: Pion mass term was $(2/3)m_\pi^2 r^2 \sin^2\!f \cos f$ (wrong — decays to 0 at large $r$, killing the mass gap). Corrected to $m_\pi^2 r^2 \cos f$ (gives correct threshold $m_\pi^2$ at large $r$).
- The 4% margin suggests additional attractive coupling could push this channel into binding.

### Physical implications
1. The pion is a quantized zero mode, not a vibrational bound state — consistent with standard Skyrme model results.
2. All K=0 (breathing) and K=1 (angular) bound state searches are definitively negative.
3. The near-binding in Channel C is the closest any mode comes to a trapped state.
4. Meson candidates (rho, etc.) would require $K \geq 2$ or modified coupling.

See `results/angular_mode_analysis.md` (if generated) or run `bin/angular_modes` for full output.

---

## Code Inventory

`proposal/hopfion_search/` — C (gcc + OpenMP) + Python visualization. `make` builds all into `bin/`.

| Directory | Contents |
|-----------|----------|
| `src/` | C source files and headers |
| `bin/` | Compiled executables |
| `data/profiles/` | Radial soliton profiles (sigma-model, rational map, finite-λ) |
| `data/scatter/` | Soliton-soliton scattering time series |
| `data/maxwell/` | Sourced Maxwell field data |
| `data/archive/` | Early exploratory runs (superseded) |
| `scripts/` | Shell batch scripts, Python visualization |
| `results/` | Analysis documents, archived run logs, figures |

See `data/README.md` for full file inventory with generator commands, parameters, and column formats.

| Source | Binary | Purpose |
|--------|--------|---------|
| `src/radial.c` | `bin/radial_solver` | B=1 shooting method solver |
| `src/rational_map.c` | `bin/rational_map_solver` | B=1–4 rational map ansatz solver |
| `src/finite_lambda.c` | `bin/finite_lambda_solver` | Coupled f-shooting + ρ-BVP Newton |
| `src/scatter.c` | `bin/scatter` | Soliton-soliton scattering (+ 8-component degenerate) |
| `src/verify3d.c` | `bin/verify3d` | 3D initialization from 1D profiles |
| `src/maxwell.c` | `bin/maxwell` | Sourced Maxwell equations |
| `src/veff.c` | `bin/veff` | Effective potential for degenerate modes |
| `src/bound.c` | `bin/bound` | Bound state eigenvalue solver (shooting + bisection) |
| `src/normal_modes.c` | `bin/normal_modes` | Hedgehog breathing mode & moment of inertia |
| `src/angular_modes.c` | `bin/angular_modes` | K=1 angular mode eigenvalue solver (3 channels) |
| `src/coupling.c` | (library) | Degenerate sector coupling terms (E_{2,D}, E_{4,C}, E_int) |
| `src/coupling.h` | (header) | Coupling term interface |
| `src/verify_coupling.c` | `bin/verify_coupling` | Coupling gradient verification |
| `src/field.c` | (library) | 3D energy functional & gradient (4th-order) |
| `src/verify.c` | `bin/verify_gradient` | Gradient verification |
| `src/main.c` | `bin/soliton_search` | 3D gradient flow (historical) |
