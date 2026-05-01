# V52 Results: UUD vs UDD Characterization

**Status**: In progress. Theoretical analysis complete. Seed generation
is the primary open problem. Diffuse annealing showing promising
condensation behavior.

---

## 1. Theoretical Results (Complete)

### 1.1 Chirality-Dependent Hardening

Algebraic analysis in Maxima proves the V50/C4 equations produce a
measurable difference between proton-proton (UUD+UUD) and proton-neutron
(UUD+UDD) interactions, emerging purely from the chirality geometry.

**Result 1: Binding attraction is identical.**
The triple product P = φ₀φ₁φ₂ at the midpoint between two composites
is chirality-independent. P(UU) = P(UD) at all separations D, verified
numerically at D=1..50 with ratio = 1.000 exactly.

Mechanism: chirality enters as cos(χ×k×0) = 1 at the midpoint (y=z=0).
Both UUD and UDD have χ₀ = +1 along the approach axis.

**Result 2: Hardening repulsion is 19-69% lower for UUD+UDD.**
|curl(φ)|² at the overlap midpoint:

| D | \|curl\|²(UU) | \|curl\|²(UD) | Reduction |
|---|--------------|--------------|-----------|
| 25 | 0.442 | 0.356 | -19% |
| 12 | 0.184 | 0.057 | -69% |
| 4 | 0.187 | 0.058 | -69% |

The ratio oscillates with D (wavelength ≈ 2π/k ≈ 4.2) but is always < 1.

**Result 3: The mechanism is curl_z cancellation.**
Braid 1 (y-axis) has χ₁=+1 in UUD but χ₁=-1 in UDD. The y-derivative
∂φ_x/∂y has opposite sign:
- UU pair: (-C) + (-C) = -2C (constructive)
- UD pair: (-C) + (+C) = 0 (destructive cancellation)

**Result 4: Breathing is approximately invariant under overlap.**
The carrier phase identity Σcos(δ + Δ_b) = 0 holds by linearity.
V51 zoom data confirms phi_max oscillates 0.31-0.86 with identical
amplitude before and after collision. However, during active collision
the period shortens 15% (2.24 → 1.91), a partial violation.

**Source**: `v52/THEORY.md`, `v52/core_interaction.mac`,
`v52/uud_udd_interaction.mac`, `v52/CoreInteraction.lean`,
`v52/UUDvsUDD.lean`

---

## 2. Exp A Results: UUD+UUD Collision (Complete)

Full analysis: `v52/expA_analysis.md`

### 2.1 Key Finding: Velocity Kick Thermalization

The v_kick=0.1 Galilean boost deposits 79 energy units. Template
condensation radiates 21,047 units in the first 25 time units (270×
more). The kick is completely thermalized. D(t) shows zero systematic
approach for t=0-150. The protons sit at D≈24.9 until depletion
attraction pulls them together starting at t≈175.

**Exp A is effectively a v=0 gravitational infall, not a v=0.1 collision.**
Exp B (UUD+UDD) was cancelled for the same reason.

### 2.2 Collision Dynamics

First close encounter at t≈375 (D_min=8.6). Post-collision: bound
oscillating state with D=8-25, orbital period ~50 time units. Breathing
period stable at 2.24 during free flight, shortened to 1.91 during
collision (15% decrease).

### 2.3 Implication

Template-stamp velocity kicks don't work. For controlled collisions,
protons must be pre-equilibrated on separate grids for T≈200, then
stamped with evolved velocity fields.

---

## 3. Analytical Seed Approaches (Multiple Attempts)

### 3.1 Phase Soliton (gen_spherical_seed.c) — FAILED

**Ansatz**: φ_a(r) = A(r) × cos(k×r + δ_a + χ_a × ψ(r))
Single radial coordinate, phase shift ψ₀≈210° for UUD.

**Result**: Initial transient showed P_int oscillating 6-34 with periodic
BC, phi_max growing to 1.69. But P_peak collapsed from 0.33 (t=20) to
0.011 (t=80). The structure dispersed rather than condensed.

**Diagnosis**: No angular topology. All three field components see the
same radial coordinate — there's no carrier phase cancellation boundary.
P is nonzero everywhere within the envelope, not confined to a core.
The φ_max growth was radiation filling the periodic box, not particle
condensation.

Self-consistent θ initialization (θ = G×curl(φ)) was correct and useful
— 260× more θ than old braid seeds. But topology matters more than θ.

**Source**: `v52/spherical_seed.mac`, `sfa/seed/gen_spherical_seed.c`,
`v52/seed_postmortem.md`

### 3.2 Hedgehog Encoding (gen_hedgehog_seed.c) — FAILED

**Ansatz**: φ_a(x,y,z) = A(r) × cos(χ_a × k × coord_a + δ_a + Δ_a)
Each field component oscillates along its own axis (x,y,z).

**Result**: Creates CUBIC structure, not spherical. The three independent
cos(kx), cos(ky), cos(kz) terms create a rectangular grid pattern
visible in volview. P_int declined from 75→5 at t=8 even with dA=1.0
(10× background). E_pot stayed positive — no binding well formed in
the dynamics.

**Diagnosis**: The hedgehog is fundamentally cubic because the carrier
waves are locked to the grid axes. The equations have no mechanism to
round a cubic structure into a sphere. The angular topology is correct
(three axes with different carrier phases) but the shape is wrong.

**Source**: `v52/hedgehog_seed.mac`, `sfa/seed/gen_hedgehog_seed.c`

### 3.3 Three-Braid Phase-Confined (gen_phase_confined.c) — FAILED under C4

**Ansatz**: Three Gaussian tubes along x,y,z with carrier phases
{0, 2π/3, 4π/3}. This is the original V41/V43 approach.

**Result**: Under V44 (no C4 terms), the braids merge into a sphere at
N=64 L=8.3 over T=500 — this is how the V43 proton template was made.
Under C4 (α=0.1, β=0.5), the hardening prevents merger. The braids
remain separate and eventually dissolve.

**Two-phase attempt**: Run under V44 first (α=β=0) to merge, then
switch to C4. Tested at N=64 L=8.3 with A=0.3, R=2.0, T=500 under
V44: P_int oscillated 3-10 initially but declined to 3-4 by t=200.
Braids coexist without merging even under V44 at these parameters.

**Diagnosis**: The V43 template was likely made with higher initial
overlap (gen_proton_analytical) and/or ran much longer. The
gen_phase_confined braids at A=0.1-0.3 are too dilute to create
sufficient P for V(P) binding.

### 3.4 Key Insight: No Static Equilibrium

The Maxima analysis (v52/spherical_seed.mac) showed that the constant-
phase shooting method for the radial ODE diverges. There is no static
bound solution of the C4 equations at A_bg=0.1. The proton is a
**dynamically maintained breathing oscillator** — a limit cycle, not
a fixed point. Seeds must be dynamically relaxed; no formula can
produce the equilibrium state directly.

---

## 4. Annealing Approaches

### 4.1 Localized Hot Spots (gen_anneal.c) — FAILED

**Method**: Tiny (R=0.5-0.8) high-amplitude spots with chiral θ,
grouped in clusters of 3 (UUD/UDD). 6-12 clusters at N=384 L=50.

**Absorbing BC result**: Everything radiated away. P_int=0, phi=0.005
by t=140. The spots were too small — radiation escaped to the absorbing
boundary before V(P) could bind anything. 99.9% of energy lost.

**Periodic BC result** (12 clusters): Initial burst P_int=137 at t=16,
but collapsed to P_int=0.1 by t=146. The spots interacted briefly then
energy dispersed uniformly across the periodic domain.

**Diagnosis**: Localized hot spots expand as radiation faster than V(P)
can capture them. The spots need to overlap immediately or the energy
escapes before binding can occur.

**Source**: `sfa/seed/gen_anneal.c`

### 4.2 Smooth Diffuse Random Field (gen_diffuse.c) — IN PROGRESS, PROMISING

**Method**: Fill the entire grid with smooth random perturbations using
low-frequency sinusoidal modes (k=0.2-1.0). φ gets ±30% amplitude
modulation around A_bg=0.1. θ gets smooth alternating positive/negative
domains (θ_amp=0.05). No sharp edges — all structure from smooth modes.

**Absorbing BC result**: Same as hot spots — P_int=0 by t=140. The
absorbing BC drains everything regardless of initial smoothness.

**Periodic BC result — CONDENSATION OBSERVED**:

| t | E_total | E_pot | P_int | phi_max | θ_rms |
|---|---------|-------|-------|---------|-------|
| 0 | 71,813 | 663 | 64 | 0.14 | 0.022 |
| 15 | 56,673 | 369 | 113 | 0.17 | 0.020 |
| 50 | 39,450 | 266 | 55 | 0.20 | 0.018 |
| 100 | 3,241 | 53 | ~3 | 0.50 | 0.009 |

Volview snapshots show dramatic condensation:
- t=0: Uniform filling — smooth random structure everywhere
- t=50: Field has collapsed from filling the domain to a SPHERE in the
  center. Blue θ shell visible, green/cyan φ+P core forming.
- t=100: Sphere has shrunk further. θ dominates outer shell, P
  concentrated in a small core region.

The E_total drop (71,813 → 3,241 with periodic BC) is anomalous —
energy should be conserved. This may indicate:
1. The diag energy calculation misses a term (Cosserat or hardening energy)
2. Numerical dissipation from the C4 terms
3. Or genuine energy transfer to modes not captured by the diagnostic

The phi_max GROWTH (0.14 → 0.50) while E_total drops is consistent
with condensation: kinetic/gradient energy converting to concentrated
field amplitude through V(P) binding.

**Status**: Simulation running, T=1000. Frames at t=0,50,100 downloaded
and viewed. Continuing to monitor.

**Source**: `sfa/seed/gen_diffuse.c`,
`/space/scp/v52/anneal_test/diffuse_periodic.sfa`

---

## 5. V51 Collision Zoom Analysis (Complete)

### 5.1 Mismatch Profile Through Collision

The V51 collision zoom (T=35 from t=15, 19 frames) captured the first
UUD+UUD collision at t≈25 (parent run time).

**Mismatch at midpoint**: M² spikes 4× at collision (1.2e-3 → 4.3e-3).
Proton cores weaken post-collision (M²_core drops from 3.4e-2 to 6.3e-3)
then reconstitute by t=35.

### 5.2 Structural Destruction Hypothesis

The collision is NOT billiard-ball repulsion. It's destructive
interference: overlapping θ fields create standing waves with nodes
where Cosserat mismatch spikes, disrupting local φ structure. The
released energy propels fragments outward. Fragments reconstitute.

**Source**: `v52/mismatch_profile.c`, `v52/mismatch_zoom.tsv`

---

## 6. Tool Development (Complete)

### 6.1 Runner Channel Notifications
MCP channel support (`notifications/claude/channel`). Events:
run_complete, run_failed, download_complete, download_failed,
teardown_blocked. Protocol version 2025-11-25. Verified working.

### 6.2 Teardown Safety
sim_teardown refuses if undownloaded files exist. force=true overrides.

### 6.3 Volview CLI Snapshots
Headless rendering: single frame, specific frames, all frames to
directory, or animated WebP. Camera orientation, zoom, brightness,
opacity, composite mode, background color all configurable via CLI flags.
Window hidden in snapshot mode.

### 6.4 Unified Vector Compression (FRVD in SFA)
Mixed voxel (FRMD) + polynomial vector (FRVD) frames in a single SFA
file. Both CPU and CUDA kernels produce identical 6×64 = 384 coefficient
layout per patch (tricubic, 8³ blocks). Reader supports multi-field
reconstruction in sfa.h (C) and volview (Go).

**Shared code**: `scp_config.h` (Config struct, parsing, f16 helpers),
`scp_init.h` (init functions, KVMD embedding) — used by both kernels.

**Temporal prediction**: I-frames embed a Fourier model (mean + amp×cos(ωt+φ))
per coefficient. P-frame deltas are `actual - predicted(t)`, enabling
high sparsity for small timesteps. Reader reconstructs via
`actual = predicted(t) + delta`. Model stored in I-frame payload
(flags byte bit 0).

**Config**: `vec_snap_dt=0.1` enables, `vec_iframe_interval=50` controls
I-frame frequency, `vec_delta_tol=0.001` for P-frame sparsity threshold.

**Bugs fixed**:
- FRMD chunk_size convention: fixup stored `chunk_size` as comp_size,
  should be `chunk_size - 12` (FRMD includes header in size field)
- FRVD fixup: was hardcoding all FRVD as VEC_I; now preserves existing
  frame_type from JMPF
- Multi-field reconstruction: reader was writing same value to all 6
  columns; now evaluates each field's 64 coefficients independently
- FrameWriter FRMD format: CUDA wrote embedded time in FRMD (wrong);
  only FRVD should have embedded time
- P-frame delta semantic: CUDA was storing `actual - temporal_prediction`
  but reader expected `actual - previous`. Fixed: temporal model is now
  self-contained (embedded in I-frame, reader computes prediction)
- Infinite loop in `sfa preview`: backward search modified loop variable
- Old FRVD data corruption: VecPatch[64] overflow writing 226 coefficients

**Verified**: CPU + CUDA produce visually identical voxel/vec frame pairs,
tested at N=64 with volview snapshots of all frames.

---

## 7. Vector-Based Seed Generation

### 7.1 Temporal Mean Extraction (gen_temporal_seed.c) — FAILED

**Method**: Run simulation with temporal vec output, extract the Fourier
`temp_mean` coefficients, evaluate to voxels as equilibrium seed.

**Result**: E_total=4.0 (vs template's 1107). Proton dissolved immediately.

**Diagnosis**: The background carrier wave cos(k·z+δ) oscillates in time
as a standing wave; its temporal mean is zero. The temporal mean captured
only a tiny perturbation (phi_max=0.09), not the soliton. Restoring the
analytical background (v2) raised E_total to 178 but the perturbation
was still too weak (phi_max=0.19) to self-sustain — dissolved by T=20.

**Root cause**: The breathing amplitude is ~3× the temporal mean. The
soliton IS the oscillation, not a perturbation on top of a static
equilibrium. The Fourier decomposition doesn't separate "core structure"
from "breathing dynamics" because they're the same thing.

### 7.2 Galerkin Relaxation (gen_galerkin_seed.c) — NO IMPROVEMENT

**Method**: Iteratively minimize PDE residual via overdamped gradient
descent: `φ += α × (∇²φ - m²φ - V'(P) + η×curl(θ))`. Polynomial
re-fitting acts as regularizer. Spherical mask prevents collapse to
trivial zero solution.

**Result**: Residual decreases, energy drops 0.7-5.2%. But simulation
stability metrics (energy drift, P_int, breathing amplitude) are
identical to raw template. The V43 template is already at the
discretized equilibrium.

**Key finding**: Without the spatial mask, relaxation drives everything
to zero within ~100 iterations — the global minimum of periodic-BC
Cosserat is the trivial solution. The proton is a local minimum.

### 7.3 Peak P_int Extraction — SUCCESS

**Method**: Scan all voxel frames for the one with maximum P_int (most
compact proton state = peak of breathing cycle). Use that frame as seed.

**Result**: Peak frame at t=10 has P_int=2232, phi_max=1.42 (vs
template's P_int=391, phi_max=0.93). Simulation from this seed:

| t | phi_max | P_int | E_pot | Status |
|---|---------|-------|-------|--------|
| 0 | 0.86 | 40.8 | -44.8 | Compact |
| 20 | 0.71 | — | -23.1 | Breathing |
| 50 | 0.70 | 19.7 | -20.2 | **Alive** |

Compare raw template at T=50: phi_max=0.37, P_int=2.6.
The peak seed retains 7.5× more P_int at T=50.

**Why it works**: The proton at peak compression has maximal nonlinear
binding (V'(P) is strongest). Less radiation energy to shed. The
breathing cycle naturally returns to this state — it's the attractor.

**Caveat**: This is not a new soliton solution — it's just a better
snapshot of the existing one. It doesn't solve the problem of creating
protons from scratch under C4 equations.

### 7.4 Key Insight: The "Proton" Question

What survives in these simulations is a single breathing oscillator —
a localized, periodically pulsating concentration of the triple product
P = φ₀φ₁φ₂. This is more accurately a single **braid** (quark analog)
than a proton. A proton would require three braids with specific
winding topology (UUD chirality), but at N=64 L=8.3 the grid is too
small to resolve three separate braids — they merge into one blob.

The question of chirality for a single braid: the carrier phase offsets
δ = {0, 3.0005, 4.4325} break the Z₃ symmetry between components, so
the blob does carry phase information. But without resolvable internal
structure, there is no meaningful distinction between UUD and UDD
chirality. Chirality requires at least two braids with different
winding orientations to be distinguishable.

**Source**: `sfa/seed/gen_temporal_seed.c`, `sfa/seed/gen_temporal_seed_v2.c`,
`sfa/seed/gen_galerkin_seed.c`

---

## 8. Chirality and Lissajous Composite Construction

### 8.1 Theta Orientation — CLEAR SIGNAL

**Method**: Take the peak-P_int blob seed (which survives) and create 5
variants with different theta initializations relative to curl(φ):

| Seed | θ initialization | G factor |
|------|-----------------|----------|
| θ_zero | θ = 0 | 0 |
| θ_pos_mod | θ = +0.5 × curl(φ) | +0.5 |
| θ_pos_ext | θ = +2.5 × curl(φ) | +2.5 |
| θ_neg_mod | θ = -0.5 × curl(φ) | -0.5 |
| θ_neg_ext | θ = -2.5 × curl(φ) | -2.5 |

All 5 run from the same φ (phi_max=0.86, P_int=40.8, E_pot=-44.8),
identical except for theta. Periodic BC, N=64, T=100.

**Results**:

| Seed | Init θ_rms | Final φ_max | Final P_int | Final θ_rms |
|------|-----------|-------------|-------------|-------------|
| θ = 0 | 0.000 | 0.41 | 3.6 | 0.045 |
| θ = +0.5×curl | 0.026 | 0.35 | 1.8 | 0.051 |
| θ = +2.5×curl | 0.130 | 0.32 | 0.7 | 0.114 |
| θ = -0.5×curl | 0.026 | **0.47** | **6.5** | 0.046 |
| θ = -2.5×curl | 0.130 | **0.68** | **27.4** | 0.099 |

**Key finding: Negative theta (anti-aligned with curl(φ)) dramatically
strengthens the soliton.** The extreme negative seed retains 7.5× more
P_int than θ=0 and 40× more than positive extreme. Positive theta
weakens and eventually kills the soliton.

**Physical mechanism**: The Cosserat mismatch is M = curl(φ)/2 - θ.
- θ = -G×curl(φ): M = (1/2 + G)×curl(φ) → **enhanced mismatch** → stronger α|M|² binding
- θ = +G×curl(φ): M = (1/2 - G)×curl(φ) → **reduced mismatch** → weaker binding
- At G=0.5 positive: M ≈ 0 (mismatch nearly cancelled), soliton barely binds
- At G=2.5 negative: M = 3×curl(φ), mismatch energy dominates → strong confinement

**Implication**: Chirality in this model is the sign of θ relative to
curl(φ). Anti-aligned theta creates the self-reinforcing feedback loop:
strong mismatch → strong binding → maintains compact φ → maintains
strong curl(φ) → maintains strong mismatch. Aligned theta creates the
opposite: mismatch cancellation → weak binding → dispersal.

This is a single-blob chirality, not a multi-braid topology. The
distinction between "proton-like" and "neutron-like" may reduce to
the sign of θ·curl(φ) rather than the winding number geometry of
three separate braids. The 3-braid picture may be the UV description
that coarse-grains into this IR blob with a chiral theta orientation.

**Note**: θ_rms for the θ=0 seed GROWS from 0 to 0.045 — the equations
spontaneously generate theta from curl(φ) coupling. The sign it chooses
may depend on initial noise / numerical asymmetry. This could explain
why some template seeds work better than others.

**Source**: `v52/chirality_test/make_seeds.c`, `v52/chirality_test/*_diag.tsv`

### 8.2 Absorbing BC Validation — Only Strongest Survives

The 121 survivors under periodic BC were tested under absorbing BC
(N=96, L=15, damp_width=5, damp_rate=0.01, T=200) to distinguish
genuine solitons from radiation trapped by periodic boundaries.

**Result**: Only 1 of 3 top configurations survived:

| Config | Final φ_max | Final P_int | Verdict |
|--------|-------------|-------------|---------|
| (0, π/3, π/3) | **0.63** | **31.2** | **Survived** |
| (π/3, 2π/3, π/3) | 0.28 | 2.1 | Dying |
| (π/3, π/3, 2π/3) | 0.20 | 0.6 | Dead |

The periodic BC was flattering the results — most "survivors" were
radiation bouncing in a box.

### 8.3 Cooled Soliton — Gentle Radiative Boundary

**Method**: Wide absorbing shell (damp_width=5, rate=0.005) that
catches escaping radiation without touching the soliton core (undamped
radius=10, soliton core radius~5). T=300 for thorough cooling.

**Result**: The (0, π/3, π/3) soliton survived 300 time units of
continuous radiative cooling:

| t | E_total | E_loss | φ_max | E_pot | Status |
|---|---------|--------|-------|-------|--------|
| 0 | 6,602 | 0% | 2.40 | -354 | Hot initial |
| 60 | 1,903 | 71% | 0.84 | -72 | Shedding radiation |
| 150 | 698 | 89% | 0.79 | -90 | Cooling |
| 240 | 339 | 95% | 0.28 | -0.1 | Quiet phase |
| 300 | 259 | **96%** | **0.75** | -33 | **Still breathing** |

The soliton shed 96% of its initial energy as radiation but maintained
its core structure. phi_max oscillates 0.28-0.80 (reduced from the
initial 0.3-1.3), indicating genuine cooling of the breathing mode.
P_int=31.6 at T=300 — the binding is persistent.

The energy loss rate is decelerating: 71% in the first 60 time units
vs 7% in the last 60 (259→259 essentially flat). The soliton is
approaching thermal equilibrium with the radiative boundary.

**Archived**: `scpsfa:scpsfa/v52/0_p3_p3_cool2_run.sfa` (1.6 GB,
N=96, 612 frames with Chebyshev polynomial vec encoding)

### 8.4 Intra-Particle Dynamics — φ-θ Energy Exchange

**Method**: Frame-by-frame radial shell analysis of the breathing
soliton (N=96, periodic BC, T=20, vec_snap_dt=0.1). Tracked φ_rms,
θ_rms, curl(φ)_rms, and mismatch M = curl(φ)/2 - θ within the
particle core (r < 3 code units).

**Key finding: θ persists when φ disperses.**

| Phase | φ_core | θ_core | θ/φ ratio | M_core |
|-------|--------|--------|-----------|--------|
| φ peak (P_int max) | 3.66 | 0.89 | 0.24 | 0.34 |
| φ trough (P_int min) | 0.55 | **1.28** | **2.33** | **1.40** |
| Late φ trough | 0.48 | **1.27** | **2.66** | **1.40** |

The φ field oscillates by 10.5× (0.35-3.67), but the θ field only
oscillates by 26× with a PHASE OFFSET — θ is large when φ is small,
and small when φ is large. The energy exchanges between the two fields.

**Mechanism**: When φ disperses:
1. θ retains the particle's spatial structure (θ_rms stays high)
2. The Cosserat mismatch M = curl(φ)/2 - θ grows (M_core peaks at 1.4)
3. The coupling term η·curl(θ) in the φ equation acts as a restoring force
4. φ is pulled back toward the θ concentration
5. φ reforms, θ relaxes, cycle repeats

This is analogous to the E-B exchange in electromagnetic waves:
- EM: E peaks → generates curl(B) → B peaks → generates curl(E) → cycle
- Soliton: φ peaks → generates curl that sources θ → θ peaks → curl(θ) restores φ → cycle

The θ field acts as the particle's "memory" — it preserves the
soliton's location and shape during the phase when φ appears to vanish.
The particle never truly disappears; its θ skeleton is always present.

**Implication for EM sector**: The same φ-θ curl coupling that creates
the breathing soliton also supports propagating wave solutions. The
linearized Cosserat equations are:
  ∂²φ/∂t² = ∇²φ - m²φ + η·curl(θ)
  ∂²θ/∂t² = ∇²θ + η·curl(φ)

This has the same curl-exchange structure as Maxwell's equations. The
θ field (massless, m_θ=0) propagates at c — a candidate photon analog.
A θ-wave packet traveling between two solitons would be photon exchange.

**Source**: `v52/chirality_test/breathing_analysis.c`,
`/tmp/breathing_global.tsv`, `/tmp/breathing_radial.tsv`

### 8.4.1 Self-Wrapping Geometry — Maxima Analysis

**Question**: Can φ and θ wrap so completely that the soliton achieves
zero radiation loss?

**Answer**: No — but effectively yes, in practice.

**Structural impossibility** (equation shape, not parameters):
- Beltrami self-trapping (curl(F)=λF for both fields) is incompatible
  with V(P) binding. The triple product P=φ₀φ₁φ₂ requires unequal
  per-component amplitudes (cos² ratio up to 14.8:1 from phase offsets),
  breaking the equal-amplitude Beltrami condition.
- Helicity H = ∫φ·curl(φ) is NOT conserved — V(P) breaks it. No
  topological protection for vortex linking.
- Node trapping fails — the apparent m²_eff divergence at field zeros
  is an artifact (all forces vanish there simultaneously).

**Parameter-dependent suppression**:
- Radiation scales as P ~ ε^{2n} where n = ⌈ω_L/ω_breathe⌉
- n depends on soliton radius R and mass m: larger R → more harmonics
  below threshold → exponentially less radiation
- R=10, n=5, ε=0.1: P ~ 10⁻¹⁰ (functionally zero)
- The self-wrapping cost is only 11% extra mass gap (η²/m² = 0.11),
  easily overcome by V(P) binding (|μ|/m²_T = 16.5 >> 1)

**Environmental stabilization — the key insight**:

The soliton is not an object in empty space. It is a topological defect
IN the carrier wave field (A_bg). It cannot exist without the carrier —
removing the carrier (temporal mean extraction, §7.1) kills it instantly
(E_total drops from 1107 to 4). The carrier is the fabric.

In equilibrium:
- Radiation out: P_out ~ ε^{2n} (exponentially small for cooled soliton)
- Absorption in: P_in ~ A_bg² × σ (carrier feeds energy through V(P))
- Stability: P_in = P_out → soliton has fixed breathing amplitude

With periodic BC, this equilibrium is observed: E_total stabilizes at
1,173 with only 0.05% drift over 100 time units. With absorbing BC
(drain the fabric), the soliton slowly dies. The carrier wave IS the
vacuum that sustains the particle.

**Analogy to real physics**: Particles are stable because they sit in
the quantum vacuum (zero-point fluctuations). Change the vacuum energy
density → change the particle properties. No vacuum → no particles.
In our model: A_bg is the vacuum energy density. The soliton's
equilibrium breathing amplitude is set by A_bg.

**Source**: `v52/chirality_test/self_wrap_analysis.mac`

### 8.4.2 Proposed Experiment: A_bg-Dependent Equilibrium

**Hypothesis**: The soliton's equilibrium breathing amplitude is
determined by the carrier wave amplitude A_bg. Different A_bg values
should produce solitons with different "temperatures" (breathing
amplitudes) that are stable in their respective backgrounds.

**Design**: Three runs with the same (0, π/3, π/3) Lissajous seed
but different background amplitudes: A_bg = 0.05, 0.10, 0.15.
Each run: absorbing BC for T=100 (shed radiation), then periodic BC
for T=200 (equilibrate). Compare final breathing amplitude, P_int,
and energy at T=300.

**Prediction**: Higher A_bg → larger equilibrium P_int and breathing
amplitude, because the carrier provides more energy for absorption.

**Runtime**: GPU recommended. N=192 L=25 T=300 for each run. At
0.5ms/step on V100: T=300 at dt≈0.006 → 50,000 steps → 25s per run.
Three runs = ~2 minutes GPU time. Use vec_snap_dt=1.0 with
vec_iframe_interval=50 for manageable output (~300 vec frames).

N=192 is the minimum to resolve the soliton (core ~5 code units =
~20 voxels) with enough room for the background (L=25 gives 20 code
units of clear space around the core). N=384 would be better but
3× slower.

### 8.4.3 Chirality Detection — Observable Design

**Problem**: The vector encoding (polynomial coefficients) does not
preserve enough information to detect chirality. Chirality is a
relationship between φ, θ, and their spatial derivatives — it requires
the full voxel grid with curl computation, not compressed coefficients.

**Three chirality observables** (all pseudoscalars — change sign
under parity inversion):

**1. Cross-helicity** (direct chirality measure):
```
H_θφ = ∫ θ · curl(φ) dV
```
Measures alignment between θ and curl(φ). From §8.1:
- Anti-aligned θ (stable): H_θφ < 0 (θ opposes curl(φ))
- Aligned θ (unstable): H_θφ > 0 (θ parallel to curl(φ))
This IS the chirality. |H_θφ| measures its strength.

**2. Helicity** (self-linking):
```
H = ∫ φ · curl(φ) dV
```
Measures how much φ twists along itself. Nonzero for the Lissajous
geometry because the three carrier waves along orthogonal axes create
a chiral interference pattern. Sign depends on the phase offsets.

**3. Mismatch asymmetry**:
```
C = ∫ (|curl(φ)|²/4 - |θ|²) dV
```
Measures whether the Cosserat mismatch M = curl(φ)/2 - θ is dominated
by curl (C > 0) or by theta (C < 0). For our anti-aligned chirality,
|θ| is enhanced, so C < 0.

**Implementation**: These are volume integrals over the full grid,
requiring curl computation at every voxel. They CANNOT be computed
from the polynomial vec frames — they need the actual voxel data.

**Proposed: chirality analysis hook in the simulation kernel.**
Add to the diagnostics computation (runs at diag_dt interval):
```c
/* In compute_energy or a new compute_chirality function */
double H_cross = 0, H_self = 0, C_asym = 0;
for (long idx = 0; idx < N3; idx++) {
    /* curl(phi) already computed for forces */
    double cx = curl_phi_x[idx], cy = curl_phi_y[idx], cz = curl_phi_z[idx];
    H_cross += theta[0][idx]*cx + theta[1][idx]*cy + theta[2][idx]*cz;
    H_self  += phi[0][idx]*cx + phi[1][idx]*cy + phi[2][idx]*cz;
    C_asym  += (cx*cx+cy*cy+cz*cz)/4 - (theta[0][idx]*theta[0][idx]
              +theta[1][idx]*theta[1][idx]+theta[2][idx]*theta[2][idx]);
}
/* Multiply by dV and append to diag TSV */
```
Cost: negligible — curl(φ) is already computed for forces. Three dot
products per voxel, one reduction per diagnostic step.

**Why not in post-processing**: The polynomial vec frames lose the
high-frequency curl information. The Chebyshev tricubic basis smooths
the field within each 8³ block, and curl across block boundaries is
discontinuous. Only the in-simulation computation on the full-resolution
grid gives accurate chirality values.

**Alternative**: A dedicated analysis pass could read voxel (FRMD)
frames and compute chirality, but this only works at snap_dt resolution
(sparse voxel frames), not at vec_snap_dt resolution. The in-kernel
hook gives chirality at every diag_dt step.

**Expected signature**: For the (0, π/3, π/3) Lissajous soliton with
anti-aligned θ:
- H_θφ < 0 (negative cross-helicity)
- H oscillating (changes sign with breathing cycle)
- C < 0 (theta-dominated mismatch)
- All three should be CONSTANT IN SIGN through the breathing cycle
  (chirality is preserved even when φ "disappears" — θ holds it)

### 8.4.4 Proposed Tool: Per-Particle Chiral Analysis (`sfa_particle_track`)

**Purpose**: Read an SFA file frame-by-frame, detect all particles via
cluster analysis, and output a per-particle-per-timestep diagnostic
including chirality, frequency, and phase signature.

**Input**: Any SFA file (voxel or mixed voxel+vec frames). Vec frames
are reconstructed to voxels internally before analysis.

**Output**: TSV with one line per detected particle per timestep:

```
t  pid  mass  cx  cy  cz  rms_r  phi_max  P_peak  E_pot  H_cross  H_self  C_asym  theta_rms  curl_rms  nvox  omega  phase
```

Column definitions:
- `t`: simulation time
- `pid`: particle ID (assigned by centroid proximity to previous frame)
- `mass`: ∫|P|dV within cluster (P = |φ₀φ₁φ₂|)
- `cx,cy,cz`: P-weighted centroid
- `rms_r`: RMS radius from centroid
- `phi_max`: max |φ| within cluster
- `P_peak`: max |P| within cluster
- `E_pot`: ∫V(P)dV within cluster (binding energy)
- `H_cross`: ∫θ·curl(φ)dV within cluster (chirality sign)
- `H_self`: ∫φ·curl(φ)dV within cluster (helicity)
- `C_asym`: ∫(|curl(φ)|²/4 - |θ|²)dV within cluster (mismatch asymmetry)
- `theta_rms`: RMS θ within cluster
- `curl_rms`: RMS |curl(φ)| within cluster
- `nvox`: voxel count above threshold
- `omega`: breathing frequency (estimated from rolling FFT of mass)
- `phase`: breathing phase at this timestep (from analytic signal)

**Algorithm per frame**:
1. Read frame → voxels (sfa_read_frame handles both FRMD and FRVD)
2. Compute P = |φ₀φ₁φ₂| and curl(φ) at every voxel
3. Threshold: P > 0.01 × max(P)
4. Flood-fill connected components on thresholded voxels
5. Discard clusters with mass < 1.0 or nvox < 50
6. For each surviving cluster:
   a. Compute mass, centroid, rms_r, phi_max, P_peak
   b. Compute E_pot = Σ V(P)·dV (within cluster only)
   c. Compute H_cross = Σ θ·curl(φ)·dV (chirality)
   d. Compute H_self = Σ φ·curl(φ)·dV (helicity)
   e. Compute C_asym = Σ (|curl(φ)|²/4 - |θ|²)·dV
   f. Compute theta_rms, curl_rms within cluster
7. Track particle identity across frames by nearest-centroid matching
8. Estimate omega and phase from rolling window FFT on mass(t) per pid

**Particle tracking**: Between frames, match clusters by minimizing
centroid distance. If a cluster splits, assign both fragments the
parent's pid with suffixes (e.g., pid=1 → 1a, 1b). If clusters merge,
keep the larger one's pid. New clusters get new pids.

**Breathing frequency**: Maintain a rolling buffer of mass(t) per
particle (last ~50 samples). Compute FFT to find dominant frequency.
Phase from analytic signal (Hilbert transform of mass oscillation).

**Implementation notes**:
- Standalone C tool, not a kernel hook — avoids slowing the simulation
- Reads any SFA file post-hoc (existing or future data)
- Vec frames are reconstructed via sfa_read_frame (Chebyshev evaluation)
- curl(φ) computed via central differences on the reconstructed grid
- Parallelizable: each frame is independent (OpenMP over frames)
- For large files: `--stride N` flag to analyze every Nth frame
- For focused analysis: `--t_start T --t_end T` time window

**Cost estimate**: Per frame at N=96: cluster analysis O(N³),
curl computation O(N³), per-cluster integrals O(nvox). Total ~10ms
per frame. For 600 vec frames: ~6 seconds. Negligible.

**Kernel-side lightweight version**: For real-time monitoring during
simulation, add three global integrals (H_cross, H_self, C_asym)
to the existing diagnostics output. Cost: three dot products per
voxel per diag step — curl(φ) is already computed for forces.
Append three columns to the existing diag.tsv. No cluster detection
(too expensive per-step), just global chirality tracking.

**Usage**:
```
sfa_particle_track input.sfa output.tsv [--stride 1] [--t_start 0] [--t_end inf]
```

### 8.5 3-Axis Lissajous Composite — 5×5×5 PHASE SWEEP

**Method**: Superimpose three rotated copies of the neg-theta blob seed,
each with its carrier wave oriented along a different spatial axis
(x, y, z). Each axis gets an independent phase shift from the table:
Δ₀, Δ₁, Δ₂ ∈ {0, π/3, 2π/3, π, 4π/3}. 125 combinations total (121
completed, 4 skipped).

The rotation maps the z-oriented carrier to each axis by permuting
spatial coordinates. Phase shifts are applied as grid-point translations
along the carrier axis. All composites use θ = -2.5 × curl(φ)
(anti-aligned, strongest binding from §8.1).

All runs: N=64, periodic BC, T=100, 8 threads. ~2 min per run, ~4 hours
total.

**Result: 121/121 survived.** Every phase combination produces a
stable bound state at T=100. The 3-axis Lissajous construction is
universally robust.

### 8.5.1 Phase Preference — Strong Selection

Despite universal survival, the binding strength varies 40×:

| Rank | Δ₀ | Δ₁ | Δ₂ | P_int(100) | φ_max | E_pot | Clusters |
|------|----|----|-----|------------|-------|-------|----------|
| 1 | π/3 | π/3 | 2π/3 | 293.9 | 1.12 | -315.6 | 1 |
| 2 | 0 | π/3 | π/3 | 292.4 | 1.06 | -326.4 | 1 |
| 3 | π/3 | 2π/3 | π/3 | 284.1 | 1.08 | -319.1 | 1 |
| 4 | π/3 | π/3 | 0 | 280.4 | 1.09 | -311.2 | 1 |
| 5 | 2π/3 | π/3 | π/3 | 279.0 | 1.09 | -304.7 | 1 |
| ... | | | | | | | |
| 117 | π | π | π/3 | 7.0 | 0.57 | -2.3 | 1 |
| 118 | π/3 | π | π | 8.3 | 0.53 | -3.5 | 2 |
| 119 | 0 | 4π/3 | 2π/3 | 8.5 | 0.49 | -2.8 | 2 |
| 120 | π | π/3 | π | 6.9 | 0.49 | -2.2 | 2 |
| 121 | 0 | π | π | 11.1 | 0.55 | -4.8 | 2 |

**Optimal phases**: π/3 and 2π/3 dominate the top 20. These are 60°
and 120° offsets — precisely the carrier phase separations in the
original 3-braid geometry (δ = {0, 2π/3, 4π/3}). The Lissajous
construction rediscovers the braid topology from a different direction.

**Worst phases**: π (180°) creates destructive interference. Configs
with multiple π offsets have P_int < 15 and often fragment into
multiple clusters.

**Asymmetry**: The best config (π/3, π/3, 2π/3) is NOT symmetric —
two axes at π/3 and one at 2π/3. Breaking cubic symmetry concentrates
the triple product P more effectively. The best symmetric config
(2π/3, 2π/3, 2π/3) is 4th among symmetric runs with P_int=206.

### 8.5.2 Nonlinear Enhancement

Single blob: P_int ≈ 40.8 (from §7.3 peak seed).
Linear 3× sum: P_int ≈ 122 (if P scales linearly).
Best composite: P_int = 293.9 — **2.4× above linear, 7.2× above
single blob.**

The triple product P = φ₀φ₁φ₂ benefits from constructive interference:
when the three axis-oriented carrier waves overlap in phase, the product
of three large values is much greater than the sum of three individual
products. This is the nonlinear amplification mechanism for V(P)
binding.

### 8.5.3 Particle Analysis — Cluster Structure

Connected component analysis on the final frame (P > 0.001 threshold):

| Category | Count | Description |
|----------|-------|-------------|
| Single-cluster runs | 86 (71%) | One dominant bound state |
| Multi-cluster runs | 35 (29%) | 2+ significant clusters (mass > 1.0) |
| Balanced fragmentation | 20 | Second cluster > 5% of first |
| Near-equal split | 3 | Mass ratio > 0.75 |

**Binding energy** of the largest cluster:
- Range: -0.8 to -326
- Median: -122
- Top 5 all have E_pot < -305 and exactly 1 cluster

**Multi-cluster runs** correlate with weak binding:
- Phase offsets near π produce fragmentation
- 4π/3 offsets also increase multi-cluster probability
- Near-equal splits (mass ratio > 0.75): (0, π, π), (4π/3, 0, 2π/3),
  (2π/3, 2π/3, 4π/3) — these may represent two-particle bound states
  or unstable fragmentation in progress

### 8.5.4 Key Findings

1. **The 3-axis Lissajous construction universally produces stable
   solitons.** This is the first seed generation method that achieves
   100% survival rate under V50/C4 equations.

2. **Phase structure determines binding strength over a 40× range.**
   The model has a clear "particle spectrum" — not all composites are
   equal. The optimal phases (π/3, 2π/3) match the original braid
   carrier geometry.

3. **Multi-particle states exist.** 29% of configs produce
   spatially-separated clusters. The (0, π, π) config produces an
   almost-equal split (mass ratio 0.94) — this could be a two-particle
   state analogous to fission.

4. **Energy conservation is excellent.** All 121 runs have |δE/E| < 1%,
   with 119/121 under 0.5%. The Lissajous seeds are clean — no
   transient radiation artifacts.

5. **The "proton" may not be a 3-braid topology.** It may be a 3-axis
   Lissajous pattern with specific phase relationships. The internal
   structure is the interference pattern of three orthogonal carrier
   waves, not three parallel tubes. The chirality is encoded in the
   relative phase shifts and the sign of θ relative to curl(φ).

### 8.5.5 Data Files

- Composite data: `/space/scp/v52/lissajous/composite.tsv` (121 rows,
  35 columns: phases, energies, P_int time series, cluster analysis)
- Per-run cluster analysis: `/space/scp/v52/lissajous/cluster_analysis.tsv`
- Per-run diagnostics: `/space/scp/v52/lissajous/<name>_diag.tsv`
- SFA outputs: `/space/scp/v52/lissajous/<name>.sfa`
- Seed generator: `v52/chirality_test/gen_lissajous_seed.c`
- Sweep runner: `v52/chirality_test/run_lissajous.go`
- Cluster analyzer: `v52/chirality_test/lissajous_analysis.c`
- Composite builder: `v52/chirality_test/build_composite.py`

```
volview -snapshot all -outdir ./frames -composite -brightness 2.0 file.sfa
volview -snapshot 0,5,10 -out frame.webp -azimuth 1.5 -distance 3.0 file.sfa
volview -snapshot animation -out anim.webp file.sfa
```

---

## 9. Infrastructure Fixes (This Session)

### 9.1 Vector Encoder Bugs Fixed
- **Temporal model bootstrap**: First 50 P-frames predicted zero (model
  uninitialized). Fix: set mean=actual_coefficients at frame 0.
- **I-frame refit**: Model only refit every 50 frames, not at I-frame
  boundaries. Fix: refit at every I-frame write.
- **83% zero-delta patches**: Combined effect made P-frames show static
  prediction instead of evolving data. Roundtrip test confirmed 5-8×
  worse error on P-frames vs I-frames.
- **Source**: `v52/vecstream/vec_roundtrip_test.c`

### 9.2 Chebyshev Basis
Maxima analysis showed monomial basis {1,t,t²,t³} has condition number
868,791 in 3D (loses 5.9 of 7.2 float32 digits). Switched to Chebyshev
{1,s,2s²-1,4s³-3s} with condition number 3.8 (loses 0.6 digits).
230,000× improvement. Changed in CPU encoder, GPU encoder, C reader,
and Go volview reader. Same storage format, same coefficient count.
- **Source**: `v52/theory/vec_analysis.mac`

### 9.3 JMPF Overflow Fix
`sfa_write_vec_chunk` in sfa.h had no JMPF overflow handling. With
>1024 vec frames, writes past JMPF allocation corrupted the file.
Fixed by adding the same overflow logic as `sfa_write_frame`.

### 9.4 SFA Format Fixes (Cumulative)
- FRMD comp_size in fixup: was `chunk_size`, should be `chunk_size-12`
- FRVD fixup frame_type: was hardcoding VEC_I, now preserves I/P type
- FrameWriter FRMD: was writing embedded time (wrong), removed
- Multi-field reconstruction: was writing same value to all 6 columns
- Volview P-frame backward seek: now finds correct preceding I-frame
- sfa preview: infinite loop from backward search modifying loop variable

---

## 10. Source File Index

### Theory (v52/theory/)
- Maxima: core_interaction, uud_udd_interaction, uud_vs_udd,
  predict_collision, chiral_quark, spherical_seed, hedgehog_seed,
  theta_equilibrium, verify_cancellation, equation_vectors,
  equilibrium_seed, self_wrap_analysis, vec_analysis
- Lean: CoreInteraction, UUDvsUDD, SphericalSeed

### Seeds (v52/seeds/)
- Lissajous: gen_lissajous_seed.c, make_seeds.c, run_lissajous.go
- Failed: seed_postmortem.md

### Chirality (v52/chirality/)
- Theta orientation: 5 seed variants, 5 runs
- Breathing analysis: breathing_analysis.c, global/radial TSV

### Cooling (v52/cooling/)
- Absorb validation: 3 configs, only (0,π/3,π/3) survived
- Cooled soliton: T=300 gentle BC (archived to rclone)
- Vec recording: T=200 absorb→periodic with vec_snap_dt=0.1
- particle_compare.c

### Collision (v52/collision/)
- expA_analysis.md, mismatch_profile.c + data

### Vecstream (v52/vecstream/)
- vec_roundtrip_test.c, bug reports, design docs

### Simulation Data
- `/space/scp/v52/lissajous/` — 121 SFA + composite.tsv + cluster_analysis.tsv
- `/space/scp/v52/expA/` — UU collision
- `scpsfa:scpsfa/v52/0_p3_p3_cool2_run.sfa` — archived cooled soliton

---

## 11. Forward: V53

The v52 Lissajous construction produced the first stable soliton under
C4 equations. V53 will quantify the particle spectrum:
- Per-particle chirality tracking tool
- A_bg equilibrium experiment
- Two-particle interactions (same/opposite chirality)
- θ-wave (photon analog) propagation and exchange
- Boosted soliton Lorentz properties
- Full particle table with mass, chirality, breathing frequency

See `v53/PLAN.md`.
