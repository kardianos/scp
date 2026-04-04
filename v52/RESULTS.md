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

```
volview -snapshot all -outdir ./frames -composite -brightness 2.0 file.sfa
volview -snapshot 0,5,10 -out frame.webp -azimuth 1.5 -distance 3.0 file.sfa
volview -snapshot animation -out anim.webp file.sfa
```

### 6.4 gen_composite Velocity Boost
Added `-velocity vx,vy,vz` for stamping templates with initial velocities.

### 6.5 Seed Generators Created
- `sfa/seed/gen_spherical_seed.c` — Phase soliton (failed)
- `sfa/seed/gen_hedgehog_seed.c` — Hedgehog encoding (failed)
- `sfa/seed/gen_anneal.c` — Localized hot spots (failed)
- `sfa/seed/gen_diffuse.c` — Smooth random annealing (promising)

---

## 7. Open Questions

1. **Will the diffuse periodic condensation produce a stable particle?**
   The t=100 frame shows a shrinking sphere with θ shell and P core.
   If it stabilizes, this is the formation pathway.

2. **Why does E_total decrease with periodic BC?** Either missing energy
   terms in the diagnostic, numerical dissipation, or a real effect.

3. **Can the condensed structure be extracted as a template?** If the
   diffuse run produces a stable particle, extract a frame and use it
   as a seed for collision experiments.

4. **UUD vs UDD comparison still pending.** The theoretical predictions
   (19-69% less hardening for UD) need simulation verification. This
   requires either pre-equilibrated templates or the diffuse annealing
   approach to generate matched proton/neutron seeds.

---

## 8. Source File Index

### Theory
- `v52/THEORY.md`, `v52/core_interaction.mac`, `v52/uud_udd_interaction.mac`
- `v52/predict_collision.mac`, `v52/hedgehog_seed.mac`
- `v52/chiral_quark.mac`, `v52/spherical_seed.mac`, `v52/theta_equilibrium.mac`
- `v52/CoreInteraction.lean`, `v52/UUDvsUDD.lean`, `lean/V52/SphericalSeed.lean`

### Analysis
- `v52/expA_analysis.md` — Collision prediction vs reality
- `v52/seed_postmortem.md` — Why spherical seed failed
- `v52/mismatch_profile.c` — Mismatch extraction tool

### Seed Generators
- `sfa/seed/gen_spherical_seed.c` — Phase soliton
- `sfa/seed/gen_hedgehog_seed.c` — Hedgehog encoding
- `sfa/seed/gen_anneal.c` — Localized hot spots
- `sfa/seed/gen_diffuse.c` — Smooth random annealing

### Simulation Data
- `/space/scp/v52/expA/` — Exp A UU collision (52 GB SFA, 349K diag)
- `/space/scp/v52/anneal_test/` — Annealing experiments
- `/space/scp/v52/spherical_test/` — Spherical seed test
- `/space/scp/v52/hedgehog_test/` — Hedgehog seed test
