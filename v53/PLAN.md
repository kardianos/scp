# V53 Plan: Quantifying Stable Particles

**Goal**: Characterize the particle spectrum of the V50/C4 Cosserat
equations. Quantify stationary particles, composite states, and
intrinsically moving (radiation) modes.

---

## Stage 1: Particle Characterization Tool

Build `sfa_particle_track` (designed in v52 §8.4.4) to enable all
subsequent stages.

### 1.1 Core Tool
- Per-frame cluster detection (flood-fill on P > threshold)
- Per-particle diagnostics: mass, centroid, radius, E_pot, phi_max
- Chirality observables: H_cross, H_self, C_asym per particle
- Particle tracking across frames (nearest-centroid matching)
- Breathing frequency + phase from rolling FFT

### 1.2 Kernel Hook (Lightweight)
- Add H_cross, H_self, C_asym as columns in diag.tsv
- Uses curl(φ) already computed for forces — negligible cost
- Enables real-time chirality monitoring during runs

### 1.3 Validation
- Run on the v52 cooled soliton (0_p3_p3_cool2_run.sfa)
- Verify H_cross < 0 (anti-aligned chirality)
- Verify chirality sign is constant through breathing cycle
- Compare per-particle vs global chirality integrals

---

## Stage 2: Stationary Particles

### 2.1 A_bg Equilibrium Experiment
Test whether soliton breathing amplitude is determined by background
carrier wave density.

**Runs**: Three (0, π/3, π/3) Lissajous seeds at A_bg = 0.05, 0.10, 0.15
- N=192, L=25 (sufficient room for soliton + clean background)
- Absorbing BC for T=100 (shed radiation), periodic for T=200
- GPU recommended: ~25s per run on V100

**Measure**: Equilibrium P_int, breathing amplitude, breathing frequency
as function of A_bg. Prediction: higher A_bg → larger equilibrium P_int.

### 2.2 Chirality Spectrum
From the v52 Lissajous sweep (121 configurations), select the top 10
surviving under absorbing BC and measure their chirality signatures.

**Question**: Do different phase offsets (Δ₀, Δ₁, Δ₂) produce
different H_cross values? If so, the Lissajous phases encode a
discrete chirality spectrum — different "particle species."

### 2.3 Mass Measurement
For each stable particle, compute the effective mass:
- E_rest = E_total(equilibrium) - E_background
- m_eff = E_rest / c² (in code units, c=1)
- Compare with the Skyrme model prediction: M = K × e/ρ₀³ × 9.1 MeV

### 2.4 Cooling Protocol
Establish a standard cooling procedure:
- Absorbing BC (damp_width=5, rate=0.005) for T=300
- Verify E_total has flattened (< 1% change in last 50 time units)
- Switch to periodic for T=100 to measure equilibrium properties
- Record as the particle's "ground state"

---

## Stage 3: Composite Particles

### 3.1 Two-Particle States
Place two cooled solitons on a large grid (N=384, L=50) at separation
D=15-25. Both with same chirality (H_cross < 0).

**Measure**:
- Do they attract or repel?
- At what separation does interaction begin?
- Does a bound two-particle state form?
- What is the binding energy (E_pair - 2×E_single)?

### 3.2 Chirality Interaction
Same as 3.1 but with OPPOSITE chirality solitons:
- Soliton A: θ = -G×curl(φ) (H_cross < 0)
- Soliton B: θ = +G×curl(φ) (H_cross > 0) — or the weakly bound version

**Question**: Does chirality affect the interaction? The v52 theory
predicts UUD+UDD has 19-69% less hardening than UUD+UUD. Does this
hold for Lissajous composites?

### 3.3 Multi-Particle Condensation
Start from the v52 diffuse annealing approach but with Lissajous-seeded
initial conditions: 4-8 solitons placed on a large grid.

**Question**: Do they spontaneously form bound composites? What is the
maximum number of solitons that can bind?

---

## Stage 4: Intrinsically Moving States (Radiation)

### 4.1 θ-Wave Propagation
The θ field is massless — it should propagate at c. Generate a θ-wave
packet and verify:
- Propagation speed = c (within discretization error)
- No dispersion (massless → ω = ck exactly)
- Absorption by a soliton (energy transfer from θ-wave to soliton)
- Emission by a perturbed soliton (soliton relaxation emits θ-wave)

### 4.2 Photon Analog
If §4.1 works, attempt photon exchange:
- Two solitons at D=30
- Perturb soliton A (add energy to its breathing mode)
- Observe: does a θ-wave packet propagate from A to B?
- Does B's breathing amplitude increase after the "photon" arrives?
- Measure propagation speed, dispersion, and energy transfer efficiency

### 4.3 Boosted Soliton
Apply a Lorentz boost to a cooled soliton:
- v = 0.1c, 0.3c, 0.5c
- Does it maintain its structure while moving?
- Does the breathing frequency Lorentz-dilate (ω' = ω/γ)?
- Does the chirality (H_cross) change under boost?

---

## Stage 5: Spectral Analysis

### 5.1 Particle Table
Compile a table of all stable particle states found:

| ID | Δ₀ | Δ₁ | Δ₂ | Mass | H_cross | H_self | ω_breathe | Stable? |
|----|----|----|-----|------|---------|--------|-----------|---------|

### 5.2 Selection Rules
From the interaction experiments (§3), determine:
- Which particle pairs attract/repel
- Which composites are stable
- Whether chirality is conserved in interactions
- Whether there are "forbidden" interactions (selection rules)

### 5.3 Comparison to SM
Map the particle table onto known physics:
- Does the chirality spectrum match quark flavors?
- Do the interaction rules match the strong/weak force?
- Is the photon analog (§4) consistent with QED?

---

## Runtime Estimates

| Stage | Grid | T | Runs | GPU Time | CPU Time |
|-------|------|---|------|----------|----------|
| 1 (tool) | — | — | — | 0 | 1 day dev |
| 2.1 | 192³ | 300 | 3 | 2 min | 30 min |
| 2.2 | 96³ | 200 | 10 | 5 min | 3 hr |
| 3.1 | 384³ | 500 | 3 | 30 min | infeasible |
| 3.2 | 384³ | 500 | 3 | 30 min | infeasible |
| 4.1 | 384³ | 200 | 3 | 10 min | infeasible |
| 4.2 | 384³ | 500 | 1 | 15 min | infeasible |
| 4.3 | 192³ | 200 | 3 | 3 min | 30 min |

**Total GPU**: ~2 hours on V100
**Total CPU**: ~4 hours (stages 1, 2.1, 2.2, 4.3 only)

Stages 3-4 require GPU (N=384 doesn't fit in reasonable CPU time).
