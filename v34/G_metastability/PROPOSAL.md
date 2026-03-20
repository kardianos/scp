# Track G: Metastability Window — Finding the Sweet Spot for m²

**Goal**: Find the smallest m² that gives practically stable vacuum +
braids, measure the force exponent n at that m², and determine whether
the depletion extends the effective gravitational range beyond the bare
Yukawa length 1/m.

**Background**: V33-C4 showed m²=0 explodes, m²=2.25 is stable, and
there's a smooth transition between them. The -3.4% energy drift at
m²=1.0 is either numerical (fixable with smaller dt) or physical
(genuine instability). This test resolves that ambiguity and maps the
full metastability landscape.

---

## Phase 1: Fine m² Scan (Single Braid Stability)

**Base code**: v33.c (standard equation, single alloc, periodic BC)

For each m² in {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25}:
- N=128, L=20, T=500, single braid, A_bg=0.1
- Record: E_total at each diagnostic step
- Measure: E_drift rate (linear fit to E(t)/E(0))
- Measure: braid survival (does the 5× threshold centroid remain coherent?)
- Save field snapshots at t=0, 100, 200, 300, 400, 500

**Output**: data/G1_scan/m{value}/timeseries.tsv + snapshots

**Success criterion**: identify the LOWEST m² where:
- E_drift < 1% over T=500 (accounting for numerical contribution)
- Braid centroid position is coherent (not jumping randomly)

## Phase 2: dt Convergence at Critical m²

For the critical m² identified in Phase 1 (expected ~1.0-1.5):
- Run at 3 different dt: 0.5×default, 1.0×default, 2.0×default
- N=128, L=20, T=200

If E_drift ∝ dt² → purely numerical → this m² is STABLE
If E_drift is constant across dt → physical instability → raise m²

**Output**: data/G2_dt/dt{factor}/timeseries.tsv

## Phase 3: Force Law at Each Stable m²

For each m² that passes Phase 1 (drift < 1%):
- Two braids, D=15, N=128, L=35, T=200
- Measure ΔD (separation change) → force estimate
- Repeat at D=10, 15, 20, 25, 30 for the best m²
- Fit: F ∝ 1/D^n → extract n

**Key question**: Does n approach 2.0 as m² decreases?

**Output**: data/G3_force/m{value}/timeseries.tsv

## Phase 4: Effective Mass Profile (m_eff map)

For the best stable m², using a settled single-braid snapshot:
- At each grid point, compute:
  P_local = φ₀φ₁φ₂
  V''(P) = μ(1 - 3κP²)/(1 + κP²)³
  m_eff² = m² + V''(P) × (coupling factor)
- Map m_eff(r) as a function of distance from braid center
- Compute effective Yukawa range: R_eff(r) = 1/m_eff(r)
- Does R_eff extend significantly in the depletion zone?

**Output**: data/G4_meff/meff_profile.tsv

## Phase 5: Rapid Expansion Stabilization (if Phase 1 finds marginal m²)

If a marginally stable m² is found (e.g., stable for T=500 but
slowly drifting), test whether cosmological expansion can stabilize it:

- Initialize at the marginal m² with a hot dense field
- Apply expansion: rescale coordinates by a(t) = (1 + H*t)
  (Hubble-like expansion, H chosen to give 2× expansion over T=1000)
- Does the expansion dilute the instability faster than it grows?
- Does the braid survive the expansion?

This tests the cosmological hypothesis: the universe started at
m≈0 (Big Bang), expanded, and the expansion stabilized the vacuum
by diluting the tachyonic modes.

**Output**: data/G5_expand/timeseries.tsv

---

## Implementation Notes

- All sims use v33.c as base (copy and modify for m² parameter)
- The -m flag in v33.c already accepts mass: `-m 1.0` sets m²=1.0
  (NOTE: v33.c computes m² = argv × argv, so `-m 1.0` gives m²=1.0)
- Phase 1 can run 9 sims in parallel (4 threads each, ~36 threads total)
  or batch them 4 at a time on 16 cores
- Each N=128 T=500 sim takes ~5-10 minutes
- Total Phase 1 wall time: ~15-30 minutes

## Expected Outcomes

**Best case**: m²≈0.75-1.0 is stable, force exponent n≈1.9-2.0,
depletion extends range by 2-3×. Gravity is approximately Newtonian
over relevant scales with Yukawa corrections at very long range.

**Worst case**: Stability requires m²≥2.0, no improvement over current.
Force exponent stays at 1.8. Would need equation modification (φ⁴ term
or background change) to get closer to Newton.

**Wild card**: The expansion test (Phase 5) shows that marginal m²
becomes stable under expansion — supporting the cosmological hypothesis
that the universe's expansion itself maintains vacuum stability.
