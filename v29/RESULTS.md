# V29 Results (In Progress)

## Tier 1 Tests Complete: T1, T2, T4

### Critical Discovery: The Mass Term

The V28 bimodal sweet spot uses m=1.50 (m²=2.25) in the dynamics.
Earlier documentation incorrectly stated "no mass term." V29-T2 clarified
the role of mass by testing all combinations:

| Config | fc | trans_l2 | torsion | |P| | What it tells us |
|--------|------|----------|---------|------|------------------|
| mass²=0, m_init=0 | 0.29 | 0.004 | 1.31 | 1.97 | Braid lives but no quadrupole |
| mass²=0, m_init=1.5 | 0.24 | 0.004 | 7.12 | 1.25 | More init energy → more torsion, still no quad |
| mass²=2.25, m_init=0 | 0.54 | 0.24 | 0.10 | 0.004 | Mass confines → quadrupole, but no torsion |
| **V28 (mass²=2.25, m_init=1.5)** | **0.93** | **0.21** | **1.02** | **0.73** | **Both — the bimodal synergy** |

The mass term plays two competing roles:
- **Confines** the field (fc↑, trans_l2↑) — good for gravity proxy
- **Suppresses propagation** (torsion↓, |P|↓) — bad for EM proxy

The V28 sweet spot is where both coexist. The initialization velocity
(m_init=1.5 → ω=√(k²+m²)) kick-starts torsion circulation that persists
despite the confining mass term.

### T1: Thermal Bath (mass²=0) — NEGATIVE

With mass²=0 and periodic BC, the braid dissolves by T=50 regardless
of noise level (fc drops to 0.088). The initialization radiation fills
the box.

### T1b: Thermal Bath (mass²=2.25 + pre-equilibration) — POSITIVE

**The braid survives in a periodic box and approaches thermal equilibrium.**

Pre-equilibration (T=300 with absorbing BC) sheds initialization energy,
then periodic BC + noise for T=700 more.

| Noise | avg fc | dE/dt | Status |
|-------|--------|-------|--------|
| 0.000 | 0.823 | -0.01 | Stable, slow radiation |
| 0.001 | 0.823 | -0.01 | Identical to control |
| 0.005 | 0.820 | -0.01 | Barely affected |
| 0.010 | 0.812 | -0.01 | Slightly degraded |
| 0.050 | 0.575 | +0.02 | Noise pumps energy in |
| 0.100 | 0.228 | +0.12 | Overwhelmed |

**Thermal equilibrium point**: dE/dt crosses zero at A_noise ≈ 0.02
(2.5% of field amplitude). At this noise level, P_in = P_out.
The braid is in dynamic equilibrium with its radiation bath.

This confirms the user's hypothesis: the braid is NOT "bleeding to death
in vacuum." With a thermal background, it reaches a steady state where
fc ≈ 0.75-0.82, |P| ≈ 0.40, E ≈ 1500-2000.

### T3: Radiation Rate (mass²=2.25, large box L=60, N=192) — STEADY

The braid radiates at a roughly constant rate in open space:

| Time | E_core | fc | dE_core/dt |
|------|--------|------|------------|
| 500 | 3537 | 0.80 | ~-0.5 |
| 1000 | 3315 | 0.83 | ~-0.5 |
| 1500 | 3114 | 0.89 | ~-0.5 |
| 2000 | 2783 | 0.93 | ~-0.5 |

Key: dE/dt ≈ **-0.5 constant** (not power-law decay, not exponential).
Vacuum lifetime ≈ 7500 time units (~190 oscillation periods).

The braid CONCENTRATES as it radiates: fc rises from 0.80→0.93.
Outer layers radiate away, leaving a tighter core ("evaporative cooling").

Combined with T1b: radiation rate of ~0.5 is balanced by thermal
background at A_noise ≈ 0.02. In a thermal universe, the braid is
**indefinitely stable** via dynamic equilibrium.

### T2: m_init Independence — SEE TABLE ABOVE

The braid structure (fc, |P|, torsion) is robust to m_init at mass²=0.
The bimodal synergy (trans_l2 AND torsion) requires mass²=2.25 AND m_init≈1.5.
"Emergent mass" claim is dead; m=1.50 is a Lagrangian parameter.

### T4: Topological Fragility (mass²=0) — ROBUST STRUCTURE, CHAOTIC WINDING

At mass²=0, braid survives perturbations up to ε=5.0 (6.25× field amplitude).
fc varies only 5%. Winding oscillates chaotically.

### T4b: Topological Fragility (mass²=2.25) — MUCH BETTER LOCALIZATION

| ε | fc (mass²=0) | fc (mass²=2.25) |
|---|-------------|-----------------|
| 0.05 | 0.224 | **0.961** |
| 0.50 | 0.222 | **0.903** |
| 1.00 | 0.219 | **0.947** |
| 3.00 | 0.213 | **0.770** |
| 5.00 | 0.215 | **0.580** |

The mass term transforms the braid: fc=0.96 (vs 0.17 after settling),
winding conserved during settling (-1.000 vs -0.000).

Structural robustness: fc>0.58 even at ε=5.0 (6.25× field amplitude).
No blowup at any amplitude. The braid structure is very robust.

Winding still oscillates dynamically during evolution (at all ε), but this
appears to be a diagnostic sensitivity issue — the phase measured along a
single z-line through the center is noisy even when the global structure
is intact. A volume-integrated topological charge would be more reliable.

## Phase B: Field Study (T10 series) — COMPLETE

The comprehensive field study reveals **structural gaps** in the model:

| Test | Result | Significance |
|------|--------|-------------|
| T10B effective metric | Mixed: attractive center, repulsive shell, birefringent | Not gravity-like |
| T10A wave scattering | NEGATIVE: no measurable scattering | Not EM-like |
| T10C dispersion | No tachyonic instability, min M²=0.82 | No mass emergence |
| T10D Kibble-Zurek | NEGATIVE: no spontaneous braids | Braids are not natural defects |
| T10E impulse response | 3.7% channel mixing at center | Mode conversion exists |
| T10F massless mediator | Mass essential; m=0 destroys braid | Fundamental tension |
| T10G complex upgrade | **POSITIVE**: braid survives, carries charge | Path to EM |
| T10H charged braid force | Charges irrelevant; Yukawa dominates | Need gauge field |

**Bottom line**: The triple-product model produces stable solitons but cannot produce
gravity (no massless mediator) or EM (no gauge symmetry). See T10_field/DISCOVERY_PHASE_B.md
for full analysis and recommended paths forward.

## Remaining Tests

| Test | Status | Priority |
|------|--------|----------|
| T1b thermal (mass²=2.25) | **COMPLETE — POSITIVE** | Done |
| T4b fragility (mass²=2.25) | **COMPLETE — ROBUST** | Done |
| T3 radiation rate | **COMPLETE — STEADY** | Done |
| T5 metric coupling | **COMPLETE — 63% quad strain** | Done |
| T6 universality | **COMPLETE — synergy in band** | Done |
| T7 two-braid | **COMPLETE — REPULSION** | Done |
| T8 boost | Not started | Lower |
| T9 substrate | Not started | Long-term |
| T10 field study | **COMPLETE — SEE ABOVE** | Done |

### T5: Self-Consistent Metric (mass²=2.25) — 63% QUADRUPOLAR STRAIN

The strain tensor h_ij at R=10 is **63% quadrupolar** (deviatoric). This is
a geometric property of the bimodal braid, present even without backreaction.
Metric coupling (α_g up to 0.05) has negligible effect on dynamics.
The braid is stable under all tested coupling strengths.

### T6: Universality — MODERATE BAND

Synergy > 1.0 at 3/10 (μ,κ) pairs: (-30,30), (-41.3,50), (-41.3,20).
This spans μ ∈ [-41, -30], κ ∈ [20, 50] — a band, not a point.
Best synergy: 1.71 at (μ=-30, κ=30), better than V28 sweet spot (1.31).
The transverse quadrupole is the fragile channel — torsion is robust.

### T7: Two-Braid Interaction (mass²=2.25) — REPULSION

Two bimodal braids at D=30 REPEL: ΔD = +3.2 (same twist), +2.5 (opposite).
Same-twist repels MORE (charge-like). This is OPPOSITE to V27's m=0
attraction. The mass term creates a short-range Yukawa tail (range ~0.67)
that doesn't reach D=30. The repulsion may be radiation pressure.
Need closer separations (D=10-15) to test for short-range attraction.
