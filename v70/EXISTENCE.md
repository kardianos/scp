# v70 EXISTENCE — Do particles "fully exist," or dissipate and reform?

**Date**: 2026-06-11. **Status**: [measured] — 7-run V100 campaign (m1, m2,
p_dphi0/90/180, b_detuned, b2_detuned; all rc=0) + spectra/phase analysis.
Question (user): do the proton/neutron analogs fully exist continuously, or are
they dynamic objects at a frequency that partly dissipate and reform — and does
that periodicity itself generate force-field effects?

Tools: `v70/analysis/existence_spectrum.py`, `make_w146_profile.py` (on-branch
ω=1.46 profile via the v69 shooter), tracker + lineout phase extraction.
Data: `/space/scp/v70/{m1_ball,m2_breather,p_dphi*,b_detuned,b2_detuned}*`,
new profile `v70/results/gprofile_w146_g005.txt` (Q=114.1, E=173.3, r_half=2.63).

## Answer in brief

The theory contains BOTH kinds of object, and they are sharply distinguished:

1. **Charged ball ("proton-analog"): existence is continuous.** Every real
   component oscillates through zero at the internal clock ω — but the
   phase-invariant density does not. Measured (m1, diag_dt=0.05, t∈[20,100]):
   |Φ|_max flat to 0.40%, s_max flat to 2.4% (a small decaying seed-ringing
   eigenmode at ω≈2.67 — NOT a clock harmonic; no spectral line at ω or 2ω),
   Q_core flat to 0.02%. Continuous existence is not incidental — it is WHY the
   object is stable (v65: blinking density ⇒ above-gap radiation ⇒ death).
2. **Neutral object ("breather", ±ω superposition): existence blinks at 100%
   depth.** m2 (same profile, equal counter-rotating parts): s_max modulation
   depth = 100.00% — the density passes through ~zero and reforms every
   ≈ 2.4 t.u. (line at 2.65 ≈ 2× its effective clock, harmonics at 5.26, 7.91)
   — while slowly dying (mean amplitude −5% per 150 t.u., visible as concentric
   radiation shells in the renders). The "dissipate and reform" particle exists
   in this theory; it is mortal, in line with the v65 theorems and the rem1
   findings. (Suggestive — no more than that — of the free neutron's mortality.)
3. **The oscillation IS a force field** (the user's second intuition,
   confirmed three ways):

## The clock-interference force, measured

(p runs: two ω=1.42 balls, same charge, D=12, η=0, g=0.05, N=144 L=24)

| Δφ (relative clock phase) | outcome | a_rel measured |
|---|---|---|
| 0 (in phase) | **FUSION at t≈57** despite Coulomb repulsion | −1.0e-3 net (contact ≈ −1.6e-3 ≈ 2.7× Coulomb) |
| π/2 (quadrature) | contact OFF — clean Coulomb separation | +5.58e-4 (Coulomb pred +5.96e-4, ratio 0.94) |
| π (anti-phase) | enhanced escape | +1.20e-3 ≈ Coulomb + contact repulsion |

The short-range force is cos(Δφ) clock interference riding on the universal
e^{−μD} tail overlap (μ = √(m²−ω²) ≈ 0.56): same particles, same distance,
same charges — the force's SIGN is set by relative phase alone. This is the
program's strong-force analog, demonstrated with the gauge sector live.

## Detuned clocks: the force beats at Δω

- **b_detuned (negative control, instructive)**: spinning a 1.42-profile at
  1.46 does NOT make a 1.46 ball — it relaxes onto the branch at ω≈1.419
  (the branch is flat in ω vs Q there), so the realized detuning was only
  ~0.002 rad/t.u. (measured Δφ drift 0.00197 — quantitative match). The pair
  stayed effectively co-phase for 300 t.u. and the contact term simply decayed
  as e^{−0.62 D} (≈ μ_t = 0.564 ✓). An apparent "phase-locking" reading was
  considered and REJECTED: branch relaxation fully accounts for it.
- **b2_detuned (genuine, both balls on-branch)**: ω=1.46 (Q=114) + ω=1.42
  (Q=311) at D=11. Relative clock phase winds at **0.0635 rad/t.u.**
  (≈ the local weff difference; 3 full beat cycles over T=300) and the force
  follows: the pair first APPROACHES (co-phase contact beats Coulomb,
  D 10.9→10.3), then surges apart at a_rel = +2.2e-3 (≈4–10× Coulomb) as Δφ
  sweeps through anti-phase at t≈50, then decays to Coulomb-only once D grows
  and the exponential overlap dies. **A force modulated at the difference
  frequency of two particles' internal clocks — the user's hypothesis,
  positively measured.** No phase-locking at this coupling: the drift runs
  free, so the locking range at D≥11 is < 0.06 in Δω.

## What this does and does not say

- In-model: "protons" (charged balls) fully exist at all times; their
  *components* blink at ω but the object does not. Neutral composites blink
  entirely and slowly die. The blink/oscillation structure is not waste heat —
  it is the origin of the short-range (strong-analog) force: coherent
  e^{−μD}·cos(Δφ) attraction/repulsion, beating at Δω for detuned particles.
  The long-range (EM-analog) force comes from the static gauge charge and is
  phase-blind.
- NOT claimed: any quantitative mapping to real nucleons. The neutron-mortality
  parallel is an analogy. (The η-drain — the model's weak-analog dissipation
  channel — was off in all these runs, η=0.)

## Caveats

- b2 windowed a(Δφ) is smeared (window 60 vs beat period 99) and the late
  windows (D>20 in L=24) are wall-image contaminated (v70 FINDINGS §4).
  The clean evidence is the early-time approach→surge sequence + the Δφ(t)
  slope.
- m1's residual 2.4% s-ripple (ω≈2.67) is unexplained in detail (seed ringing:
  amplitude eigenmode and/or trace ±ω contamination from the f32→f16 seed);
  it decays and carries no clock-harmonic line.
- E-drift in b2 was −5.6% (lighter ball relaxation radiating through the
  sponge) vs ≤0.05% in other runs; Q dropped 0.7%. Does not affect the phase
  analysis (phases measured directly per frame).
