# T12 Detailed Characterization: M0, M2, M4, M7

## Setup
N=128, L=30, T=500, A_bg=0.1, mass²=2.25 (BIMODAL params)
50 time snapshots, 100 radial bins, 6 full field dumps per mechanism.
60 output files total across 4 mechanisms.

## THE KEY RESULT: Depletion is Power-Law, Not Yukawa

ALL four mechanisms (including control M0) show the depletion profile
follows a POWER LAW δρ ~ A/r^α, NOT an exponential:

| Mechanism | α (power law) | λ (exponential) | Better fit |
|-----------|--------------|-----------------|-----------|
| M0 (control) | 1.887 | 6.14 | Power law (residual 0.37 vs 0.88) |
| M2 (back-pressure) | 1.886 | 6.15 | Power law (residual 0.37 vs 0.88) |
| M4 (saturating) | 2.109 | 5.36 | Power law (residual 0.47 vs 1.41) |
| M7 (two-component) | 1.730 | 7.38 | Power law (residual 0.47 vs 0.72) |

This is significant: despite the mass m=1.5 creating Yukawa tails for
individual field excitations (e^{-mr}/r), the DEPLETION of the background
follows a power law. The depletion is a collective/nonlinear effect that
bypasses the single-particle mass gap.

The exponent α ≈ 1.7-2.1 is between 1/r (Newtonian gravity) and 1/r²
(dipolar force). For true gravity we need α → 1. M7 is closest (1.73).

## Mechanism Comparison

### M0 (Control) vs M2 (Back-Pressure)
M2 is IDENTICAL to M0 at this resolution. The gradient force (β=0.01)
is too weak to affect the dynamics. A much larger β might help but risks
destabilizing the braid.

**Verdict**: M2 adds nothing. Drop it or increase β dramatically.

### M4 (Saturating Potential)
- Steepest falloff (α=2.11) — furthest from 1/r
- Lowest fc (0.40) — braid degrades more than M0
- Highest |P| (0.38) and more total energy preserved (6257 vs 6152)
- The V_depl = λ(ρ₀-½Σφ²)² restoring force pushes field amplitudes
  toward ρ₀, which interferes with the braid core

**Verdict**: Self-limiting works but damages the braid and steepens the
profile. Not ideal for gravity.

### M7 (Two-Component)
- BEST braid preservation (fc=0.93 vs 0.51 for M0)
- Shallowest power law (α=1.73) — CLOSEST to 1/r
- Highest |P| (0.48)
- More energy in the far field (403 at r>20 vs 325 for M0) — B field
  maintains structure at large r because it's a separate dynamical field

**Verdict**: Clear winner. The S/B separation protects the braid while
allowing the background to develop a nearly 1/r depletion profile.

## Why α ≈ 1.7 and Not 1.0?

The depletion exponent α=1.73 (M7) is between 1/r and 1/r². Possible
explanations:

1. **Finite time**: The 1/r profile may not have fully developed by T=500.
   The depletion wavefront is still propagating outward. At later times,
   the far-field depletion may deepen, moving α toward 1.

2. **Cylindrical geometry**: The braid is a cylinder along z (periodic BC),
   not a point source. A cylindrical source in 3D creates a 1/r_perp
   (perpendicular distance) profile, which when z-averaged gives 1/r_perp.
   But our measurement is radial in the xy-plane → expect 1/r_perp ≈ 1/r.
   The observed α=1.73 > 1 may be from the finite tube radius smearing.

3. **Nonlinear effects**: The triple-product coupling creates a structured
   depletion (not simple radial decay). The profile may be a sum of
   1/r (monopolar) and higher multipoles that steepen the apparent α.

4. **Background dynamics**: The background waves scatter off the braid,
   creating standing wave patterns that modulate the depletion.

## Energy Flow Analysis (M7)

| Shell | E(T=0) | E(T=500) | Change | Interpretation |
|-------|--------|----------|--------|----------------|
| r < 5 (core) | 5437 | 1727 | -3710 | Braid radiates |
| 5-10 | 1141 | 1023 | -118 | Slight depletion |
| 10-15 | 584 | 1429 | +845 | Energy ACCUMULATES |
| 15-20 | 397 | 1108 | +711 | Energy accumulates |
| 20-25 | 271 | 403 | +132 | Slight accumulation |

The braid core loses energy (radiation), which accumulates in the
r=10-20 shell. This is the OPPOSITE of depletion — it's CONCENTRATION
at intermediate radii. The "depletion" measured in the differential
(braid vs control) is actually the braid REDISTRIBUTING background
energy toward intermediate radii, leaving the far field slightly depleted.

This is consistent with Mechanism 5's conservation picture: the braid
redistributes, not consumes. The depletion at large r is because the
braid pulls energy inward from the far field, concentrating it at r≈10-15.

## Conceptual Summary

### What the equations do:

**M4**: V_depl = λ(ρ₀ - ½Σφ²)² adds a PENALTY for field amplitude
deviation from ρ₀. This is Lagrangian (energy-conserving) and acts like
a spring pulling |φ| back to ρ₀. It self-limits depletion but also
limits the braid amplitude, degrading fc.

**M7**: L = L_S(S) + L_B(B) - g×(ΣS²)(ΣB²) separates structure from
background. The coupling converts B→S at the braid core and S→B via
radiation. Self-limiting because depletion of B reduces the conversion
rate. Braid-safe because S dynamics are nearly independent of B.

### Possible combinations:

**M4+M7**: Apply the saturating potential ONLY to B fields:
V_depl = λ(ρ₀_B - ½ΣB²)². This limits B depletion without touching
the braid (S). Best of both: braid protection + depletion regulation.

**M7 with variable g**: Make the coupling strength g depend on the
local B density: g_eff = g × B²/(B² + B₀²). When B is depleted,
g_eff drops → consumption slows → natural feedback.

## Recommendations

1. **Use M7 as the base model** going forward
2. Test at larger domain (L=60-100) and longer time (T=1000-2000) to
   see if α evolves toward 1.0
3. Test M4+M7 combination (saturating potential on B fields only)
4. Measure the RATE of depletion development: when does δρ(r) stabilize?
5. Two M7 braids: does the depletion gradient create attraction?

## Data Files
- charMX_timeseries.tsv: 50-row time series (t, fc, |P|, wind, E, E_shells)
- charMX_profiles.tsv: radial profiles at all snapshots
- charMX_rho_xy_tNNN.tsv: xy energy density slices (6 per mechanism)
- charMX_zaxis_tNNN.tsv: z-axis field profiles (6 per mechanism)
- charMX_fit.tsv: power-law and exponential fit parameters
