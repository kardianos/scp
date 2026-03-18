# T10: Field Study — The Field Itself, Not the Particle

## Core Question
What does the field DO? How does it propagate, scatter, curve, and
self-interact? Gravity, EM, and mass are field properties — they must
be visible in the field's own behavior, independent of any soliton.

## T10A: Wave Scattering
**Question**: How do two wave packets scatter in the triple-product field?

Two Gaussian wave packets collide head-on. Measure:
- Transmission vs reflection coefficients
- Angular dependence (is there a preferred scattering angle?)
- Polarization mixing (does field 0 scatter into field 1?)
- Does the scattering amplitude depend on relative phase? (→ charge-like)
- Cross-section vs energy (→ force law)

If the triple-product V(P) causes wave-wave scattering that depends on
the relative phase of the three fields, this IS the EM interaction at
the field level. No braid needed.

## T10B: Effective Metric Extraction
**Question**: Does a background field configuration create an effective
metric for test waves?

1. Set up a STATIC background: a frozen braid profile (no time evolution
   of the background, only of test waves)
2. Send small-amplitude test wave packets through the background
3. Measure: speed, bending angle, time delay
4. Extract the effective metric g_ij(x) that reproduces the observed
   wave propagation: ds² = g_ij dx^i dx^j
5. Does this metric have Schwarzschild-like structure near the braid?

This directly tests whether the braid CURVES the field around it.
If test waves slow down and bend toward the braid, that IS gravity
at the field level.

## T10C: Dispersion Relation
**Question**: What is ω(k) for small perturbations?

Linearize the EOM around:
(a) The vacuum φ=0
(b) A uniform condensate φ_a = φ₀
(c) The braid background (spatially varying)

For (a): ω² = k² + m² (trivial, just three Klein-Gordon fields)
For (b): the triple-product coupling modifies the dispersion. If there
are modes with ω² < k², the field has an effective negative mass² →
tachyonic → spontaneous structure formation.
For (c): the braid background creates a band structure (like electrons
in a crystal). Gaps in the band structure = mass generation.

This tells us whether mass can EMERGE from the field configuration
rather than being put in by hand.

## T10D: Phase Transition / Defect Formation
**Question**: Does the field spontaneously form braided structures
when cooled from a high-energy state?

1. Start with uniform high-temperature random field (all modes excited)
2. Slowly cool (reduce kinetic energy via thermostat)
3. Monitor for spontaneous formation of localized structures
4. If braids form spontaneously: they ARE the natural defects of this field
5. If NOT: braids are artifacts of our specific initialization

This is the Kibble-Zurek mechanism applied to our field. It connects
to T9 (substrate) — if the field has a phase transition that produces
braids, the substrate must support that transition.

## T10E: Field Response Function (Green's Function)
**Question**: How does the field respond to a localized impulse?

1. Start with vacuum (φ=0) or condensate
2. Apply a sharp impulse: δφ_a(x₀, t₀) = ε
3. Measure the response φ_a(x, t) at all later times and positions
4. Extract the retarded Green's function G_R(x-x₀, t-t₀)
5. Fourier transform → spectral function

The spectral function tells us EVERYTHING about what propagates in
the field: the speed, the tensor structure, the decay rate. If there's
a massless spin-2 pole, that's a graviton. If there's a massless spin-1
pole, that's a photon.

## Priority Order
1. **T10B** (effective metric) — most direct test of gravity at field level
2. **T10A** (wave scattering) — most direct test of EM at field level
3. **T10C** (dispersion) — mass generation mechanism
4. **T10E** (Green's function) — complete spectral characterization
5. **T10D** (phase transition) — spontaneous defect formation

## IMPORTANT: Dynamics Mass
Use mass²=2.25 (m=1.50) in all tests unless specifically testing mass
dependence. The braid background for T10B uses BIMODAL params as-is.

## Grid & Runtime
- T10A: N=128, L=30, T=200 — ~5 min per scattering config
- T10B: N=128, L=20, T=200 — ~5 min per test wave direction
- T10C: N=64 (linear analysis, fast) — ~1 min per mode
- T10D: N=96, L=20, T=1000 — ~15 min per cooling rate
- T10E: N=128, L=30, T=300 — ~8 min per impulse direction
