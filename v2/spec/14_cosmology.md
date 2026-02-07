# 14 — Cosmology

This chapter examines CHPT's account of the large-scale universe: the Big Bang, expansion, the cosmic microwave background, dark matter, dark energy, and black holes. Cosmology provides some of the most stringent tests of any fundamental theory because it involves extreme conditions and precision measurements.

Depends on: [02_energy_and_density.md](02_energy_and_density.md), [04_knots_and_particles.md](04_knots_and_particles.md), [07_gravity.md](07_gravity.md), [11_relativity.md](11_relativity.md)

---

## The Big Bang

### Standard Cosmology

The universe began approximately 13.8 billion years ago in an extremely hot, dense state and has been expanding and cooling since. The evidence:

1. Hubble expansion (galaxies receding with velocity proportional to distance).
2. Cosmic Microwave Background (CMB) — relic radiation from ~380,000 years after the Bang.
3. Primordial element abundances (H, He, Li ratios match Big Bang nucleosynthesis predictions).

### CHPT Account

The initial state is a field with extremely high, nearly uniform density — far above the background rho_0. At such extreme density, no stable knots can form (the field is too energetic / too turbulent).

As the field dilutes (expands/spreads), it passes through density thresholds where:
1. The simplest knots can first form (lightest particles).
2. Progressively more complex knots become stable.
3. Composite knots form (nuclei, then atoms).

This sequence matches the standard timeline: quark-gluon plasma -> hadrons -> nuclei -> atoms.

### What "Expansion" Means in CHPT

This is a critical point. In standard cosmology, space itself expands (the metric changes). In CHPT (flat spacetime), expansion must mean something different:

**Option A — Field Dilution**: The total density is conserved but spreads over an increasing volume. The background density rho_0 decreases with time. This is conceptually simple but raises questions:
- What drives the spreading? (Pressure? Inertia from initial conditions?)
- If rho_0 changes, does c change? (Since c may depend on field properties.)
- If rho_0 decreases, the vacuum is different at different epochs — violating the assumption that rho_0 is a fixed constant.

**Option B — Effective Metric Expansion**: The density field creates an effective metric (see [11_relativity.md](11_relativity.md)) that is expanding. Physical distances don't change, but the effective metric that null-rotors experience is expanding. This is mathematically equivalent to standard cosmology in GR but philosophically different.

**Option C — Knot Recession**: Knots formed in the early field carry outward momentum from the initial high-energy state. They are still moving apart due to inertia, not because space is expanding. This is a "Milne universe" type model. It works for the Hubble law but fails to explain the CMB power spectrum and the accelerating expansion.

### Recommended Investigation

Option B is the safest, as it is explicitly designed to match GR's predictions. Option A is the most interesting for CHPT specifically — a decreasing rho_0 would have unique consequences. Option C is likely ruled out by observations.

---

## Cosmic Microwave Background (CMB)

### What Must Be Explained

The CMB is a near-perfect blackbody spectrum at T = 2.725 K, isotropic to 1 part in 100,000, with specific small-amplitude fluctuations that encode the initial density perturbations of the universe.

The fluctuation power spectrum (angular distribution of temperature variations) has been measured to extraordinary precision (Planck satellite). Its shape encodes:

- The geometry of the universe (flat, to 0.4%).
- The matter-radiation ratio.
- The baryon-to-photon ratio.
- The spectral index of initial perturbations.

### CHPT Account

When the field cooled enough for atoms to form (recombination), null-rotors were released from their last interaction with chiral knots. These null-rotors have been freely propagating since, cooling as the field dilutes. The current temperature (2.725 K) reflects the ratio of the field density at recombination to the current density.

The fluctuations in the CMB map the density variations in the field at the time of recombination. These density variations are the seeds of galaxy formation.

### Critical Test

The CMB power spectrum has specific acoustic peaks at angular scales of ~1 degree, ~0.5 degrees, ~0.33 degrees, etc. These correspond to sound waves in the early plasma that froze in at recombination. The peak positions and heights depend on precise cosmological parameters.

CHPT must reproduce this power spectrum. If the effective metric approach (Option B above) is adopted, the calculation is identical to the standard one. If CHPT proposes a different mechanism, the power spectrum becomes a stringent test.

### Unknown

- **Inflation**: The standard cosmological model includes a period of exponential expansion (inflation) in the first ~10^-32 seconds. Inflation solves the horizon problem (why the CMB is so uniform) and the flatness problem (why the universe is so close to spatially flat). Does CHPT need inflation? If so, what drives it? If not, how are the horizon and flatness problems solved?

---

## Dark Matter

### Observational Evidence

- Galaxy rotation curves: Stars at the edges of galaxies orbit too fast for the visible mass.
- Galaxy cluster dynamics: Clusters are more massive than their visible content.
- Gravitational lensing: Mass maps show more matter than is visible.
- CMB power spectrum: The heights of acoustic peaks require ~5x more non-baryonic matter than baryonic.
- Large-scale structure: Galaxy distribution matches simulations with cold dark matter.

### CHPT Options (from [07_gravity.md](07_gravity.md))

**Option A — Modified gravity (extended depletion halos)**:
The density gradient from knot aggregates extends further than Newtonian gravity predicts, mimicking extra mass.

**Problem**: MOND-like modifications work for galaxy rotation curves but fail for:
- Bullet Cluster (two clusters colliding; the mass map from lensing is offset from the visible matter, proving the dark matter is a separate substance).
- CMB acoustic peaks (require a pressureless cold matter component, not modified gravity).
- Large-scale structure formation (simulations require cold dark matter particles, not modified gravity).

**Assessment**: Option A is almost certainly insufficient. The Bullet Cluster alone is widely considered to rule out purely gravitational alternatives to dark matter.

**Option B — Dark knots (non-luminous stable particles)**:
Some stable knots are achiral (no charge) and do not emit/absorb null-rotors (no EM interaction). They have mass (density excess) and therefore gravitate. This IS dark matter — a specific type of CHPT knot.

**Assessment**: This is the conservative and likely correct approach. CHPT predicts dark matter as part of its knot spectrum. The question becomes: which specific knot configuration is the dark matter particle? Its mass, cross-sections, and abundance must match observations (WIMP-like with mass ~10-1000 GeV, or axion-like with much lower mass, etc.).

### Recommendation

Adopt Option B. Do not claim to "solve" dark matter — instead, predict that dark matter is a specific knot type and attempt to determine its properties from the field equation.

---

## Dark Energy

### Observational Evidence

The expansion of the universe is accelerating (discovered 1998 via Type Ia supernovae). This requires a component with negative pressure — "dark energy" — constituting ~68% of the total energy budget. The simplest model is a cosmological constant Lambda.

### CHPT Options

**Option A — Vacuum null-rotor pressure**: If the nonlocal hidden-variable extension (see [12_quantum_phenomena.md](12_quantum_phenomena.md)) is adopted, the field may have residual nonlocal fluctuations even in the "ground state" that produce a net pressure at cosmological scales.

**Problem**: This directly contradicts the definition in [02_energy_and_density.md](02_energy_and_density.md), which defines the uniform vacuum as having zero dynamical energy. If the vacuum has fluctuation-driven pressure, it has nonzero energy, and the vacuum energy magnitude must be explained — recreating the cosmological constant problem (QFT predicts 10^120 times too much). CHPT must either (a) keep a truly zero-energy vacuum and find another mechanism for dark energy, or (b) allow a nonzero vacuum energy and explain its tiny value. **This is an unresolved internal contradiction.**

**Option B — Field dilution dynamics**: The field is not in a static state; it is dynamically evolving (diluting). The rate of dilution is not constant and naturally produces an apparent acceleration at late times.

**Problem**: This requires a specific dynamical mechanism. In standard cosmology, a cosmological constant is constant in time. Observations are consistent with constant dark energy (w = -1.0 +/- 0.1). A dynamically varying dark energy (w != -1) is constrained but not ruled out.

**Option C — No dark energy needed**: CHPT modifies gravity at cosmological scales in a way that produces apparent acceleration without dark energy. This is a bold claim that requires demonstration.

### Assessment

Dark energy is poorly understood in ALL frameworks. CHPT is not uniquely disadvantaged here. The honest statement: CHPT has no specific prediction for dark energy and must wait for the field equation to determine what the theory says about cosmological-scale dynamics.

---

## Black Holes

### Standard Physics

A black hole forms when mass is concentrated enough that the escape velocity exceeds c. The event horizon is the surface from within which nothing can escape.

Key properties:
- Singularity at the center (infinite density in GR — universally considered unphysical).
- Event horizon: one-way causal boundary.
- Hawking radiation: quantum effect producing slow evaporation.
- No-hair theorem: characterized only by mass, charge, angular momentum.

### CHPT Account

The original proposal claims: "No event horizons — gradients refract info outward."

As discussed in [11_relativity.md](11_relativity.md), this claim may be premature. If the density field produces an effective metric, and if the density can become large enough, the effective metric can have a horizon. The question is whether the field dynamics prevent density from reaching the critical value.

**Possible CHPT scenario**:
1. Extreme knot concentration creates very deep density depletion zone.
2. The density gradient becomes so steep that null-rotors are severely refracted.
3. At some threshold, null-rotors are refracted into closed orbits — this IS a photon sphere.
4. If the density concentration is sufficient, null-rotors cannot escape — this IS an event horizon.
5. Inside, the density may reach a maximum (field saturation) rather than a singularity.

**CHPT advantage over GR**: GR predicts a singularity (infinite density). CHPT could avoid singularities if the field has a maximum density (saturation). This would resolve the singularity problem while preserving the horizon. This is genuinely interesting.

**CHPT prediction**: Black holes have horizons (matching observations) but no singularities (replacing the infinite-density point with a maximum-density core). The internal structure is a smooth field configuration at saturation density.

### The Information Paradox

If black holes have horizons and eventually evaporate (Hawking radiation), what happens to information about the matter that fell in? GR + QM produces a paradox: the radiation appears thermal (random), suggesting information is lost, violating quantum mechanics.

CHPT with Hawking radiation: The internal field configuration at saturation density contains the full information about the in-fallen knots. As the black hole evaporates, this information is encoded in the correlations of the emitted null-rotors (Hawking radiation). No information is lost — it is just scrambled and difficult to decode.

This is essentially the "black hole complementarity" or "final-state projection" resolution, expressed in CHPT language. It is not novel, but it is consistent.

### Unknown

- **Does CHPT produce Hawking radiation?** Hawking radiation is a quantum effect (particle creation near the horizon). CHPT must derive an equivalent from the nonlocal field dynamics.
- **What is the maximum density?** Is there a natural saturation scale? This determines the internal structure of black holes.

---

## Baryogenesis (Matter-Antimatter Asymmetry)

The observable universe contains far more matter than antimatter. The Sakharov conditions for generating this asymmetry require:
1. Baryon number violation.
2. C and CP violation.
3. Departure from thermal equilibrium.

CHPT must satisfy all three at some epoch in the early universe.

- Condition 1: If baryon number is topologically conserved, it is NEVER violated, and the matter-antimatter asymmetry must have been built into the initial conditions. This is unsatisfying.
- Condition 2: CP violation requires the field dynamics to treat matter and antimatter differently. See [10_conservation_laws.md](10_conservation_laws.md).
- Condition 3: The expanding/cooling field provides departure from equilibrium naturally.

### Unknown

Baryogenesis is a major unsolved problem in cosmology. CHPT does not obviously help, and topological baryon number conservation may actually make it harder.

---

## Summary

| Cosmological Observation | CHPT Account | Status |
|-------------------------|-------------|--------|
| Big Bang (hot dense origin) | High-density field state | Qualitative match |
| Expansion (Hubble law) | Dilution, effective metric, or inertia | Three options, none developed |
| CMB blackbody spectrum | Relic null-rotors from recombination | Standard picture, reworded |
| CMB fluctuation spectrum | Density perturbations at recombination | Requires effective metric approach |
| Primordial element abundances | Knot formation during field cooling | Standard picture, reworded |
| Dark matter | Dark knots (achiral, non-EM) | No specific prediction |
| Dark energy | Unknown mechanism | No specific prediction |
| Inflation | Unknown | Not addressed |
| Black holes | Horizons yes, singularities no (maybe) | Interesting but unproven |
| Baryogenesis | Difficult if baryon number topological | Potential problem |

### Honest Assessment

CHPT's cosmological sector is largely a re-narration of standard cosmology in CHPT vocabulary. The theory adds no new predictions in this domain until the field equation is specified. The most interesting potential contribution — singularity resolution inside black holes — requires demonstrating that the field has a maximum density.
