# 07 — Gravity

This chapter describes how CHPT produces gravitational effects from field density gradients, assesses the proposal against known gravitational phenomenology, and identifies where it succeeds, where it struggles, and where it is silent.

Depends on: [02_energy_and_density.md](02_energy_and_density.md), [04_knots_and_particles.md](04_knots_and_particles.md)

---

## The Back-Pressure Mechanism

CHPT's gravitational mechanism:

1. A knot (particle) concentrates density locally, depleting the surrounding field (see [02_energy_and_density.md](02_energy_and_density.md)).
2. The surrounding field, at higher density than the depletion zone, exerts inward pressure on the depleted region.
3. Another knot in the vicinity also has its own depletion zone.
4. The overlapping depletion zones create a net inward pressure that pushes the two knots toward each other.

This is functionally identical to saying: regions of lower field density (near mass) have lower field pressure, so the higher-pressure surrounding field pushes objects toward each other. Gravity is a pressure differential.

### Why ~1/r^2

In 3D space, the depletion zone around a spherically symmetric knot spreads over the surface of a sphere. The density deficit per unit area at distance r from the knot scales as:

    delta_rho(r) ~ M / (4 pi r^2)

where M is the total density excess of the knot. The force (pressure gradient times area) therefore scales as 1/r^2. This is just the geometric dilution of a conserved quantity (density deficit) over a sphere — the same reason Coulomb's law and Newton's gravity both go as 1/r^2.

### Why Gravity Is Weakest

In CHPT, gravity is the cumulative back-pressure from the depletion zone — an indirect, geometric effect of the density redistribution. The other forces (EM, strong) arise from direct interactions between knot structures (chiral repulsion/attraction, density overlap). Direct interactions are inherently stronger than indirect pressure effects.

Quantitatively: the depletion zone of a knot extends to infinity but becomes vanishingly small at macroscopic distances. The density perturbation is tiny compared to rho_0. Meanwhile, the chiral interaction at close range involves density perturbations comparable to the knot's own density — orders of magnitude stronger.

### Why Gravity Is Universal

Every knot depletes the surrounding field, regardless of its chirality, spin, charge, or internal structure. The depletion depends only on the total density excess (mass). Therefore, every knot creates a gravitational field, and every knot responds to one. This is the universality of gravity.

### Why Gravity Is Only Attractive

The depletion zone always reduces local density below rho_0. Higher surrounding density always pushes inward. There is no mechanism for a knot to create a density EXCESS in its surroundings (that would violate density conservation — the knot's excess must come from somewhere). Therefore, gravity is always attractive.

**Exception?**: If a configuration could exist that depletes its own interior and enhances its surroundings ("anti-knot" or "density hole"), it would repel other knots gravitationally. Whether such configurations are possible depends on the field equation. In standard physics, negative mass is not observed.

---

## Comparison with General Relativity

General relativity (GR) describes gravity as spacetime curvature caused by mass-energy. CHPT describes gravity as density gradients in a flat-spacetime field. These are deeply different frameworks.

### What GR Predicts That CHPT Must Match

| GR Prediction | Experimental Status | CHPT Account |
|---------------|-------------------|--------------|
| Newton's 1/r^2 (weak field) | Confirmed | Yes (geometric dilution) |
| Gravitational redshift | Confirmed (Pound-Rebka, GPS) | See below |
| Light bending by mass | Confirmed (Eddington, lensing) | See below |
| Gravitational time dilation | Confirmed (clocks at altitude) | See below |
| Gravitational waves | Confirmed (LIGO 2015) | See below |
| Frame-dragging | Confirmed (Gravity Probe B) | Unknown |
| Perihelion precession | Confirmed (Mercury) | Unknown |
| Black holes / event horizons | Strong evidence (EHT, LIGO) | Contested |

### Gravitational Redshift

In GR: photons lose energy climbing out of a gravitational well because spacetime is curved.

In CHPT: A null-rotor propagating from a high-density region (near a knot) to a low-density region (far from a knot) passes through a density gradient. If the propagation speed c depends on local density (even slightly), the null-rotor's frequency changes as it transitions between regions. This produces a redshift.

**Alternatively**: If c is strictly constant everywhere, the redshift must come from the energy cost of the null-rotor escaping the density gradient. The null-rotor loses energy to the gradient, reducing its frequency.

Either way, the qualitative effect is correct. The quantitative match (delta_f/f = g*h/c^2 for height h in field g) must be derived from the field equation.

### Light Bending

A null-rotor passing through a density gradient (near a massive knot) is deflected because the gradient acts as a refractive medium — the "optical density" of the field varies with position. This is gravitational lensing.

In GR, the bending angle for grazing incidence at the Sun is 1.75 arcseconds. CHPT must reproduce this number. The calculation requires knowing how the density gradient around a solar-mass knot aggregate affects null-rotor propagation — which requires the field equation.

### Gravitational Time Dilation

In CHPT, a knot's internal oscillation frequency (which constitutes its "clock") depends on the local field density. In a depletion zone (near a massive object), the reduced density slows the internal dynamics. Clocks run slower near mass.

This is conceptually clean and matches the GR prediction qualitatively. Quantitative verification requires the field equation.

### Gravitational Waves

If a massive knot (or pair of knots) accelerates, its depletion zone changes dynamically. This change propagates outward at speed c as a ripple in the density field — a gravitational wave. The wave carries energy and produces oscillating density gradients that squeeze and stretch test objects.

**Critical question**: GR gravitational waves are tensor perturbations (spin-2). They have two polarization states (+ and x). CHPT density waves are scalar perturbations (if the field is scalar). A scalar gravitational wave would have one polarization state, not two. LIGO observations are consistent with two polarizations. This is a potential falsification point: if the field is scalar, CHPT predicts the wrong polarization count for gravitational waves.

**Resolution**: If the field is valued in a geometric algebra (as suggested by the null-rotor concept), gravitational perturbations could be tensor-valued, supporting two polarizations. This again constrains the field type.

### Frame-Dragging

In GR, a rotating mass drags spacetime around it. This has been measured by Gravity Probe B.

In CHPT: A rotating knot aggregate would create a dynamic density pattern that is not spherically symmetric — it has a rotating component. Whether this produces the correct frame-dragging effects is entirely unknown. It requires solving for the density field around a rotating mass distribution.

**Status**: Not addressed. Open problem.

### Perihelion Precession

Mercury's orbital perihelion precesses by 43 arcseconds per century beyond Newtonian predictions. GR explains this exactly.

CHPT must produce a correction to Newtonian gravity at the post-Newtonian level. Whether the back-pressure mechanism naturally produces the correct 1/c^2 correction terms is unknown. This is a specific quantitative test that the theory must pass.

**Status**: Not addressed. Open problem.

---

## Historical Precedent: Le Sage Gravity

The back-pressure/push-gravity idea has historical precedent in Le Sage's theory of gravitation (1748): gravity arises from the pressure of ultramundane particles pushing objects together. This theory was analyzed and rejected for several reasons:

### Le Sage Problems and CHPT Responses

**Problem 1 — Heating**: In Le Sage gravity, objects absorb gravitational mediators, which should heat them. Calculation shows Earth would be heated to stellar temperatures.

**CHPT Response**: There are no discrete mediators being absorbed. Gravity comes from a static density gradient, not from absorption of incoming particles. The gradient is a self-consistent equilibrium state of the field and does not continuously deposit energy. No heating.

**Problem 2 — Shielding**: In Le Sage gravity, body A should partially "shadow" body B from the gravitational flux coming from the direction of A. This would mean gravity between B and a third body C (on the far side of A) would be reduced. This is not observed — gravity is perfectly linear (no shielding).

**CHPT Response**: The density gradient is a property of the field configuration as a whole, not a stream of particles being absorbed. Two knots each create their own depletion zones, and the depletion zones superpose linearly (for weak fields). There is no shielding. However: this linearity holds only in the weak-field limit. In the strong-field regime (very massive/dense objects), nonlinear effects could produce deviations from superposition. Whether these deviations are consistent with observations is TBD.

**Problem 3 — Speed**: Le Sage mediators must travel much faster than c to avoid gravitational aberration effects (the force should point toward where the source IS, not where it WAS). Laplace showed the mediator speed must exceed 10^7 c.

**CHPT Response**: The density gradient is a static (or quasi-static) property of the field. It is established at the speed of light and then maintained self-consistently. For slowly moving sources, the gradient adjusts quasi-statically without aberration issues. For rapidly moving/accelerating sources, the gradient lags — and this lag IS the source of gravitational waves and post-Newtonian corrections. This is actually the correct physics (GR also predicts gravitational radiation from accelerating masses).

### Assessment

CHPT avoids the worst problems of Le Sage gravity because it uses a continuous field gradient rather than discrete particle absorption. However, the detailed quantitative behavior (post-Newtonian effects, strong-field regime) remains unverified.

---

## Dark Matter

Galaxy rotation curves show that stars at large radii orbit faster than Newtonian gravity from visible matter predicts. The standard explanation is dark matter — invisible mass providing additional gravitational pull.

CHPT offers two alternative explanations:

### Option A — Extended Density Halos

The depletion zone around a galaxy's knot aggregate (stars) extends far beyond the visible matter. At large radii, the density gradient is still significant, providing additional back-pressure that flattens the rotation curve. This would mean standard Newtonian analysis underestimates the gravitational effect because it doesn't account for the field depletion extending to large distances.

**Problem**: This is effectively MOND-like (Modified Newtonian Dynamics). MOND works for galaxy rotation curves but fails badly for galaxy cluster dynamics, the CMB power spectrum, and large-scale structure formation. Any CHPT alternative to dark matter must match ALL of these observations, not just rotation curves.

### Option B — Non-Luminous Stable Knots

Some stable knot configurations may be achiral (uncharged) and non-interacting via null-rotors (dark). They would still have mass (density excess) and create depletion zones (gravitational effects). This IS dark matter — just composed of CHPT knots rather than unknown particles.

**This is effectively the standard dark matter hypothesis**, rephrased in CHPT language. It doesn't explain what dark matter is in terms of the theory's fundamentals unless it predicts which specific knot configurations are dark and what their properties would be.

### Recommended Approach

Don't claim to solve the dark matter problem until the theory can make specific, testable predictions. Be honest: CHPT currently has no more insight into dark matter than standard physics does.

---

## Summary

| Gravitational Phenomenon | CHPT Mechanism | Match to Observation |
|-------------------------|---------------|---------------------|
| Newton's 1/r^2 | Geometric dilution of depletion | Qualitative yes |
| Universality | All knots deplete field | Qualitative yes |
| Attractive only | Depletion always reduces density | Qualitative yes |
| Weakest force | Indirect geometric effect | Qualitative yes |
| Redshift | Frequency change in gradient | Not quantified |
| Lensing | Refractive bending in gradient | Not quantified |
| Time dilation | Slowed internal dynamics in gradient | Not quantified |
| Gravitational waves | Propagating density ripples | Polarization TBD |
| Frame-dragging | Rotating density asymmetry? | Not addressed |
| Perihelion precession | Post-Newtonian density corrections? | Not addressed |
| Black holes / horizons | Extreme density gradient | Contested (see [14_cosmology.md](14_cosmology.md)) |
| Dark matter | Extended halos or dark knots | No unique prediction |
