# 11 — Relativity

This chapter addresses how CHPT reproduces the experimentally verified predictions of special and general relativity. Relativity is among the most precisely tested theories in physics; any failure to reproduce it is fatal.

Depends on: [01_field_axioms.md](01_field_axioms.md), [03_propagation.md](03_propagation.md), [07_gravity.md](07_gravity.md)

---

## Special Relativity

### What Must Be Reproduced

Special relativity (SR) rests on two postulates:
1. The laws of physics are the same in all inertial reference frames.
2. The speed of light c is the same in all inertial frames.

From these follow:
- Lorentz transformations between frames.
- Time dilation: moving clocks run slow by factor gamma = 1/sqrt(1 - v^2/c^2).
- Length contraction: moving objects are shorter by 1/gamma along the direction of motion.
- Mass-energy equivalence: E = mc^2.
- Relativistic momentum: p = gamma * m * v.
- Addition of velocities: v_total = (v1 + v2) / (1 + v1*v2/c^2).

All of these are confirmed experimentally to extraordinary precision (particle accelerators routinely rely on relativistic kinematics; GPS satellites correct for time dilation daily).

### CHPT Approach

CHPT posits a field that fills space with a propagation speed limit c (Axiom 5). The field dynamics are fundamentally **Process-Based**, not Object-Based (see [03_propagation.md](03_propagation.md)).

Because particles (knots) are self-sustaining patterns *of* the field rather than objects *in* the field, they have no "rest frame" relative to a medium. Their internal dynamics scale exactly with their translational motion through the field's process-time structure. Lorentz Invariance is therefore a **mathematical consequence** of the field equation, not an imposed constraint or an emergent approximation.

### Time Dilation from Knot Dynamics

In CHPT, a knot's "clock" is its internal oscillation. The frequency of this oscillation determines the knot's experienced time.

For a knot moving at velocity v through the field:
- Part of the knot's internal dynamics must now "keep up" with the translational motion.
- The internal oscillation frequency decreases because the available field dynamics are shared between translation and internal vibration.
- The factor by which the oscillation slows is exactly gamma — IF the field equation is Lorentz-invariant.

This is not a new prediction; it's a restatement of time dilation in field language. The point is that Lorentz-invariant field equations automatically produce this behavior for any localized solution (knot) moving through the field. This was proven by Einstein and independently by many others.

### Length Contraction

Similarly, a Lorentz-invariant field equation ensures that the spatial extent of a moving knot is contracted by 1/gamma in the direction of motion. This is a mathematical consequence of Lorentz covariance applied to localized field configurations.

### E = mc^2

If the field equation is Lorentz-invariant, the energy of a stationary knot of mass m (integrated density excess) is:

    E_rest = m * c^2

This is a theorem, not an assumption: for any Lorentz-invariant system, rest energy equals mass times c^2. The derivation is standard and applies directly to CHPT knots.

### Relativistic Momentum

For a knot moving at velocity v:

    E^2 = (pc)^2 + (mc^2)^2
    p = gamma * m * v

Again, automatic for Lorentz-invariant systems.

### Assessment

CHPT adopts a Lorentz-invariant field equation (as derived in `spec/math/03_dynamics.md`). All of special relativity is automatically reproduced. There is nothing to derive — it's built into the mathematical framework. This guarantees correctness but means CHPT makes no new predictions in the SR sector (no Lorentz violations).

---

## General Relativity

### What Must Be Reproduced

General relativity (GR) extends SR to include gravity. Key predictions:

1. Gravitational redshift (confirmed).
2. Light bending by mass (confirmed).
3. Gravitational time dilation (confirmed).
4. Perihelion precession of Mercury (confirmed).
5. Gravitational waves (confirmed by LIGO).
6. Frame-dragging (confirmed by Gravity Probe B).
7. Black hole event horizons (strong evidence from EHT, LIGO ringdowns).
8. Expansion of the universe (confirmed by Hubble observations, CMB).

### CHPT vs. GR: Fundamental Difference

GR describes gravity as spacetime curvature. Mass-energy curves spacetime, and curved spacetime tells matter how to move.

CHPT describes gravity as density gradients in a flat spacetime field. Mass (knots) depletes the field, and the resulting density gradient pushes matter toward the depletion zone.

These are conceptually different but can be mathematically equivalent in certain regimes. The question is whether they are equivalent in ALL regimes.

### Weak-Field Regime (Newtonian Limit)

Both GR and CHPT reduce to Newtonian gravity (F ~ Mm/r^2) in the weak-field, slow-motion limit. This is guaranteed as long as CHPT's density gradient produces a 1/r^2 force law, which it does by geometric construction.

### Post-Newtonian Regime

The interesting tests are the post-Newtonian corrections — the first-order deviations from Newtonian gravity. GR predicts specific corrections that have been precisely measured:

- **Perihelion precession**: GR predicts 43"/century for Mercury. Matched to 0.1%.
- **Shapiro delay**: Light passing near the Sun is delayed by ~200 microseconds. Matched to 0.001%.
- **Nordtvedt effect**: Different bodies fall at the same rate in a gravitational field (strong equivalence principle). Tested by lunar laser ranging.

CHPT must produce these exact corrections. In the post-Newtonian expansion, GR has specific coefficients (parameterized by the PPN parameters beta and gamma, both equal to 1 in GR). Any deviation from beta = gamma = 1 is detectable and would rule out CHPT.

### Unknown: PPN Parameters

The Parameterized Post-Newtonian (PPN) framework characterizes any metric theory of gravity by 10 parameters. GR predicts specific values. CHPT must either:

1. Be reformulable as a metric theory (density gradients equivalent to an effective metric) with PPN parameters matching GR.
2. Not be a metric theory, in which case it must still reproduce all observations to the same precision.

**This is a hard constraint.** The PPN parameters are measured to parts-per-million precision in some cases. CHPT's density-gradient gravity must be quantitatively equivalent to GR at this precision, or the theory is ruled out.

### Can Density Gradients Mimic Spacetime Curvature?

There is a well-known formal analogy between wave propagation in a non-uniform medium and wave propagation in curved spacetime. This is the "analogue gravity" program (Unruh, 1981):

- A fluid with non-uniform flow velocity creates an effective metric for sound waves.
- Sound in a flowing fluid behaves exactly like light in curved spacetime.
- A sufficiently extreme flow can create a "sonic horizon" analogous to a black hole event horizon.

CHPT's density field could function as such an analogue: null-rotors propagating in a non-uniform density field experience an effective metric determined by the density distribution. If this effective metric matches the Schwarzschild metric (for a static mass) and the Kerr metric (for a rotating mass), then CHPT reproduces all of GR's predictions in the corresponding regimes.

**This is the most promising path for CHPT's gravitational sector**: show that the density field creates an effective metric for null-rotors, and that this effective metric satisfies Einstein's field equations (or a close equivalent).

### Strong-Field Regime (Black Holes)

The original proposal claims "no event horizons — gradients refract info outward." This is a strong claim that contradicts GR and the growing body of observational evidence for event horizons.

**Evidence for horizons**:
- LIGO/Virgo gravitational wave observations of binary black hole mergers match GR predictions for objects WITH horizons.
- The Event Horizon Telescope (EHT) images of M87* and Sgr A* show features consistent with the photon sphere of a black hole.
- No surface emission has been detected from black hole candidates (expected if they had a surface instead of a horizon).

**CHPT's challenge**: If the density gradient model produces an effective metric, and if that metric has a coordinate singularity (horizon), then CHPT DOES predict horizons — contrary to the original claim. The "no horizon" claim may be based on the intuition that a continuous density field can't have a discontinuity, but horizons in GR are not discontinuities — they are surfaces where the effective metric has a specific property (g_tt = 0). A continuous density field can certainly produce this.

**Recommended approach**: Do not claim "no horizons" as a prediction until the effective metric has been computed. It may turn out that CHPT predicts horizons after all, in which case the theory is consistent with observations. If it genuinely predicts no horizons, that may be a falsification point given current evidence.

---

## The Equivalence Principle

The equivalence principle (EP) has three forms:

1. **Weak EP**: All bodies fall at the same rate in a gravitational field. (Galileo, confirmed to 10^-15 by MICROSCOPE satellite.)
2. **Einstein EP**: In a small region, the effects of gravity are indistinguishable from acceleration. (The basis of GR.)
3. **Strong EP**: The gravitational interaction energy itself gravitates. (Tested by lunar laser ranging.)

### CHPT and the Weak EP

In CHPT, every knot is a density excess in the same field. The density gradient (gravitational field) acts on every knot proportionally to its density excess (mass). The "acceleration" is:

    a = F/m = (gradient pressure * mass) / mass = gradient pressure

The mass cancels. All knots experience the same acceleration in a given gradient. The weak EP holds by construction.

### CHPT and the Einstein EP

The Einstein EP requires that gravity and acceleration are locally indistinguishable. In CHPT, an accelerating knot experiences changes in its internal oscillation (from the field dynamics in the accelerating frame) that are identical to those produced by a density gradient. Whether this holds exactly depends on the field equation. If the field equation is Lorentz-covariant and the density gradient produces an effective metric, the Einstein EP follows from the equivalence of the effective metric to an accelerating-frame metric (Rindler metric).

### CHPT and the Strong EP

This requires that the gravitational binding energy of a composite knot also gravitates. In CHPT, binding energy is part of the total density distribution, so it contributes to the depletion zone and therefore to gravity. The strong EP should hold if all forms of energy (kinetic, potential, binding) are on equal footing in the density field.

---

## Summary

| Relativistic Phenomenon | CHPT Status | Risk Level |
|------------------------|-------------|------------|
| Lorentz invariance | Built-in (Process Ontology) | None |
| Time dilation | Automatic from Lorentz-invariant knots | None |
| Length contraction | Automatic from Lorentz-invariant knots | None |
| E = mc^2 | Automatic | None |
| Newtonian gravity 1/r^2 | Geometric dilution | None |
| Gravitational redshift | Density gradient effect on null-rotors | Needs quantification |
| Light bending | Refractive bending in density gradient | Needs quantification |
| Perihelion precession | Post-Newtonian density corrections | Unknown, high priority |
| Gravitational waves (2 pol.) | Depends on field type (scalar vs. tensor) | Risk: wrong polarization count |
| Frame-dragging | Rotating density asymmetry | Not addressed |
| Event horizons | Effective metric may produce them | Do not pre-commit to "no horizons" |
| Weak equivalence principle | Mass cancellation | Built in |
| Strong equivalence principle | All energy gravitates equally | Plausible, not proven |
