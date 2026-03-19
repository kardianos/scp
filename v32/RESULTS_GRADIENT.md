# V32 Gradient Coupling: (∇ρ/ρ)·∇φ

## The Idea

The standard wave equation ∇²φ assumes uniform medium. The correct
equation on a non-uniform medium is the divergence form:

    ∂²φ/∂t² = (1/ρ)∇·(ρ∇φ) - m²φ - V'(φ)
             = ∇²φ + (∇ρ/ρ)·∇φ - m²φ - V'(φ)

The gradient coupling (∇ρ/ρ)·∇φ should make the braid drift along ∇ρ.

## Results: REPULSION + BLOWUP

| α | D(0) | D(300) | ΔD | E(300)/E(0) | Status |
|---|------|--------|-----|-------------|--------|
| 0.0 | 20.0 | 21.7 | +1.7 | 0.9× | Neutral, stable |
| 0.5 | 20.0 | 27.9 | +7.9 | 13× | Repulsion + energy growth |
| 1.0 | 20.0 | 30.3 | +10.3 | 460000× | Blowup |

## Root Cause: Self-Interaction Dominates

The gradient term (∇ρ/ρ)·∇φ is dominated by the braid's OWN gradients:
- ∇ρ at the braid edge: huge (ρ transitions from ~3 to ~0.03)
- ∇φ at the braid edge: also huge (field oscillates in the braid)
- Product: massive self-force that pushes the braid's own field outward
- The other braid's gradient at D=20 is tiny by comparison

This is the same self-interaction problem that plagues all single-field
approaches. The braid's own structure overwhelms any external signal.

## The Fundamental Tension

For gravity to work in a single field:
- The braid must respond to EXTERNAL gradients (from other braids)
- But NOT respond to its OWN gradient (self-interaction)
- In a single field, there's no way to separate "self" from "other"

The M7 split (S/B) solves this by construction: S doesn't feel its own
effect on B, only the other braid's effect. But the split is artificial.

The divergence form equation is MATHEMATICALLY correct for a wave in
a varying medium, but PHYSICALLY the braid IS the medium — it doesn't
propagate through itself. The gradient coupling only makes sense for
the gradient of the AMBIENT field, not the braid's own structure.

## What This Means

The gradient coupling is the RIGHT PHYSICS but applied to the WRONG QUANTITY.
It should couple to the gradient of the LARGE-SCALE density (ambient field),
not the SMALL-SCALE density (braid structure).

This brings us back to scale separation: gravity is a LONG-RANGE effect
that depends on the large-scale field configuration, while the braid is
a SHORT-RANGE structure. Any coupling that doesn't separate these scales
will be dominated by self-interaction.

Possible fixes:
1. Smooth ρ over a large scale before computing ∇ρ (but user rejected smoothing)
2. Use the SPH approach where particle clustering naturally separates scales
3. Frequency filtering (separate fast braid oscillation from slow ambient)
4. The M7 split (works but is artificial)
5. Compute ∇ρ from a time-averaged ρ (the braid oscillation averages out,
   leaving only the slow ambient gradient)
