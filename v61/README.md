# v61 — Nonlinear/curved gravity backreaction + the EW-vev home (R1)

**Date**: 2026-05-26 (kickoff)
**Parents**: [`../v60/lagrangian/CLOSEOUT.md`](../v60/lagrangian/CLOSEOUT.md),
[`../v60/lagrangian/LAGRANGIAN_v60.md`](../v60/lagrangian/LAGRANGIAN_v60.md)

v60's dynamical-Lagrangian loop produced a verified **linearized** dynamical field
theory (gravity ⊕ matter, EP-exact, stable spectrum, genuine time evolution) and
isolated four residual value-conjectures. v61 takes the two recommended next steps
from the v60 closeout:

1. **Nonlinear / curved-space gravity backreaction** — go beyond the linearized,
   flat-space results: the full (nonlinear) Einstein equation sourced by `ρ_grav`,
   curved solutions (Schwarzschild, interior, cosmological), and genuine
   relativistic observables.
2. **The EW-vev home (R1)** — `v = 784 a²` is the last residual without a dynamical
   origin. Find a dynamical mechanism that produces the full-rank democratic
   `End(L)` condensate (`‖·‖²_F = 784 a²`) from the one scale `a_ℓ`.

**Methodology** (unchanged from the v60 loop): one new aspect per generation, each
verified with at least two of {SymPy, Lean, Maxima, C}. Lean modules in `v61/lean/`
build against the v59 Mathlib (`lake env lean`). See [`LOOP_LOG.md`](LOOP_LOG.md).

The bar: either derive the missing structure/observables, or a sharp, verified
statement of why it cannot be derived (a value-conjecture), as in v60.
