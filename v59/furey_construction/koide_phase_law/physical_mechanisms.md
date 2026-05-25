# Physical mechanisms for  φ = Q/3  (why, not just that)

The law `φ = Q/3` (≡ `9·arg ξ = 1 + 2|ξ|²`) is empirically tight (10⁻⁵) but has no
derived cause.  This document collects *physical* mechanisms that could coerce it —
each stated so it makes a falsifiable prediction beyond "fits 2/9", with an honest
verdict.  Mathematical exclusions (not π-rational, not the Weinberg cosine, not a
naive Z₃ potential) are in `CONCEPT.md §3`; here we ask what dynamics could produce it.

## The physical object

ξ is the **off-diagonal element of a circulant (Z₃-family-symmetric) charged-lepton
mass matrix** — i.e. the order parameter of family-symmetry breaking (a "flavon" /
"yukawaon" VEV in the Koide–Harari–Haut–Weyers tradition).  Two real data:
its **magnitude** |ξ|² = 1/2 (Koide) and its **phase** arg ξ = 2/9 (this law).
A mechanism must explain why a single complex VEV lands on the locus
`9·arg ξ = 1 + 2|ξ|²` and, on it, at the point |ξ|² = 1/2.

Crucially the law is **lepton-specific** (quarks miss), so the mechanism must invoke
something the leptons have and the quarks do not — in v59 that is the **color
singlet**: the lepton sector is the `ℂ`-line `{0,7}` carrying the complex structure
`J_c = γ₀γ₅` (`lean/ColorSU3.lean`), while quarks live in the color triplet `3⊕3̄`.

---

## M1 — Spontaneous CP in a flavon potential  (most conventional)

*Idea.*  arg ξ is the spontaneously-chosen CP phase at the minimum of a Z₃-symmetric
flavon potential V(ξ).  The magnitude is fixed by the radial minimum (|ξ|²=1/2), the
phase by an angular term.

*Test / problem.*  A pure Z₃ angular term ∝ cos 3φ extremises at 3φ ∈ πℤ, i.e. φ ∈
{0, π/3} — *not* 2/9 (proven, `PhaseExclusions.cos6_potential…`).  To land at 2/9 the
potential must mix radial and angular pieces so that the minimum traces the locus
`9φ = 1+2|ξ|²`.  That requires a *specific* coupling — generically tunable, hence not
predictive unless the coupling is itself structural.

*Verdict.*  Possible but not predictive as stated.  Would become predictive only if the
radial–angular coupling were fixed by the same `dimG₂/dimSpin7` data that fixes Q.
**Open lead, low confidence.**

## M2 — Geometric (Berry) phase of the generation cycle  (best "why the /3")

*Idea.*  The Z₃ shift S transports the family index around a closed cycle.  arg ξ is the
*geometric phase per generation* accumulated in this transport, so the total winding is
3φ.  The law 3φ = Q says **the closed-cycle winding equals the Koide ratio** — and the
"/3" is automatic (geometric phase shared over the three generations).

*Why attractive.*  It is the only candidate that *predicts the factor 1/3 = N_gen*
without inserting it.  It also explains lepton-specificity if the relevant Berry
connection lives on the lepton `ℂ`-line `{0,7}` (where `J_c` defines the complex
phase), with the quark triplet carrying a *different* (non-abelian) holonomy.

*Test / problem.*  Needs the connection whose holonomy over the Z₃ cycle is Q.  Q is a
*dimension ratio* (dimG₂/dimSpin7), not obviously a holonomy.  The bridge would be: the
Berry curvature integrated over the generation cycle = (a Chern-Simons / index density)
that evaluates to dimG₂/dimSpin7.  Concretely this is the **J∘Z₃ computation** flagged
as the open T1.1 target: orient `J_c` relative to the generation shift `S` and compute
the transported phase.  Not yet done (needs the generation↔ideal Witt map).

*Verdict.*  **Most promising physically.**  Predicts the 1/3 and the lepton-specificity;
the missing piece is a clean derivation that the cycle holonomy = Q.  Pursue first.

## M3 — ξ on the silent-direction sphere with a phase-fixing constraint

*Idea.*  In v58/v59 the Koide condition is |ξ|²=1/2 — ξ lies on the "silent-direction"
sphere `S³ ⊂ ℍ` of radius 1/√2 (`lean/SilentDirection.lean`).  The law adds a *phase*
constraint, pinning ξ to a single point on that sphere.  If the silent direction has an
*intrinsic* phase reference (a preferred `J`-axis), arg ξ is the angle between ξ and
that reference, and `9φ = 1+2|ξ|²` is the equation of the orbit.

*Test.*  The reference axis is exactly `J_c = γ₀γ₅ ∈ Λ²` (now pinned, `ColorSU3`).  So
arg ξ would be the angle of the lepton ξ relative to `J_c`'s plane.  Computing that
angle and checking it is 2/9 is, again, the J∘Z₃ computation (overlaps M2).

*Verdict.*  Same missing ingredient as M2 (the Witt map / explicit ξ-embedding); M2 and
M3 are two readings of one mechanism.  **Pursue with M2.**

## M4 — Hierarchy-from-phase: the phase keeps the electron light

*Idea.*  φ is not arbitrary because it sets the *lightest* mass.  At φ near π/12 the
electron mass vanishes (`√m_e/a = 1+√2 cos(φ+2π/3) → 0`).  The observed φ = 2/9 sits
just inside this edge (gap ≈ 0.04, `study.py §4`), so the electron is light but nonzero.
A dynamical principle that drives m_e toward (but not to) zero — e.g. an approximate
chiral protection of the first generation — would fix φ near π/12, with the offset
2/9 − π/12 set by the protection-breaking.

*Test / problem.*  This explains *why φ is in the small-m_e regime* but not why exactly
2/9 = Q/3 rather than π/12 (a geometric angle we already excluded!).  Note π/12 *is*
π-rational, 2/9 is not — so the true value is shifted off the chiral-limit angle by a
non-geometric amount, consistent with M2 (a holonomy = Q correction).

*Verdict.*  A *constraint* (φ lies in the light-electron window) more than a derivation;
useful as a consistency check and a hint that the value is "near a special angle, shifted
by Q".  **Supporting evidence for M2, not a standalone mechanism.**

## M5 — RG fixed point / self-consistency

*Idea.*  ξ is generated (radiatively or by a gap equation) rather than a free VEV, and
the self-consistent solution sits at arg ξ = 2/9.

*Verdict.*  Generic; no concrete realization tying the fixed point to Q.  **Lowest
priority** — record only to keep the option list complete.

---

## Synthesis & recommended path

> **UPDATE 2026-05-24 — M2/M3 is RULED OUT** (see `holonomy_result.md`).  The
> inter-generational connection was built (the Brannen ring) and the holonomy test run:
> the flux `arg(ξ³)=3φ` is *not* forced to `Q` by any structural curvature (the v58
> connection is an exact 1-form ⇒ zero holonomy; discrete/center fluxes are π-rational;
> the associator gives integers; the color `Z₃` fixes the lepton singlet).  The clincher
> (machine-checked `PhaseExclusions.koide_not_pi_rational`): a structural holonomy is
> π-rational, but `Q = 2/3 rad` is not π-rational.  So `φ = Q/3` is **not** a geometric
> phase.  The leading candidate is dead.

Original assessment (now superseded for M2/M3): the physically-motivated candidates
converged on **M2/M3 (a geometric phase of the Z₃ generation cycle)** as the only family
predicting the **1/3 = N_gen** and the **lepton-specificity**.  Building it closed it.
M4 explains why the value lies near a chiral-limit angle; M1/M5 are fallbacks.

**Post-exclusion state:** no *predictive* mechanism survives.  The only remaining
structural avenue is a **non-geometric algebraic/spectral relation** fixing `cos 3φ` to the
transcendental `cos(2/3)` (M-new, low confidence — faces the same not-nice-number hurdle);
otherwise `φ = Q/3` is a precise (`~10⁻⁵`) but mechanism-less empirical law.

**Concrete next step (the experiment for this study):** carry out the J∘Z₃ computation —
embed the lepton ξ on the `J_c`-line `{0,7}`, act with the generation shift `S`
(`lean/CyclicShift.lean`), and compute the transported phase over the Z₃ cycle, testing
whether its holonomy is forced to Q = dimG₂/dimSpin7.  Positive ⇒ first derivation of
`φ = Q/3`; negative ⇒ a clean exclusion of the geometric-phase mechanism, leaving the
law as a sharp unexplained regularity.  Either outcome is a result.  This requires the
generation↔minimal-ideal (Witt) map deferred in `7D_Algebra/SevenDAlgebra.lean`.

**The connection M2/M3 needs is supplied by v58** (see `FORWARD_PROPOSAL_v58_dynamics.md`).
The v58 unified field equation produces a *bivector (grade-2 = complex/chiral) connection*
`ω ∝ ∇(M M̃)/|M|²` on a Higgs-like vacuum manifold `M M̃ = v²` — which is exactly the Koide
surface `|ξ|²=1/2`.  So the "missing connection whose Z₃ holonomy is Q" is not a free
postulate: it is the v58 connection, derived from a dynamical law.  The experiment above
becomes **compute the holonomy of the v58 `ω` around the generation Z₃ cycle on the
`|ξ|²=1/2` manifold and test `= Q`** — turning the mechanism search into a concrete
v58+v59 calculation.
