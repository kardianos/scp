# Multivector build, first result: grade-balance gives `t²=½` (right value, right type — not yet forced)

*2026-05-24.  First computation of the multivector-representation program (`09_…md`): treat the
mass operator as a Cl-multivector and look for a **geometric-product magnitude condition** fixing
the Koide amplitude.  **Result: found one that gives `t²=½` exactly** — grade-balance between the
scalar (grade-0) and non-scalar parts of `M`.  It is the *type-correct* condition the audit
predicted (magnitude from the geometric product, not symmetry), and the **first thing to land the
right value**.  Honest caveat: it is a well-typed *reframing* that hits `t²=½`, **not yet a
forced derivation** — the forcing ("why is the lepton grade-balanced?") is the sharpened residual.*

## The condition

Write the lepton mass operator as a multivector `M/a = I + ζ`, `ζ = ξS + ξ̄S²` (the off-diagonal /
non-scalar "hopping" part).  Decompose by grade:

- **scalar / grade-0 part** = `I` → norm² = 1 (per the identity coefficient);
- **non-scalar part** = `ζ` → norm² = `|ξ|² + |ξ̄|² = 2t²`.

**Grade-equipartition** — the scalar and non-scalar parts carry *equal* norm:

> `1 = 2t²  ⟺  t² = ½  ⟺  Koide Q = 2/3.`

Verified numerically (Frobenius norms of the diagonal vs off-diagonal of `M`): the parts are equal
**only** at `t²=½`, exactly where `Q=2/3`.  Equivalently: the scalar part carries exactly **half**
the total norm of `M/a` (`scalar fraction = 1/(1+2t²) = ½ ⟺ t²=½`).

## Why this matters — it vindicates the audit's diagnosis

Every *symmetry*-based attempt at the amplitude failed with a wrong value: maximal-mixing over the
generation group gave `t²∈{0,1}` (`08_…md`), holonomy/energy-min/G₂-inert all failed. The audit
(`09_…md`) predicted the reason — symmetry is magnitude-blind — and that the fix must be a
**magnitude-carrying multivector (geometric-product) condition**. Grade-balance is exactly that:
it is built from the geometric-product **grade decomposition** (scalar vs non-scalar), not from any
group, and it gives `t²=½` *exactly*. **This is the first condition of the right type to produce
the right number** — direct evidence the multivector reframe is the correct direction.

It also recovers, in the right language, the earlier geometric reading: `Q=2/3 ⟺ √m-vector at 45°
to the democratic axis ⟺` equal power in the S₃-trivial and standard reps.  Grade-balance is the
same 45°/equipartition, now phrased as a **grade** equipartition — the type-correct form.

## Lepton-specific (a feature, with the same caveat as before)

The quark sectors are **not** grade-balanced: `t²_d=3/5` (scalar fraction `5/11≈0.45`),
`t²_u=7/9` (scalar fraction `≈0.39`). Only the lepton sits at the balanced point `½`. So the
lepton is the **grade-balanced critical state**; quarks are off-balance. (Consistent with the
lepton being the clean `10⁻⁵` case and `D_lepton=2·dimG₂`; quark Koide is soft.)

## Honest scope — what this is NOT (yet)

- **It is algebraically equivalent to `t²=½`.**  So as stated it is a **reframing**, not a
  derivation — its value lies in being *type-correct and motivated* (a grade-equipartition /
  critical point), where the symmetry attempts were type-wrong.  To become a derivation it needs a
  **reason the lepton is grade-balanced** (`|ζ|²=1`, off-diagonal norm = diagonal norm).  Candidates:
  (i) **criticality** — grade-balance is the critical point where the operator is equally "local"
  (diagonal) and "hopping" (off-diagonal); (ii) the **Furey idempotent** — the lepton is a
  primitive idempotent, and idempotency may force the balance (idempotents carry the "½" signature:
  the Furey vacuum block `(1+iω)/2` is genuinely idempotent with trace-fraction ½ — to be connected);
  (iii) the **color singlet** — the lepton is the simplest/most-symmetric state.
- **The phase `φ=2/9` is untouched** — grade-balance constrains the scalar/non-scalar *norm*, not
  the internal phase of `ζ`. Still open (and not π-rational, so likely a separate story).

## Net

The multivector reframe produced its first positive result: **a geometric-product grade-balance
condition gives the Koide amplitude `t²=½` exactly** — the right value from the type of structure
the audit said was missing. It is not yet a forced derivation (it is a type-correct reframing of
`t²=½`), but it converts the open question from the vague "why `t²=½`?" into the sharp **"why is
the lepton mass operator grade-balanced (`|off-diagonal| = |diagonal|`)?"** — with concrete
candidate forcings (criticality / idempotent / color-singlet) to test next. First forward motion
after the string of negatives, and it is on the audit-predicted track.
