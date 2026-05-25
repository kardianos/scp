# Does anything FORCE grade-balance (t²=½)? — the three candidates + all combinations

*2026-05-24.  Per request, test whether the Koide grade-balance `t²=½` (`10_…md`) is *forced* by
any of three candidate conditions — (1) Furey idempotent, (2) criticality, (3) color singlet — in
isolation and in every combination (1&3, 2&3, 1&2, 1&2&3).  Careful, step-by-step.  **Result: NONE
force `t²=½`.** Grade-balance gives the right value but remains a type-correct *reframing*, not a
forced derivation.  One residual thread survives (C1 via the still-missing generation map).*

## Operational definitions (made precise so each can pass or fail)

- **C1 (Furey idempotent):** the lepton mass operator is (twice) a primitive idempotent,
  `M/a = 2P`, `P²=P` — equivalently the off-diagonal `ζ=ξS+ξ̄S²` satisfies `ζ²=I`.
- **C2 (criticality):** `t²=½` is an extremum / special point of a natural functional of the mass
  spectrum (entropy, purity, det, spread, …).
- **C3 (color singlet):** the lepton is the SU(3)_c singlet (color-invariant mass operator).

*(A representation note found en route: in the 8×8 irrep the pseudoscalar `ω₇=±I`, so grade-`k`
and grade-`(7−k)` blades are proportional — grade-projection reversion is ill-defined there.  All
tests below are done in the clean generation space + analytic internal facts, unaffected.)*

## Isolation tests

**C1 — NO.**  `ζ² = ξ²S² + ξ̄²S + 2t²I`; exact `ζ²=I` requires `2t²=1` **and** `ξ²=0`
(impossible).  So `M/a=2P` is **unachievable** for any `(t,φ)`.  The *closest* point
(`min‖ζ²−I‖`) is at **`t²=⅓`** (analytic), **not ½**.  The scalar-part-only version `⟨ζ²⟩₀=1`
does give `t²=½` — but that is exactly the grade-balance restated (**circular**).  The genuine
idempotent "½" lives **internally** (`(1+u)/2`, `u²=−1`, has reversion-norm ½) — but its link to
the *generation* amplitude `t²` is the **still-missing generation map**.

**C2 — NO.**  Scanning natural functionals at `φ=2/9`, the interior extrema fall at
`t²≈0.54` (√m-entropy), `0.87` (m-entropy), `0.95` (purity), none for det/spread — **none at ½**.
No natural spectral functional is critical at grade-balance.

**C3 — NO (vacuous).**  Color acts on the internal 8-spinor; the color `Z₃` fixes the lepton
singlet for *any* internal content (`ColorSU3.lean`), and the generation amplitude `t²` is
orthogonal to color.  Color-invariance is satisfied for **every** `t²` — it constrains nothing.

## Combinations

Because **C3 is vacuous for `t²`**, it adds nothing in any combination:

| case | reduces to | forced `t²` | verdict |
|---|---|---|---|
| **1&3** | C1 | min-defect ⅓ | **NO** |
| **2&3** | C2 | extrema 0.54/0.87/0.95 | **NO** |
| **1&2** | criticality *of* the idempotent-defect = `min‖ζ²−I‖` | **⅓** | **NO** |
| **1&2&3** | 1&2 | ⅓ | **NO** |

The one non-trivial joint condition — **1&2**, "be as close to idempotent as possible" (the
natural marriage of C1's measure and C2's criticality) — extremizes at **`t²=⅓`**, the analytic
minimum of `‖ζ²−I‖²=3(2t²−1)²+6t⁴`.  Not ½.

## What the grind found

1. **None of {1, 2, 3} or their combinations force `t²=½`.**  Grade-balance gives the right value
   (`10_…md`) but is **not forced** by idempotency, criticality, or color — it stays a
   type-correct *reframing* of `t²=½`.
2. **The recurring value when a condition actually bites is `t²=⅓`, not ½** (the idempotent-closest
   point).  `t²=⅓` ⟺ the mass operator nearest to twice-an-idempotent — and note `1−Q_lepton=⅓`
   and the user's `−1/3` thread.  The idempotent structure *prefers* ⅓, not the Koide ½.
3. **C3 (color) is provably vacuous for the magnitude** — color ⊥ generation (reinforces that the
   color-`3` and generation-`3` are distinct, cf. Task #1; an early *negative* indicator for it).
4. **The single surviving thread:** the internal idempotent genuinely carries a "½" (reversion-norm
   of `(1+u)/2`), and C1 *would* force `t²=½` **iff** the generation↔internal map carries that
   internal ½ to the generation amplitude.  That map is the same one unbuilt in `08_…md` — so C1 is
   not fully closed, only **conditional on the missing map**.

## Honest net

After testing symmetry (magnitude-blind, `09_…md`) **and** the multivector grade structure via
idempotent / criticality / color (this note), **nothing in the available octo-space structure
forces `t²=½`.**  Grade-balance is the right-typed, right-valued condition, but it is not *derived*
— and the natural forcings prefer `⅓` over `½`.  The weight of evidence now points to **`t²=½`
(Koide) being a genuinely external magnitude**, *unless* the one untested link — the generation↔
internal map carrying the idempotent's internal ½ — rescues C1.  That map is therefore the single
remaining make-or-break: build it (closes C1, derives Koide) or show it cannot carry the ½
(closes the question negatively, `t²` external).
