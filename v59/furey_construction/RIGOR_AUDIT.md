# Rigor audit: derived vs. numerically-relational vs. conjecture

*2026-05-24.  Honest per-claim status, prompted by the question "are these fully derived or
just numerically relational?"  Short answer: **only the algebraic skeleton is theorem-grade;
the physics predictions are empirical/structural matches hanging on a few load-bearing
conjectures.**  The suspicion about `v_Higgs` (and the others) is well-founded.*

## The key distinction

Three very different statuses get blurred when we say "predicts":
- **Theorem** — a Lean-checked identity (e.g. `dimG₂/dimSpin7 = 2/3`, `5 = h∨(Spin7)`).  These
  are *mathematics*; they are certain but say nothing physical by themselves.
- **Numerically relational** — a formula, *found by matching* a structural number to data, with
  **no derived mechanism** for why that formula holds.
- **Conjecture** — an ansatz (often with a structural integer in it) that is not derived at all.

A crucial meta-point: **even the Koide relation `Q=2/3` is not "derived" as physics** — it is an
empirical fact (`~10⁻⁵`) that we *match* to the ratio `dimG₂/dimSpin7`.  "Fully derived" applies
to the math identities; every physical prediction is ultimately an empirical match to a
structural number.  So the right question is not "derived?" but "**(a) how tight is the match,
and (b) is there a mechanism, or just a fitted formula?**"

## Per-claim verdict

| claim | mechanism? | match | verdict |
|---|---|---|---|
| grade structure, `J_c`, lepton=L, dual-Coxeter 5, dim identities | n/a (pure algebra) | exact | **theorem** ✓ |
| Koide `Q=2/3` (lepton) | none (Koide is unexplained) | `10⁻⁵` | tight match, **no mechanism** |
| `φ = Q/3` (lepton phase) | `/3` = sedenion S₃ ✓; magnitude `Q` no | `10⁻⁵` | `/3` structural, **magnitude residual** |
| **quark Koide `Q_d=11/15, Q_u=23/27`** | **selection rule `D_N∈{28,35,63}` UNDRIVEN** | **~0.3%, scale-dependent** | **numerically relational — weak** |
| **`α(M_Z)=25/(324π²)`** | none — `g_W²=5√α` shown *equivalent*, not more primitive (2026-05-24) | `0.03%` | **dimensionless value-match (the `g_W²=5√α` "law" is this value in disguise)** |
| `g_W²=5√α` | the `5` is `h∨` ✓; the `√α` form **excluded as a law** (2026-05-24) | `0.28%` | **NOT a derivable law — = `α(M_Z)` value in disguise** |
| **`v_Higgs=28²·a_l²`** | **none — looks fitted** | `0.07%` | **numerically relational — suspicious** |
| `G_e=(21/16)α²¹`, α(0) | none | varies | **conjecture** |

### Quark Koide (11/15, 23/27) — you're right to doubt it
It is **not derived**.  Two independent weaknesses: (1) the assignment `D_0=28` (lepton),
`D_1=35` (d), `D_2=63` (u) is exactly the **undriven Z₂×Z₂ selection rule** — memory: "the
mapping `N→(Bit-L,Bit-F)` is not derived."  (2) Unlike the lepton case, the *quark* Koide
relation is **scale- and scheme-dependent** (quark masses run; Koide is not RG-invariant) and
`u,d` masses carry large uncertainties — so "0.3%" is a soft, scale-conditioned match, far weaker
than the lepton's `10⁻⁵`.  Status: a suggestive numerical pattern, **not** a derivation.

### α(M_Z)=25/(324π²) — derived, but only *given* a conjecture
The chain is honest arithmetic: `4πα = g_W²·sin²θ_W` (standard SM relation `e²=g_W²sin²θ_W`),
with `sin²θ_W=2/9` (**now derived** from Pati-Salam `(5,2)`) and `g_W²=5√α`, solves to
`α=25/(324π²)`.  The **load-bearing weak link is `g_W²=5√α`** — the `5` is the dual Coxeter
(theorem), but the *relation* `g_W² ∝ √α` has **no Lagrangian / no derivation**.  So `α(M_Z)` is
"derived modulo `g_W²=5√α`," not authoritative.

### v_Higgs=28²·a_l² — the most suspicious, agreed
This ties the **electroweak scale** to the **lepton mass scale** by the pure integer `28²` (=
`dim(L)²`).  There is **no mechanism** — it reads like "we found `v_Higgs/a_l² ≈ 784` and noticed
`784=28²`."  It is also **load-bearing for the "1 input" headline**: that claim *requires*
`v_Higgs` to collapse onto `a_lepton`.  If `28²` is coincidence, `v_Higgs` is a *second* free
scale and the lepton/EW block has **2 scales, not 1**.  This deserves a real derivation or a
demotion.

## Where you're (slightly) off base

**`a_lepton` cannot be "derived authoritatively" — it is dimensionful.**  No pure-mathematics
construction can output a number with units of mass; `a_lepton` *is* the choice of scale (the
unit).  One dimensionful input is unavoidable in any theory.  So the rigorous targets are the
**dimensionless** quantities:
- **`α`** (dimensionless) — a legitimate derivation target; this is where to push (the α
  conjectures).  If `α` becomes structural, the lepton/EW block runs on the *single
  dimensionful* `a_lepton` plus pure structure.
- **dimensionless ratios** like `v_Higgs/a_l²`, `a_lepton/M_Planck`, `a_lepton·G_F^{1/2}` — these
  *can* be structural, and `v_Higgs=28²a_l²` is precisely a claim about one such ratio.

So your instinct (nail down `a_lepton`↔`α`) is right **if** rephrased as: *derive `α` and the
dimensionless ratio linking `a_lepton` to the EW scale.*  You don't derive `a_lepton`; you derive
everything dimensionless around it.

## Priority order for the hard work (load-bearing first)

1. **`v_Higgs=28²·a_l²`** — most load-bearing for "1 input," weakest mechanism.  Either find a
   Higgs-potential / mass-sum-rule derivation, or demote it to coincidence (then count = 2 scales).
2. **`g_W²=5√α`** — load-bearing for `α(M_Z)`, `m_W`, `m_Z`.  Needs a Lagrangian-level reason for
   the `√α` form; the `5` is already secured.
3. **`α` itself** (the α(0)/α(M_Z) conjectures) — the only dimensionless input; if derived, drops
   to "1 dimensionful input."
4. **quark Koide / the Z₂×Z₂ selection rule** — lower priority (soft match, RG-dependent); worth
   honesty-flagging as "pattern, not prediction."

## Honest bottom line

You are correct: **most of the "dependent" relations are numerically relational or
conjecture-dependent, not derived.**  What is genuinely theorem-grade is the *algebraic
skeleton* (grades, dimensions, dual Coxeter, the kernel algebra, `sin²θ_W=2/9` given `(5,2)`,
the sedenion S₃).  The *physical* reductions — quark Koide, `α(M_Z)`, `v_Higgs`, `g_W²=5√α`,
`G_e` — hang on conjectures (chiefly `g_W²=5√α` and `v_Higgs=28²a_l²`) that have **not** been
derived.  The "≈1 free parameter" headline is therefore **conditional on those conjectures**;
rigorously it is "1 dimensionful scale `a_lepton` + `α` + a stack of unproven structural
ansätze."  Working through that stack — starting with `v_Higgs` and `g_W²=5√α` — is exactly the
right next phase, and yes, it is a lot of work.
