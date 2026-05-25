# C1 attempt: `g_W²=5√α` from the v58 `Ω²` quadratic — RESULT (negative, clean)

*2026-05-24.  Per the steer "attempt C1 from v58's `Ω²` quadratic."  Outcome: **the `√α` form
is NOT derivable from `Ω²` (or any gauge Lagrangian) as a dynamical law** — it is incompatible
with the definitional SM relation `e=g_W sinθ_W` except at a single point, and the `Ω²`
mechanism is excluded three independent ways.  What survives is a *value*-match for `α(M_Z)`,
not a derivation.  Computation: `bridge_vhiggs_cl7.py`'s sibling check (inline below).*

## The core obstruction: `g_W²=5√α` is a single value, not a law

The SM relation `e = g_W sinθ_W` is **definitional** (it *defines* the mixing angle); it gives,
with `sin²θ_W = 2/9` (v59-derived),

> `g_W² = 4πα/sin²θ_W = 18π·α`   — **linear in α, holds at every scale** (with running `sin²θ_W`).

The conjecture `g_W² = 5√α` is `∝√α`.  A line and a square-root **meet at exactly one point**:

    18π·α = 5√α   ⟺   √α = 5·sin²θ_W/(4π) = 5/(18π)   ⟺   α = 25/(324π²) = 1/127.91.

(verified — identical to `25/(324π²)`).  So **`g_W²=5√α` is algebraically the single value
`α(M_Z)=25/(324π²)` in disguise**, given the definitional `e=g_W sinθ_W`.  It is *not* an
independent functional law.  The running check makes this explicit:

| `α⁻¹` | `g_W²` (SM `=18πα`) | `5√α` | agree? |
|---|---|---|---|
| 137.0 | 0.4127 | 0.4271 | no |
| **128.0** | **0.4418** | **0.4419** | **yes (EW scale)** |
| 120.0 | 0.4712 | 0.4564 | no |
| 100.0 | 0.5655 | 0.5000 | no |

The two forms coincide only at `α⁻¹≈128` — the EW scale.  As couplings run, the SM linear law
holds everywhere; `5√α` crosses it once.  **`g_W²=5√α` is a single-scale coincidence, not a law.**

## Why the `Ω²` mechanism cannot supply the `√α` (three independent failures)

The proposal (`INTEGRATION_PLAN_v58 §2`) was: "a coupling entering EM linearly enters the
squared/`√` form in the cross-channel," via the quadratic `λΩ²+μ⟨Ω,Ω⟩`.  Steelmanned, this is
the balance `DΩ_EM ~ λ·Ω_weak²` ⇒ `e ~ λ·g_W²` ⇒ `g_W² ∝ √α`.  Each step fails:

1. **Contradicts v58's own derived structure.**  v58 *derives Maxwell* as `DΩ_EM = J_χ` — EM is a
   **primary** sector with its own kinetic term sourced by the chiral current, **not** sourced by a
   weak `Ω²`.  The balance the mechanism needs (EM sourced *by* the weak self-interaction) is the
   opposite of what v58 proved.
2. **v58 has no weak sector and the coefficient is out of band.**  v58 is a classical 2D
   gravity+EM model — there is no `SU(2)_L`, no `Ω_weak`.  And the steelman balance needs the
   self-coupling `λ ≈ 1/h∨ = 0.2` to land the prefactor `5`; v58's **safe band is `|λ| ≤ 0.005`**
   (required for its classical limits).  `0.2` is **40× outside** it — using it would break v58's
   own Newtonian/Maxwell derivations.
3. **Standard gauge embedding gives linear, not `√`.**  In genuine Yang–Mills, embedding
   `U(1)_em ⊂ Spin(7)` gives `e = g_W sinθ_W` (the photon is the Cartan combination), so
   `α = e²/4π ∝ g_W²` — **linear**.  The non-abelian self-interaction `A∧A` (the real "`Ω²`")
   contributes the `g_W²·W⁺W⁻` quartic to the photon, but still through `e ∝ g_W` — giving the SM
   linear relation, never `√α`.  No projection of a Spin(7) connection's `Ω²` yields `g_W² ∝ √α`.

So neither v58's equation nor standard v59 gauge structure produces the `√α` scaling.  The `√` —
the one thing that makes the number work (replacing `√α` by `α` gives `1/11`, not `1/128`) — is
exactly the un-mechanized part.

## What survives (and the honest restatement)

- **Proven and untouched:** `5 = h∨(Spin(7))` (`GaugePrefactorDualCoxeter`) and `sin²θ_W = 2/9`
  (`WeinbergPatiSalam`).  The `5` was never the problem.
- **The real content is a value-form for `√α`:**

  > `√α(M_Z) = h∨(Spin7)·sin²θ_W/(4π) = 5·(2/9)/(4π) = 5/(18π)`  ⇒  `α(M_Z)=25/(324π²)` (0.032%).

  This is clean (built from two proven integers and `4π`) **but it is a match of a structural
  number to the measured `α(M_Z)`, not a derivation** — and it predicts **`√α`, with the square
  root un-derived.**  Stating it as "`g_W²=5√α`" dresses a value-coincidence as a gauge law; the
  honest form is "`√α(M_Z)` equals this structural number," a numerical-relational claim.

## Consequence for the parameter budget

This **demotes** `g_W²=5√α` from `RIGOR_AUDIT`'s "load-bearing conjecture needing a Lagrangian" to
**"not a derivable law — excluded as `Ω²`-mechanism; it is `α(M_Z)`'s value repackaged."**  Hence:

- `α(M_Z)` is **not** "derived modulo `g_W²=5√α`."  The two are *equivalent*; the single
  conjectural content is the dimensionless **value** `α(M_Z)=25/(324π²)` (0.032% match).
- So **`α` remains a genuine conjectural input** (a value-match), exactly as `RIGOR_AUDIT`'s
  "rigorous-only" count had it: the lepton/EW block runs on `a_lepton` **+ `α`**.  The optimistic
  "1 input" picture does **not** get the `α(M_Z)` reduction for free — there is no mechanism, only
  a tight number.
- `m_W, m_Z` (which use `α(M_Z)`, `sin²θ_W`, `v_Higgs`) inherit this: structural up to the one
  `α(M_Z)` value-match.

## Net

`g_W²=5√α` **cannot be derived from v58's `Ω²`** — the `√α` form is incompatible with the
definitional `e=g_W sinθ_W` as a law (single-point coincidence at the EW scale), and the `Ω²`
balance fails on v58's own structure (EM is primary, no weak sector, `λ` 40× out of band) and on
standard gauge embedding (which gives linear `g_W²∝α`).  The result is a clean, machine-anchorable
**negative**: the `5` and `sin²θ_W=2/9` are real; the `√α` is an un-mechanized value-form for
`α(M_Z)`.  *Recommendation:* record `α(M_Z)` honestly as a dimensionless value-match (0.032%), drop
the "gauge-law" framing, and treat `α` as the genuine remaining input.  C2/C3 (per `ALPHA_SCOPING`)
add nothing derivable.  **The `α` derivation target is, on this analysis, closed negative — `α` is
a value-conjecture, not derivable from v58/v59 as currently constituted.**
