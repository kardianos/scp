# α scoping audit — derivable vs. asserted, and what v58/v59 can supply

*2026-05-24.  Same honest treatment just applied to `v_Higgs` (`integration_v58/
03_higgs_bridge_result.md`), now for the `α` program **before** committing to a derivation.
Sources audited: `lean/{AlphaZero,GaugePrefactorDualCoxeter,ScaleBridge}.lean`,
`03_alpha_from_u1.py`, `04_alpha_prediction.py`.  Numbers re-verified.*

> **UPDATE (2026-05-24): C1 was then attempted and came back NEGATIVE**
> (`integration_v58/04_gW_sqrt_alpha_result.md`).  The `√α` form is **not a derivable law** —
> given the definitional `e=g_W sinθ_W` it is the single value `α(M_Z)=25/(324π²)` in disguise
> (linear `g_W²=18πα` meets `5√α` at one point), and the `Ω²` mechanism is excluded 3 ways
> (v58 makes EM primary; no weak sector + `λ` 40× out of band; standard embedding gives linear
> `g_W²∝α`).  So the recommendation below ("attack C1") is **discharged negative**: all three
> α-conjectures are now mechanism-poor, and **`α` is a genuine value-conjecture**, not derivable
> from v58/v59 as constituted.

## There are THREE separate α-conjectures (not one)

They live at different scales and use **different templates** — conflating them inflates the
headline.

| # | conjecture | scale | template | precision (re-verified) |
|---|---|---|---|---|
| **C1** | `g_W²=5√α ⇒ α(M_Z)=25/(324π²)=(5/18π)²` | `M_Z` (~1/128) | **algebraic** (`√α`) | `α(M_Z)`: **0.032%** vs PDG |
| **C2** | `−ln α + 2α = π²/2` ⇒ `α(0)` | `0` (1/137) | **instanton-exp** | **3.6×10⁻⁵** *(but see below)* |
| **C3** | `G_e=(21/16)α²¹` ≈ `exp(−21π²/2)` ⇒ electron grav. coupling | — | instanton-exp²¹ | **0.25%** vs `Gm_e²/ℏc` |

## C1 — `g_W²=5√α` (→ α(M_Z)).  Best target; one clean gap.

- **Proven:** `5 = h∨(Spin(7))` — the **dual Coxeter number** of `so(7)=B₃`, a genuine Lie
  invariant (`GaugePrefactorDualCoxeter.gW_prefactor_is_dualCoxeter_spin7`), and the natural
  β-function normalization (`C₂(adj)=2h∨`).  The arithmetic `α(M_Z)=25/(324π²)` then follows
  from `4πα=g_W²·sin²θ_W` (SM tree) + `sin²θ_W=2/9` (derived) — all machine-checked
  (`ScaleBridge`).  *The `5` is NOT the gap.*
- **The gap (single, sharp):** the **`√α` scaling itself** — why `g_W² ∝ √α` (not `∝α`, not a
  constant).  `GaugePrefactorDualCoxeter` is explicit: "does NOT derive … the `√α` scaling …
  still require[s] the gauge-embedding Lagrangian."
- **v58/v59 home:** the **most plausible of the three.**  v58's gauge field *is* the bivector
  (grade-2) connection `Ω`; the proposal (`INTEGRATION_PLAN_v58 §2`) is that `√α` comes from the
  **quadratic self-interaction `λΩ²+μ⟨Ω,Ω⟩`** (a coupling entering linearly in one channel
  appears as `√` in the cross-channel).  *Unworked, but concrete and falsifiable* — there is an
  actual equation (`⟨DΩ+λΩ²+μ⟨Ω,Ω⟩⟩` ) with an actual `Ω²` term to project.
- **Verdict:** **derived modulo the `√α` form; the `5` is secured.**  Closing `√α` closes
  `α(M_Z)`, `m_W`, `m_Z`.  Highest leverage, most tractable.

## C2 — `−ln α + 2α = π²/2` (→ α(0)).  The tight match is partly **fitted**.

- **The mechanism (instanton):** `04_alpha_prediction.py` proposes a **BPST instanton** on the
  Furey constraint surface `G₂\Spin(7) ≈ S⁷`, action `S=8π²/g²` with `g²=16=dim Cl(3,1)`, giving
  `S=8π²/16=π²/2` and `α=e^{−S}`.  `π²/2 = 8π²/16` is machine-checked (`AlphaZero.rhs_eq_structural`).
- **BUT the mechanism-backed form is only 1.4%:** the pure instanton gives `−ln α = π²/2`, i.e.
  `α⁻¹ = e^{π²/2} = 139.04` vs `137.036` — **1.46% off.**
- **The tight `3.6×10⁻⁵` requires the `+2α` term, which is REVERSE-ENGINEERED.**
  `04_alpha_prediction.py` **Part 3** literally computes "what correction makes π²/2 exact?"
  (`= 0.0146`) and scans candidates (`1/70, 1/64, 1/8π², α, 2α`) — `2α` fits.  So `+2α` is a
  **post-hoc fudge with no instanton justification.**  The "4×10⁻⁵ / 0.004%" headline (memory,
  `RESUME`) rides entirely on it.  *(This is the α analog of the v_Higgs "equipartition is wrong"
  finding: a tight number propped by an unmotivated step.)*
- **Two further foundational holes in even the 1.4% form:** (a) **why `g²=16=dim Cl(3,1)`** —
  asserted, not derived; (b) **`π₃(S⁷)=0`** (cf. the v6/gravity work, memory `Avenue C`) ⇒ there
  are **no topological instantons with `S⁷` as base** — so `S=8π²/g²` is the BPST *number*
  borrowed, with **no constructed instanton, no 4D base, no computed topological charge.**
- **v58/v59 home:** v59 *has* the `S⁷` surface but **no dynamics/instanton sector**; v58 is
  classical 2D with no instantons.  **Neither can currently host the mechanism.**
- **Correction to flag:** `AlphaZero.lean`'s comment "actual gap ≈2.4×10⁻³" is **wrong** — the
  true gap (with `+2α`) is `3.6×10⁻⁵`.  (The relation is tighter than the doc says, but only via
  the fitted term.)
- **Verdict:** **numerically striking but mechanism-poor.**  Honest claim = "instanton template
  gives α(0) to 1.4%; the `0.004%` needs a fitted `+2α`."  Lower priority than C1.

## C3 — `G_e=(21/16)α²¹` (→ electron gravitational coupling).  Numerology.

- Matches `Gm_e²/ℏc = 1.75×10⁻⁴⁵` to **0.25%**, with `21=dim Spin(7)`.  An alternative form
  `G_e=e^{−21π²/2}` (`04` Part 4, "21-fold instanton") is in the same ballpark but **not
  identical** — two inequivalent parametrizations of one number.
- **No mechanism, wild exponent (`α²¹`), no home in v58 or v59.**
- **Verdict:** pure numerical match.  Flag as pattern, do not pursue.

## What v58/v59 can actually supply (the honest map)

| conjecture | needs | v58 can give? | v59 can give? |
|---|---|---|---|
| **C1 √α** | a Lagrangian where `g_W² ∝ √α` | **maybe** — the `Ω²` quadratic term in the bivector connection (concrete, unworked) | the `5=h∨` (done); `sin²θ_W=2/9` (done) |
| **C2 instanton** | a grounded instanton, `g²=16`, no `+2α` | no (classical 2D, no instantons) | the `S⁷` surface only; **no dynamics, `π₃(S⁷)=0`** |
| **C3** | anything | no | no |

## Recommendation (mirrors `RIGOR_AUDIT` priority)

1. **Attack C1 (`g_W²=5√α`).**  Only one gap (the `√α`), the `5` is proven, and it has a
   *concrete v58 mechanism to test* (the bivector-connection `Ω²` quadratic).  High leverage:
   closing it closes `α(M_Z)`, `m_W`, `m_Z`.  This is the real "derive α" target.
2. **Demote C2's headline.**  Re-state α(0) honestly: instanton template → 1.4%; the `0.004%`
   is a fitted `+2α`.  Fix the `AlphaZero.lean` "2.4×10⁻³" comment.  Pursue only if a grounded
   instanton (real base, charge, `g²=16` reason) can be built — currently blocked by no v58/v59
   instanton sector and `π₃(S⁷)=0`.
3. **C3: flag as numerology, drop.**

**Bottom line:** the α program is **one tractable conjecture (C1, gap = the `√α` form, with a
concrete v58 home) plus two mechanism-poor numerical matches (C2 partly fitted, C3 numerology).**
The dimensionless-`α` endgame runs through **C1**, exactly as `RIGOR_AUDIT` ordered.
