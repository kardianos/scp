# Gap G2/G3 — the dimensionless couplings (α and the √α gauge form)

**Folder:** `v59/gaps/alpha_couplings/` · **Date:** 2026-05-25 · **Status:** background + solution-space + verification.

This folder works the two *dimensionless-coupling* gaps of the v59 unified theory. They are
entangled (see §4), so they are treated together.

> **Status tags** (same convention as `UNIFIED_THEORY.md`):
> **[thm]** machine-checkable identity (certain mathematics);
> **[emp≈X]** empirical match of a structural number to data at precision X;
> **[conj]** ansatz with a structural integer, not derived;
> **[free]** genuine input;
> **[dead]** investigated and excluded as a derivation.

---

## 1. G3 — the fine-structure constant α itself

v59 conjectures **two** structural forms for α, at two scales, using two different templates.
Both are conjectures; α is currently a genuine dimensionless input.

### IR form (α at zero momentum, 1/137)
$$-\ln\alpha + 2\alpha \;=\; \frac{\pi^2}{2} \;=\; \frac{8\pi^2}{\dim Cl(3,1)}, \qquad \dim Cl(3,1)=16.$$

- The RHS identity `8π²/16 = π²/2` is **[thm]** (`AlphaCouplingIdentities.ir_rhs_structural`).
- The *template* is a BPST-instanton action `S = 8π²/g²` with `g² = 16 = dim Cl(3,1)`, so
  `α = e^{−S}`. The pure form `−ln α = π²/2` gives `α⁻¹ = e^{π²/2} = 139.05`, off by **1.47%**
  from `137.036` (`verify_forms.py [3a]`). This is the mechanism-backed number. **[conj]**
- The tight match (**3.6×10⁻⁵**, `verify_forms.py [3b]`) requires the **`+2α` term**, which
  `04_alpha_prediction.py` Part 3 **reverse-engineered**: it computes the residual
  `π²/2 + ln α₀ = 0.014559` and scans `{1/70, 1/64, 1/(8π²), α, 2α}` — `2α = 0.014595` fits best.
  The `+2α` has **no instanton justification**; it is a one-parameter fudge. **[conj/fitted]**
- Two foundational holes even in the 1.47% form: (a) *why* `g²=16` is asserted, not derived;
  (b) `π₃(S⁷)=0` ⇒ **no S⁷-based topological instanton exists**, so `S=8π²/g²` is a borrowed
  number with no constructed instanton, no 4D base, no computed topological charge.

### EW form (α at the Z mass, 1/128)
$$\alpha(M_Z) \;=\; \frac{25}{324\pi^2} \;=\; \Big(\frac{5}{18\pi}\Big)^2.$$

- The squared identity `(5/(18π))² = 25/(324π²)` is **[thm]** (`AlphaCouplingIdentities.ew_form_sq`).
- Numerically `α(M_Z)⁻¹ = 127.910` vs PDG `127.951` → **[emp≈0.032%]** (`verify_forms.py [1]`).
- The integers have prior provenance: `5 = h∨(Spin7)` (see §2), `18 = 2·9` where `9` is the
  denominator of the Brannen phase `sin²θ_W = 2/9`, and the `4π` is the SM `e²=4πα`. So this is
  not a free-floating fit — but see the overfitting discussion (§5).

---

## 2. G2 — the √α gauge-coupling form `g_W² = 5√α`

$$g_W^2 = 5\sqrt\alpha, \quad g_R^2 = 5\sqrt\alpha, \quad g_{B-L}^2 = 2\sqrt\alpha,
\quad \tfrac1{g'^2}=\tfrac1{g_R^2}+\tfrac1{g_{B-L}^2}\Rightarrow g'^2=\tfrac{10}{7}\sqrt\alpha.$$

- **The `5` is theorem-grade.** `5 = h∨(Spin(7))`, the **dual Coxeter number** of `so(7)=B₃`,
  the natural β-function normalization (`C₂(adj)=2h∨`). It equals `dim Spin(7) − dim Cl(3,1)
  = 21 − 16` and `N−2 = 7−2`, all coinciding. **[thm]**
  (`GaugePrefactorDualCoxeter.gW_prefactor_is_dualCoxeter_spin7`;
  `AlphaCouplingIdentities.five_eq_dim_diff / five_eq_dualCoxeter / five_readings_agree`.)
- **The `√α` form is a conjecture — and an audit (`04_gW_sqrt_alpha_result.md`) found it is NOT
  a derivable law.** The Killing/harmonic combinations `1/5 + 1/2 = 7/10` and the mixing-angle
  index `2/(5+2·2)=2/9` are clean **[thm]**, but the `√` itself has no Lagrangian origin.

---

## 3. What is theorem-grade vs conjectural (one table)

| object | status | where |
|---|---|---|
| `5 = h∨(Spin7) = 21−16 = 7−2` | **[thm]** | `GaugePrefactorDualCoxeter`, `AlphaCouplingIdentities §1` |
| `8π²/16 = π²/2` | **[thm]** | `AlphaZero.rhs_eq_structural`, `…ir_rhs_structural` |
| `(5/(18π))² = 25/(324π²)`, `√ = 5/(18π)` | **[thm]** | `ScaleBridge.sqrt_alpha_MZ_form`, `…ew_form_sq` |
| `sin²θ_W = 2/(5+2·2) = 2/9`, `cos²=7/9`, `1/5+1/2=7/10` | **[thm]** given (5,2) | `ScaleBridge`, `…§5` |
| `g_W² = 5√α` ⟺ `α(M_Z)=25/(324π²)` (entanglement) | **[thm]** | `…gW_sqrt_alpha_is_alphaMZ` |
| `α(M_Z) = 25/(324π²)` matches PDG | **[emp≈0.032%]** | `verify_forms.py [1]` |
| `−ln α + 2α = π²/2` matches CODATA | **[emp≈3.6e-5]** *(needs fitted +2α)* | `verify_forms.py [3]` |
| `−ln α = π²/2` (pure instanton) | **[conj, emp≈1.47%]** | `verify_forms.py [3a]` |
| the `√α` *scaling* has a Lagrangian origin | **[dead]** | `04_gW_sqrt_alpha_result.md` |
| `g²=16` for the instanton; the `+2α` correction | **[conj/fitted]** | `04_alpha_prediction.py` |
| α is derivable from v59 as constituted | **[dead]** (value-conjecture) | `ALPHA_SCOPING.md` |

---

## 4. The "√α in disguise" finding — why G2 and G3 are entangled

The SM relation `e = g_W sinθ_W` is **definitional** (it *defines* the mixing angle). With
`sin²θ_W = 2/9` it gives the **linear** law `g_W² = 4πα/sin²θ_W = 18π·α`, exact at every scale.
The conjecture `g_W² = 5√α` is `∝√α`. **A line and a square-root meet at exactly one point:**

$$18\pi\,\alpha = 5\sqrt\alpha \iff \sqrt\alpha = \frac{5}{18\pi} \iff \alpha = \frac{25}{324\pi^2} = \frac{1}{127.91}.$$

So `g_W² = 5√α` is **algebraically the single value `α(M_Z) = 25/(324π²)` in disguise**, not an
independent functional law. The two forms coincide only at `α⁻¹≈128` (the EW scale) and diverge
elsewhere (`verify_forms.py [2]` tabulates this). **G2 carries no information beyond the G3 EW
value.** Machine-checked both directions: `AlphaCouplingIdentities.gW_sqrt_alpha_is_alphaMZ` and
`…alphaMZ_satisfies_consistency`. Consequence: closing one closes the other; there is only *one*
dimensionless content here (the EW α value), not two.

---

## 5. The overfitting risk — how many forms can hit α?

Quantified in `overfit_scan.py`. Honest, two-sided result:

- **Within a fixed 2-parameter family, the v59 forms are surprisingly tight.**
  - `α = (p/(mπ))²`, `p≤30, m≤60`: the v59 `(5,18)` is the **only** reduced form beating 0.1%;
    **zero** competitors beat 0.03% in range.
  - `−ln α + c·α = a·π²/b`, `a≤4, b≤12, c∈[−5,8]`: the v59 `(a=1,b=2,c=2)` is the **unique**
    small-integer hit to 0.1%.
  So the *parametrization* is not loose — a modest point in v59's favour that earlier prose
  understated.
- **The real overfitting axis is the choice of TEMPLATE (look-elsewhere across families).** When
  you let yourself pick *which* simple template to use, `overfit_scan.py [D]` finds that **5 of 8**
  arbitrary, physically-unmotivated templates (`i·π²−j`, `i·ln(j)`, `(i/j)·e^π`, …) already reach
  0.03% on `α₀⁻¹` with small integers. So "α has *a* clean structural form" is **cheap**.
- **Therefore precision alone is not evidence.** The discriminator is the **shared provenance** of
  the integers: do `5`, `2/9`, `16` come from the *same* algebra that fixes Koide, `sin²θ_W`, the
  generation count? For the EW form, `5` and `2/9` do (moderate weight). For the IR form, `16`
  does but the `+2α` is free (weak). For `g_W²=5√α`, nothing new beyond the EW value (§4).

---

## 6. Files in this folder

| file | purpose |
|---|---|
| `README.md` | this background (both gaps, theorem vs conjecture, the disguise finding, overfitting) |
| `ALTERNATIVES.md` | solution avenues beyond the current conjectures, with falsification tests |
| `verify_forms.py` | verifies every conjectured form vs CODATA/PDG to full precision (numpy/stdlib) |
| `overfit_scan.py` | quantifies overfitting: within-family uniqueness + cross-template look-elsewhere |
| `AlphaCouplingIdentities.lean` | Lean (Mathlib) formalizing every clean identity; `sorry` only on the genuine open steps |
| `FINDINGS.md` | per-gap verdict (live/dead/open), most promising avenue, what a real derivation needs |
