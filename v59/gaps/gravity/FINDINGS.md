# Gravity gaps — FINDINGS (G8 magnitude, G9 scalar-vs-tensor)

**Date**: 2026-05-25 · **Cluster**: gravity · **Scope**: `/v59/gaps/gravity/`
Verification: `g8_exponent_test.py`, `g9_polarization_test.py`,
`g9_spin2_route_test.py` (all numpy, self-contained, runs clean);
`G8G9_Gravity.lean` (`import Mathlib`, **written not built** this run — the shared
`furey_construction/lean` project must not be built concurrently).

---

## Verdict on G8 — gravity magnitude mechanism

**Is `G_e = (21/16) α²¹` a mechanism, or numerology with a striking exponent?**

> **It is a VALUE conjecture with a theorem-adjacent EXPONENT and a plausible-
> but-unbuilt PRODUCT mechanism. The exponent is real; the prefactor and the
> Lagrangian are not.**

Findings, in order of strength:

1. **The exponent 21 is data-forced [strong, ~thm-adjacent].** The exponent that
   makes the prefactor O(1) is `n = 21.0005 = dim Spin(7)` to **0.0024%**. A
   falsifier scan over `{14, 16, 21, 28}` and `{18..23}` shows **only n=21** lands
   the prefactor in the O(1) window — neighbours are off by `10^-15` to `10^15`. So
   `f_g ~ α^(21/2)` is forced by the data: the V6 `~10^40` overshoot is *exactly* the
   missing `α^21 ≈ 10^-45`.

2. **`α^(21/2)` is the natural value of a "product over 21 generators" [structural,
   LIVE].** Three equivalent realizations give the power 21 cleanly:
   - **determinant** `det(√α · I_21) = (√α)^21 = α^(21/2)` (proved
     `det_scalar_matrix_pow`);
   - **top-form/Haar volume** of the 21-dim Spin(7) algebra, each leg `√α`;
   - **instanton-like** `e^{-S}` with `S = 21·ln(1/α)` — "21 units, one per
     generator." Notably `−ln(α_G(e)) = 103.06` and `21·ln(1/α) = 103.33`, so the
     **entire prefactor is the 0.27 sub-leading correction `ln 1.31`** to a
     21-instanton action.
   This is the structural dual of the **additive** gauge index `g_W² = 5√α`
   (`5 = 21−16`, a *sum* over a subset). Multiplicative-over-all vs additive-over-
   subset is a clean, consistent contrast that *would* explain why gravity is
   uniquely weak.

3. **The mechanism is NOT derived [open].** No Lagrangian maps the graviton coupling
   to that determinant / top-form. The best candidate is a Wess–Zumino-like Haar-
   volume `Λ^21` term on the (odd-dim) Spin(7) group manifold — unbuilt. **Decisive
   subtlety**: only odd-degree top-forms / volume terms can give a 21-form;
   even-degree `Tr(F^k)` characteristic classes cannot. So the falsifier is sharp.

4. **The prefactor is weak [open].** `21/16 = 1.3125` (0.25%) is **beaten by
   17/13 = 1.3077 (0.12%)**; `4/3` is 1.8%. And 21 is over-determined
   (`= 28−7 = 14+7 = 35−14`). So "21/16 = dim Spin(7)/dim Cl(3,1)" is *not* uniquely
   selected. A genuine mechanism must fix the prefactor as a one-loop fluctuation
   determinant — until then it's a free O(1) rational.

**Net G8**: the most law-like gravity number in v59 (exponent forced to dim Spin(7)),
on a plausible determinant/top-form/instanton skeleton, but with no Lagrangian and a
non-unique prefactor. **Mechanism: partial. Honest label: striking-exponent value
conjecture.**

---

## Verdict on G9 — scalar vs tensor (the decisive gap)

**Can the v59 structure carry a spin-2 (`h = ±2`) mode, as LIGO requires?**

> **No — not as currently formulated. The propagating gravity mode is a SCALAR
> (`h = 0`); the internal `so(8)/Λ²` bivector does NOT supply spacetime spin-2.
> The spin-2 *representation* exists in the algebra but cannot reach a spacetime
> tensor without machinery v59 does not have. This is the decisive obstruction,
> and it is fatal for LIGO consistency regardless of the magnitude.**

The chain of results (helicity along z; `g9_polarization_test.py`,
`g9_spin2_route_test.py`):

1. **The v59 gravity mode is `h = 0` [thm].** `□Ω_grav = f_g ρ_grav` with
   `ρ_grav = Tr(M†M)` a **Lorentz scalar** ⇒ scalar carrier ⇒ helicities `{0}`.

2. **The internal index is inert [thm].** Helicity is a *spacetime* little-group
   label; the `so(8)` index is internal and adds **no** spacetime helicity
   (`internal_index_inert`: `n·{0}` has no ±2). So "gravity lives in the 28-dim
   bivector `L`" does **not** give spin-2.

3. **A spacetime bivector is spin-1, not spin-2 [thm].** An antisymmetric 2-form
   carries `{±1}` (`spacetime_bivector_is_spin1`); a massless 4D 2-form has **1 DOF
   and is dual to a scalar** (the "notoph", `h=0`). So even promoting `Ω` to a
   *spacetime* 2-form fails. Only a **symmetric** rank-2 tensor carries `{±2}`.

4. **The spin-2 representation IS in the algebra, but stranded [thm + analysis].**
   `Sym²(so(8) adjoint)` contains symmetric-traceless rank-2 (spin-2) reps
   (verified in the so(4) analog: 3 j=2 multiplets in `Sym²(adjoint)`). But it
   cannot become the LIGO graviton because:
   - **(i) no soldering**: nothing identifies 4 of the 8 internal `so(8)` directions
     with spacetime Lorentz, so the spin-2 rep stays *internal* (a graviton charged
     under so(8), not propagating in spacetime). v59 has `G₂⊂Spin(7)⊂Spin(8)` but no
     vielbein/Plebański soldering form;
   - **(ii) composite/quadratic**: realized as `T_μν = Ω_μ^a Ω_ν^a` it is the gauge-
     field stress tensor — `O(α)`-suppressed, short-range, a *source*, not a
     fundamental massless graviton.

**Net G9**: the algebra *has* the spin-2 rep but no way to put it on spacetime as a
massless field. **The structure cannot currently carry spin-2.** This matches (with
a different algebra) the historical V6/null-rotor cautionary result: long-range
modes were `h = 0, ±1` only — no `h = ±2`.

---

## Most promising avenue

**Prioritize G9 over G8** — a correct magnitude on a scalar force is still ruled out
by LIGO, so G9 is the gate.

The single most promising line is **G9-A: an induced / Plebański metric**. Recast the
OBE as a *constrained 2-form (BF/Plebański) theory* in which `Ω ∈ Λ²` is the
fundamental self-dual 2-form `B` and the metric/`h_μν` is **derived** (this is how GR
*is* a constrained 2-form theory). The make-or-break tests are concrete:
- **soldering**: can a simplicity/Plebański constraint on the so(8) 2-form pick out a
  4D Lorentz frame (a vielbein), mapping the `Sym²(adjoint)` spin-2 rep to a
  *spacetime* symmetric-traceless tensor?
- **DOF count**: does the constrained 2-form propagate **exactly 2** TT DOF (the
  graviton), not 0/1/5?

If that branching + DOF test passes, gravity gets a genuine `h = ±2` mode and G8's
`α^(21/2)` sets its strength — a real unification. If it fails (no soldering, wrong
DOF, or only the `O(α)` composite), then **gravity's spin-2 sector is not inside the
Furey algebra**, and v59 must either add a fundamental soldered `h_μν` (G9-C —
conceding unification: v59 predicts `G_N`, not the graviton) or abandon the LIGO
claim. The branching/DOF computation in a real `so(8)→Lorentz` embedding is the
recommended next concrete step.

For **G8** in parallel: compute the **one-loop fluctuation determinant** around the
candidate 21-generator top-form / instanton and check whether it equals `21/16`
(vs `17/13`). That is the test that would upgrade the prefactor from "free O(1)
rational" to "derived," and with it the whole magnitude from value-conjecture to
mechanism.

---

## File manifest

- `README.md` — precise statement of both gaps, tagged [thm]/[emp]/[conj].
- `ALTERNATIVES.md` — solution space: G8 (det/top-form/instanton/transmutation +
  prefactor-as-determinant); G9 (induced metric, `Sym²` branching, fundamental
  graviton fallback, F-grade route) — each with test + falsifier.
- `g8_exponent_test.py` — exponent pinning, product/det/top-form counting,
  falsifier scan, instanton/transmutation alternatives.
- `g9_polarization_test.py` — helicity decomposition of all candidate carriers;
  internal-index inertness; only sym-TT tensor has h=±2.
- `g9_spin2_route_test.py` — 4D massless-2-form DOF no-go; `Sym²(adjoint)` spin-2
  rep exists but is stranded (no soldering, composite).
- `G8G9_Gravity.lean` — Lean formalization of the structural identities (written,
  not built this run; one flagged `sorry` on the empirical magnitude match).
