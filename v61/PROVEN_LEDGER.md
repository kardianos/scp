# v59–v61 Proven Ledger — What Is Actually Established

**Date**: 2026-05-27 · **Purpose**: a cross-version audit answering one question —
*of everything claimed across v59, v60, v61, what has **actually been proven**, and in
what sense?* This is an adversarial re-read of the closeouts, the Lean modules, and the
verification scripts, separating genuine theorems from arithmetic identities, empirical
matches, curve-fits, textbook recitation, and hypothesis-as-theorem.

It is deliberately stricter than the per-version closeouts. Where it disagrees with a
closeout's framing, the disagreement is stated.

## Verification basis and its limits

- **Lean was NOT rebuilt for this audit** — `lake`/`lean`/`elan` are not in PATH in the
  audit shell. The v59 closeout claims `lake build` succeeds (8278 jobs) with `AxiomCheck`
  clean; v60/v61 claim 7+5 modules build against v59 Mathlib. Those build claims are
  **taken on report, not re-verified here**. What *is* assessed below is the **content of
  the theorem statements**, which is independent of whether the prover ran.
- A machine-checked theorem proves *its statement*. It does **not** prove the physical
  interpretation attached to the statement. Most genuine math here proves an algebraic
  fact; the identification with the electron, the Higgs vev, or gravity is conjecture
  carried on top.

**Status tags** (extending the v59 closeout set):
**[thm✓]** substantive theorem · **[arith]** true but content-free (arithmetic / dimension
count) · **[emp≈X]** empirical match at precision X (a fit, not a proof) · **[fit]** a free
exponent/constant solved to match data · **[txt]** standard textbook result, true
independent of this theory · **[hyp→]** stated as `hypothesis ⟹ conclusion` with the
physics in the unproven hypothesis · **[selfcon]** true by construction / consistency
check · **[null]** tested and found false/absent.

---

# PART I — MATH: genuinely proven

## M-A. Substantive theorems (real algebra / representation theory)

These would stand in a math paper with no physics attached.

- **M-A1. `Cl(7)_even` grading and isomorphism** **[thm✓ / arith for the sum]**
  `Cl(7)_even ≅ Cl(6) ≅ ℂ⊗ℍ⊗𝕆`, graded `Λ⁰⊕Λ²⊕Λ⁴⊕Λ⁶`. The decomposition is genuine;
  the headline `1+21+35+7 = 64` part is arithmetic.
  *Files*: `v59/furey_construction/lean/SpinDimension.lean`, `Predictions.lean`.

- **M-A2. Blade square-sign dichotomy** **[thm✓]** — *the load-bearing theorem*
  A simple even k-blade squares to `(−1)^{k(k+1)/2}·𝟙`, so `Λ²,Λ⁶` are complex structures
  (`B²=−𝟙`) and `Λ⁴` is an involution (`B²=+𝟙`). Genuine induction + ring algebra,
  representation-independent.
  *File*: `v59/furey_construction/lean/BladeSquareSign.lean`
  (`prod_sq`, `even_grade_complex_structure_dichotomy`).

- **M-A3. Color `su(3)` closure and lepton-`J` pinning** **[thm✓, partly by `decide`]**
  Eight explicit generators close to `su(3)` (28 brackets); a skew complex structure
  commuting with the color Cartans is pinned on the lepton singlet to `±[[0,1],[−1,0]]`
  (Schur-type). Closure is verified by exhaustive `decide` over explicit 8×8 matrices —
  true, but a case-bash rather than a conceptual proof.
  *Files*: `ColorSU3.lean` (`colorSU3_closes`, `colorInvariant_pins_lepton`),
  `7D_Algebra/LeptonComplexStructure.lean`.

- **M-A4. Koide algebraic identity** **[thm✓]**
  For `M = a(𝟙+ξS+ξ̄S²)`, `S³=𝟙`: eigenvalues `λ_k=a(1+2t cosθ_k)` give
  `Q ≡ Σm/(Σ√m)² = (1+2t²)/3`, and `Q=2/3 ⟺ t²=1/2`. Clean trig/algebra
  (the `Σcosθ_k=0`, `Σcos²θ_k=3/2` identities). **Proves the formula, not that `t²=1/2`
  is physical.**
  *File*: `BrannenKernel.lean` (`Q_value`, `koide_iff_constraint`).

- **M-A5. `dim End(so(8)) = 784` via Burnside/Jacobson** **[thm✓ + num✓]**
  The `so(8)` adjoint is absolutely irreducible, so the algebra it generates is all of
  `M₂₈(ℝ)`, dimension `28²=784`. Confirmed numerically by Gram–Schmidt on `ad(so(8))`
  (and shown to *detect* reducibility: `so(4)` gives 18<36). A genuine theorem. **It does
  not establish that the EW scale equals `784a²`** — that identification is conjecture (see
  P-emp and N4).
  *Files*: `v59/synthesis/obe_options_2_5.py` (num✓); `EwScaleBridge.lean` (`endL_dim`).

## M-B. True but content-free (arithmetic / dimension counts)

Correct, machine-checked, and **physically empty** — the physics is entirely in *why the
integers appear*, which is not proven by these:

- `5 = 21 − 16` (`gauge_index_additive`), `784 = 28²`, `1+21+35+7 = 64`,
  `(28a)² = 784a²`, vacuum-manifold dim `= 783`, `sin²θ_W = 2/9` *as a number*.
  All discharged by `decide`/`norm_num`/`ring`. Calling `5 = 21−16` a "structural
  identity" overstates it: it is arithmetic, not a derivation of why 5, 21, 16 enter
  physics.

## M-C. Genuine negative / structural results (proven — high value)

- **M-C1. `g_W² = 5√α` is not an independent law** **[null, thm-backed]**
  Given the SM-definitional `4πα = g_W² sin²θ_W` with `sin²θ_W=2/9`, the form collapses to
  the single value `α(M_Z)=25/(324π²)`. It carries no information beyond `α`.
  *Files*: `v59/gaps/alpha_couplings/verify_forms.py`, `AlphaCouplingIdentities.lean`.

- **M-C2. The Koide phase is not a geometric angle** **[thm✓]**
  `2/9` and `2/3 rad` are proven **not π-rational** (`phase_not_pi_rational`,
  `koide_not_pi_rational`); `cos(2/3)` is transcendental and equals no Lie/Casimir/√3
  quantity (best near-miss 11/14, off 1.7×10⁻⁴ ≫ data). Rules out the holonomy/Berry-phase
  mechanism.
  *Files*: `v59/gaps/lepton_phase/`, `PhaseExclusions.lean`, `PhaseAmbiguity.lean`.

- **M-C3. `g² = h∨√α` fails for the strong sector** **[null, thm-backed]**
  `g_s² = 3√α ≈ 0.26` vs measured `≈1.48`; a `√α` form cannot produce an O(1) coupling.
  *File*: `v59/gaps/frontier_sectors/`, `FrontierSectors.lean` (machine-checkable falsifier).

---

# PART II — PHYSICS: what is actually established

**Headline: almost no physics is *proven*.** The physical content is empirical matches,
self-consistency, and textbook GR. There is exactly one solid physics *conclusion*, and it
is negative.

## P-emp. Empirical numerical matches — conjectures with a fit, NOT proofs **[emp]**

Striking, worth keeping, but each is a formula chosen to hit data, subject to
look-elsewhere (v59 N10: 5 of 8 arbitrary templates already hit `α` to 0.03%):

- Koide `Q=2/3` (≈10⁻⁵), Brannen `φ=2/9` (≈7×10⁻⁶), `sin²θ_W=2/9`, `m_W`/`m_Z`
  (0.02–0.04%), `v_Higgs=784a²` (0.07%). The v59 audit's own tags are `[emp]`/`[conj]`,
  and it states the rigorous input set is **`a_ℓ + α`**, not `a_ℓ` alone.

## P-fit. `G_e = (21/16)α²¹` is curve-fitting, not derivation **[fit]**

The exponent was obtained by *solving* `n = log(α_G/(21/16))/log(α)` against CODATA — it is
fit to the answer (the script says so: "a VALUE conjecture, not a derivation"). The match
to `dim Spin(7)=21` is post-hoc, and 21 is over-determined (`=28−7=14+7=35−14`).
*File*: `v59/synthesis/gravity_magnitude_test.py`.

## P-selfcon. Gravity-sector "results" are construction or textbook

- **charge/mass = 1 universal** **[selfcon]** — true because `ρ_grav ≡ Tr(M†M) = Σm` is
  *defined* as the mass, not a derived equivalence principle.
  *File*: `v59/synthesis/gravity_charge_test.py`.
- **`Φ~1/r`, force `~1/r²`** **[txt]** — standard potential theory (massless kernel +
  nonzero monopole). Correct, not specific to this theory.
  *File*: `v59/synthesis/obe_radial_test.py`.
- **v61's "three classic GR tests"** (Schwarzschild, deflection `4GM/b`, perihelion
  43″/cy) **[txt]** — these are textbook GR *sourced by* `M=Σm`. Any metric theory yielding
  Schwarzschild reproduces them; they validate *consistency given a Schwarzschild metric*,
  not a theory-specific prediction. v61's own "honest limitations" concede this: "they
  confirm the gravity sector *is* GR with the right source, not a new gravitational
  dynamics."
  *Files*: `v61/01_curved_backreaction.py`, `v61/05_perihelion.py`.

## P-solid. The one solid physics conclusion is NEGATIVE **[null]**

Within this algebra, gravity carried by the `so(8)/Λ²` bivector is a **Lorentz scalar
(helicity 0 only)** and **cannot supply the `h=±2` tensor mode LIGO observes**; the
internal index is inert under the spacetime little group. This (v59 N3) is a genuine
rigorous result and outranks the magnitude work — a perfectly tuned scalar force still
fails LIGO.
*File*: `v59/gaps/gravity/g9_polarization_test.py`.

---

# PART III — Claimed-but-NOT-proven (where the closeouts oversell)

- **"G9 resolved" (v60/v61)** **[hyp→]** — `G9_resolved_if_soldering` is literally
  `IF SolderingExists THEN h=±2`, and `SolderingExists` is an **unconstructed** `Prop`
  (`∃ _ : Unit, True`). The soldering that does the physical work is *assumed*, not built.
  G9 (P-solid) is **reframed, not resolved**.
  *File*: `v60/lean/G9Soldering.lean`.
- **"GEN3 derives Koide Q=2/3" (v60)** — a potential was **reverse-engineered** to have the
  Brannen vacuum; the Hessian/spectrum is then a consistency check on a posited `V`. The
  constraint `(1−Q)D=28/3` is imposed to match the empirical `Q`, not derived.
  *Files*: `v60/lagrangian/14_spectrum.py`, `15_selection_rank.py`.
- **"R1 home for `v=784a²`" (v61 GEN3)** **[hyp→]** — the O(784)-symmetric hat is
  degenerate (vacuum manifold `S^783`), so the democratic vacuum is **not selected**;
  v61's own text downgrades this to "home + isolated conjecture." `R1_frobenius_is_scale`
  carries the physics in its hypothesis `hR1 : v_phys = fro2 Y`.
- **Quark Koide (`Q_d=11/15`, `Q_u=23/27`)** **[null as prediction]** — RG-convention
  artifact; at a single scale the gap is 2–5%, and 12 simple rationals fall in each RG band
  (v59 N2).
- **`α` value, Brannen `φ=2/9`, gravity magnitude `f_g`** — all four "residual conjectures"
  are inputs/values, explicitly **not derived** (v59 N10, and the v60/v61 closeouts agree).

---

# PART IV — One-line verdict

| | Genuinely proven | Claimed but not proven |
|---|---|---|
| **Math** | Blade-sign dichotomy (M-A2); `Cl(7)_even` grading (M-A1); color `su(3)` + lepton-`J` pinning (M-A3); Koide formula `Q=(1+2t²)/3` (M-A4); `dim End(so(8))=784` (M-A5); negatives M-C1–3 | that the integers 5/21/16/784/64 *carry* their physical meaning (M-B) |
| **Physics** | the gravity-is-scalar / no-`h±2` result (P-solid) — and it is negative | Koide/φ/Higgs/`sin²θ_W`/`m_W,Z` matches (P-emp, = conjecture+fit); `α²¹` (P-fit); "GR tests passed" (P-selfcon/txt, theory-independent); "G9 resolved" and "R1 home" (hyp→); `α`, `f_g`, quark Koide |

**Net**: the durable output is a small set of clean **algebra theorems** (M-A, M-C) plus a
set of striking-but-unexplained **numerical coincidences** anchored on two genuine inputs,
`a_ℓ` and `α`. The v59 audit is honest about this; the **v60/v61 closeouts oversell** —
"derives Q=2/3" means a potential was fitted to that vacuum, "passes all GR tests" means
standard GR was re-derived, and "G9 resolved" means resolved *conditional on an unproven
existence claim*. No physics prediction has been *proven*; the strongest physics result is
the negative LIGO obstruction.
