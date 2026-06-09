# Gap G7 — the origin of the lepton mass-phase magnitude

**Folder:** `v59/gaps/lepton_phase/` · **Date:** 2026-05-25 ·
**Status tags:** [thm] machine-checked Lean · [emp≈X] empirical match at precision X ·
[conj] structural ansatz, not derived.

This is the background / problem statement for the "gap attack" on G7. It states
the gap precisely, separates what is **structural** from what is the **residual**,
and lists the **ruled-out** mechanisms with their Lean theorems. The proposed new
(non-geometric) avenues are in `ALTERNATIVES.md`; the verdict is in `FINDINGS.md`.

---

## 1. The gap, stated precisely

The charged-lepton masses are reproduced by the **Brannen circulant mass kernel**
on the 3-generation flavour space:

> `M = a (I + ξ S + ξ̄ S²)`,  `ξ = t·e^{iφ}`,  `S` = the 3-cycle (Z₃ generation shift, `S³=I`).

Its eigenvalues are `√m_k = a(1 + 2t·cos(φ + 2πk/3))`, `k=0,1,2`. There are two
dimensionless shape parameters: the **amplitude** `t = |ξ|` and the **phase**
`φ = arg ξ`. Empirically both are fixed by the single Koide ratio `Q`:

| parameter | relation | lepton value | precision |
|---|---|---|---|
| amplitude | `t² = (3Q−1)/2` (this **is** Koide: `Q=(1+2t²)/3`) | `t²=1/2` | [thm] for the algebra; [emp≈10⁻⁵] for the fit |
| **phase** | `φ = Q/3` (the law studied here) | `φ = 2/9` | [emp≈7×10⁻⁶, 0.89σ] |

Equivalently, the law locks the phase to the modulus: `9·arg(ξ) = 1 + 2|ξ|²`.

**The decomposition of G7.** The factor `Q = 14/21 = 2/3 = dim G₂ / dim Spin(7)`.
The law `φ = Q/3` factors into two pieces with different status:

- **The "/3" is structural [thm-verified].** It is the generation count = the
  sedenion `S₃ ⊂ Aut(𝕊) = G₂ × S₃` automorphism `ψ` (order 3). `ψ` is a genuine
  sedenion automorphism (`ψ³=I`, verified on all 16×16 basis products) acting as a
  120° rotation of the two Im-octonion blocks (`sedenion_s3.py`,
  `sedenion_S3_analysis.md`).
- **The magnitude `Q = 2/3` is the residual.** It is the *same* `Q` as the Koide
  amplitude (`t²=1/2 ⟺ Q=2/3`). The phase invariant `cos(3φ) = cos(2/3)` is what
  needs an origin. The factor `9 = 3·3`: one 3 from `Q=(1+2t²)/3`, one from `φ=Q/3`.

So **G7 = "why is the coupling magnitude exactly `3φ = Q = 2/3`?"** — a mass-sector
input, decoupled from the (explained) generation count.

---

## 2. Why the magnitude `Q` is genuinely extra content [thm]

The Koide ratio is **phase-independent**: `Q(a,t,φ) = Q(a,t,φ')` for every phase
(`BrannenPhase.Q_phase_independent` [thm]). The phase enters observables **only**
through `cos(3φ)`, via the third power sum
`Σ s_k³ = 3a³(1 + 6t² + 2t³·cos 3φ)` (`BrannenPhase.sum_s_cube` [thm]). Nothing in
the kernel forces `3φ = Q`; the relation is a separate empirical fact, as tight as
Koide itself (both `m_τ`-limited, `~0.9σ`; `LeptonPhaseEmpirical.phase_as_tight_as_koide`
[thm], `generation_count_pinned` [thm]).

**The genuine invariant is `cos(2/3)`, not the rational `2/9`.** The masses fix the
phase only up to the generation `S₃` (`φ → φ + 2π/3`, `φ → −φ`); the only invariant
is `cos(3φ) = cos(2/3) ≈ 0.78589` (`PhaseAmbiguity.invariant_is_cos_two_thirds`
[thm]). A naive PDG fit returns `φ ≈ 2.3166 = 2/9 + 2π/3`, a full Z₃-shift away —
so "`φ = 2/9`" is a principal-branch *representative*, and a structural derivation
must target the **transcendental** `cos(2/3)`, not the convention-nice `2/9`.

---

## 3. RULED-OUT mechanisms (the design space already closed)

The geometric route is **dead**. These exclusions are machine-checked.

| # | ruled-out mechanism | Lean theorem | why it fails |
|---|---|---|---|
| E1 | **π-rationality / any geometric angle.** A holonomy / rotation / root angle is a rational multiple of π. | `PhaseExclusions.phase_not_pi_rational`, `koide_not_pi_rational` [thm] | `φ=2/9` and `Q=2/3` are rational *in radians*, hence **not** π-rational. So the value is not any rotation/holonomy/root angle. |
| E2 | **"nice fraction of π"** `φ = π/n`. | `PhaseExclusions.phase_ne_pi_div_nat` [thm] | corollary of E1: `φ ≠ π/n` for all `n`. |
| E3 | **naive Z₃ vacuum-alignment** (lowest invariant phase potentials `cos 3φ`, `cos 6φ`). | `BrannenPhase.z3_potential_does_not_select_2_9`, `PhaseExclusions.cos6_potential_does_not_select_2_9` [thm] | neither harmonic has a critical point at `2/9` (`sin(2/3)>0`, `sin(4/3)>0`). |
| E4 | **the Weinberg-cosine identity** `cos 3φ = cos²θ_W = 7/9`. | `PhaseAmbiguity.phase_invariant_ne_cos_sq_thetaW` (via `cos_two_thirds_gt_seven_ninths`) [thm] | `cos(2/3) > 7/9` strictly; gap ≈ 0.008 ≈ 100× the data precision. The two `2/9`'s share only the *value* `Q/3`, not the cosine. |
| H | **inter-generational holonomy** (the Brannen-ring Wilson loop `arg ξ³ = 3φ`, the v58 connection `ω ∝ ∇log MM̃`). | `holonomy_result.md` + `koide_not_pi_rational` [thm] | the loop flux is a *free* Aharonov-Bohm phase: the v58 connection is an exact 1-form (zero holonomy); discrete/center fluxes are π-rational; the associator gives integers; color-Z₃ fixes the lepton singlet. No structural curvature forces it to `Q`. |
| S | **sedenion S₃ symmetry fixing the value.** | `sedenion_S3_analysis.md` [computed] | `ψ`'s intrinsic phases are `{0, ±2π/3}` (π-rational); the Brannen `φ` is the *coupling* phase `arg ξ`, not an eigenphase. S₃-covariance forces only that `M` is circulant — **any `ξ` is allowed**. |
| M | **max-mixing over the generation structure fixing `t²`.** | `MaximalMixingKoide.lean` HONEST-SCOPE note + `08_generation_map_result.md` [retracted] | the `Aut(𝕊)=G₂×S₃` generation symmetry leaves `t²` free (max-mixing over the *generation* structure gives `t²∈{0,1}`, never `1/2`). So `(D−dimG₂)/D` is a **so(8)-grade dimensional coincidence, not a generation-level mechanism** (see §4 — it is still the best magnitude candidate). |

**Net of §3:** the magnitude must come from a **non-geometric** algebraic / spectral
/ extremal relation that produces the *transcendental* `cos(2/3)`. This is the only
structural avenue left, and it faces the "not-a-nice-number" hurdle.

---

## 4. The one surviving thread (and its honest limit)

The maximal-mixing / G₂-content reading gives the **magnitude** `Q=2/3` cleanly,
*non-geometrically*, as a Lie-**dimension** ratio (no angle):

> `t²_sector = (D − dim G₂)/D`, with the `G₂ = Aut(𝕆)` core inert. For the lepton,
> `D = dim L = 28 = 2·dim G₂`, so the non-G₂ fraction is **exactly 1/2** ⇒ `Q=2/3`.
> (`MaximalMixingKoide.lepton_half_from_double_core`, `koide_from_maximal_mixing` [thm
> for the arithmetic]; the quark sectors give the "ugly" `t²=3/5, 7/9` ⇒ `Q=11/15,23/27`.)

This is a **theorem for the arithmetic** and grounds `dim G₂ = 14` (octonion
derivation algebra). Its **honest limit** (retraction M above): the *bridge* — that
the maximally-mixed L-vacuum's non-G₂ weight **is** the Brannen `t²` — does not close
at the generation level. So `Q=2/3` is established as the **so(8)-grade dimension
ratio that numerically equals `t²=1/2`**, not as a derived dynamical `t²`.

This thread is pursued and extended in `ALTERNATIVES.md`.

---

## 5. Files in this folder

- **`README.md`** (this file) — background, structural vs residual, ruled-out list.
- **`ALTERNATIVES.md`** — non-geometric candidate origins (max-mixing/G₂-content,
  Casimir ratios, equipartition/extremal, the skewness/"second-Koide" reading, the
  √3 of `ψ`), each with a test and a falsifier.
- **`search.py`** — numerical search of v59 algebra invariants (Casimirs, dimension
  ratios, max-mixing weights, extremal points) for `2/3` and `cos(2/3)`; sensitivity
  of `Q=2/3` vs nearby fractions. Run: `python3 search.py`.
- **`second_invariant.py`** — the skewness reframing: `cos(3φ)` = standardized third
  moment of `√m`; tests whether any clean `f(Q)` (non-radian) fits the data. Run it.
- **`LeptonPhaseMagnitude.lean`** — Lean: the maximal-mixing Koide derivation
  (re-stated for this gap) + the clean identities found (equipartition `2t²=1` at
  `t²=1/2`; the skewness identity `cos 3φ = (4/(2t)³)·⟨x³⟩`). `import Mathlib`;
  open steps flagged with `sorry`. **Written, not built this run** (shared project).
- **`FINDINGS.md`** — verdict: is there a live non-geometric origin for `Q=2/3`?
