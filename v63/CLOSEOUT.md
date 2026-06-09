# v62→v63 Closeout — The Number-Type Map of the SCP Phase Sector

**Date**: 2026-05-28
**Scope**: the v62 number-type reframing of the v59–v61 residuals, and the v63
program that tested every route to a *parameter-free* derivation of the Brannen
phase. **Outcome**: the phase sector is fully mapped — what's proven, which routes
are closed, what reduced, and the single irreducible residual (`Q = 2/3`).
**Verification**: 7 self-checking Python scripts (all pass) + 5 Lean modules
(build against the v59 Mathlib; one cited Lindemann `sorry`, the rest `sorry`-free).

---

## 0. The thesis (the organizing discriminant)

Every target quantity is either **algebraic** (a root of a ℚ-polynomial — closed
under `+−×÷`, the output of any dimension/Casimir/character/eigenvalue
construction) or **transcendental** (not). An algebraic machine can only emit
algebraic numbers, so:

> A transcendental target is unreachable by any algebraic construction.

The program's wins are the algebraic column (`784`, `14`, `Q=2/3`, `sin²θ_W=2/9`,
gauge `5/2/9`); its chronic frustration is the transcendental column (the phase
invariant `cos(2/3)`, `α`, `f_g`). The frustration was a category error.

---

## 1. PROVEN (machine-checked / self-verifying)

| # | result | status | artifact |
|---|---|---|---|
| P1 | `cos(2/3)`, `cos(2/9)` transcendental ⇒ the phase `φ=2/9` is producible by **no** algebraic map (character, Casimir, J-eigenphase) | **thm** (1 cited Lindemann `sorry`) | `v62/lean/PhaseNoGo.lean`, `v62/no_go/` |
| P2 | the 01–06 fixed phase is a corollary (amplitude maps can't rotate an algebraic eigenphase to a transcendental); the "07 cut" is barred a priori | **verified** | `v62/no_go/transcendence_nogo.py` |
| P3 | **flat-direction law**: every Z₃-symmetric mass invariant (`Q`, `e₂`, `Σm`) is φ-independent; φ lives only in the `cos(3φ)` class | **thm/sym** | `v62/flat_direction/` |
| P4 | the EW democracy (`S^783`) and the phase are the *same* phenomenon (magnitude fixed, angle flat) | **verified** | `v62/flat_direction/`, `EwVevHome.lean` |
| P5 | the democracy vacuum is **selected** 3 independent ways (Sₙ-fixed-point, max-entropy, explicit-breaking Hessian) | **verified** | `v62/residual_audit/democracy_selection.py` |
| P6 | `φ=2/9` is **not** π-rational; the canonical G₂ character at a triality element is **algebraic** (`χ₇=2`); `Λ⁴(ℝ⁷)=1+7+27` (no 14) | **thm/verified** | `v63/lean/PhaseOffLattice.lean`, `dh_localization.py`, `g2_loop_precheck.py` |
| P7 | **bounds on the phase parameter** `r=c₃/c₆`: hard `r∈[−4,4]` (φ∈[0,π/3]); loop `|φ|≲14α`. Lepton `2/9 > 14·α(M_Z)`; `(2/9)/α(M_Z)>20`; `φ_d<14·α(M_Z)` | **thm** | `v63/lean/LoopBound.lean`, `loop_bound.py` |

---

## 2. CLOSED — every route to a *parameter-free* `2/3` is ruled out

| route | why it fails | artifact |
|---|---|---|
| algebraic invariant / character / J-eigenphase | `cos(2/3)` transcendental (P1) | `PhaseNoGo.lean` |
| `c_eff` = BLV core velocity `v_min/c ≈ 2/3` | parameter-dependent (`0.670` massless vs `0.625` massive), not `2/3` | `ceff_test.py`, `CeffNotStructural.lean` |
| `c_eff` in the Cl(3,1) gravity sector | induced metric is conformal-class only (Urbantke degree-3); every magnitude open/fitted | `conformal_check.py` |
| holonomy / Wilson loop `ξ³` | structural holonomies are π-rational; `Q=2/3` is not (v59 already ran this) | v59 `holonomy_result.md`, `PhaseOffLattice.lean` |
| canonical (DH/Liouville) localization | canonical = commensurate = algebraic; `2/3` needs a generic off-lattice point | `dh_localization.py` |
| density-weighted integration (soliton ρ) | inherits soliton parameter-dependence; extracts but doesn't derive | `integrate_density.py` |
| G₂-equivariant finite-sum loop (Step-6) | the "14" is in `Λ²` not `Λ⁴`; and `N_X` is an overall prefactor that **cancels** in `c₃/c₆` | `g2_loop_precheck.py` |
| parameter-free 1-loop ratio `−c₃/(4c₆)` | `−7.7`/`+38` (out of `[−1,1]`, form-sensitive) → commensurate only | `loop_ratio.py` |

**Unifying dichotomy** (every route converges here): the transcendental phase
`2/3` ⟺ an off-lattice (non-π-rational) angle ⟺ a **free continuous parameter**.
`2/3` is a *bare* radian ratio (not π·rational), the signature of a dynamical
coupling value, not a geometric/algebraic angle.

---

## 3. REDUCED — the partial wins (the residual ledger tightens)

- **Quark Brannen phases → α-downstream.** `φ_d ≈ 14·α(M_Z)`, `φ_u ≈ −10·α(0)` (N
  = `dim G₂`, `dim Spin(5)`), inside the loop ceiling `≲14α`, matching to **<1%**.
  They collapse into `α` + a structural integer — **not independent inputs**.
- **EW democracy `v=784a²` → closable** (P5): selected by Sₙ-equivariance /
  max-entropy / explicit breaking. Removes `v=784a²` from the residual list (mod
  the scale `a`).
- **`f_g ~ α^{21/2}`**: exponent `21/2 = dim Spin(7)/2` — if forced, `f_g` is
  α-downstream, not independent (open: confirm the exponent is derived not fit).

---

## 4. THE one irreducible residual — the lepton coupling magnitude `Q`

The lepton **sector** is maximally structural — `D_lepton = 28 = 2·dim G₂ =
dim Spin(8) = dim Λ²+Λ⁶(ℝ⁷)`, `L=Λ²⊕Λ⁶`, `Q=2/3=dim G₂/dim Spin(7)=14/21`,
`(1−Q)D=28/3`, `784=28²=dim End(L)` — and the phase `φ=Q/3` decomposes:

- **the `/3` is STRUCTURAL** — the sedenion `S₃ ⊂ Aut(𝕊)=G₂×S₃` (a verified order-3
  automorphism = exactly three generations);
- **the magnitude `Q=2/3` is RESIDUAL** — the *coupling* phase `arg(ξ)` (the S₃
  only forces `M` circulant; any `ξ` is allowed), not symmetry-fixed, not a
  holonomy (π-rational obstruction).

The value is fully structural (`(14/21)/3 = 2/9`, exact to `~10⁻⁵`); the mechanism
is absent. The `N≈28=D_lepton` α-reading is a coincidence, not a mechanism (1.5%
vs exact `Q/3`; conflicts with the quark `N=dim G₂` rule).

> **The entire residual reduces to one sharp question:**
> *why is the coupling magnitude `3φ = Q = dim G₂/dim Spin(7) = 2/3`?*
> Its only structural target — the transcendental `cos(2/3)` — matches no candidate
> (Lie/Casimir/√3/character/holonomy). It is mechanism-less in exactly the way the
> Koide relation itself was for decades: a tight (`~10⁻⁵`) law with a structural
> value and no derivation.

---

## 5. Honest limitations

- P1's transcendence of `cos q` (Lindemann–Weierstrass) is a cited `sorry` —
  Mathlib ships only its analytical part. The proof is airtight; the corollaries
  are machine-checked.
- "Closed" routes are closed for a **parameter-free** `2/3`. The phase still exists
  as a dynamical/free input; nothing here makes the theory inconsistent.
- The quark `α`-reduction (§3) is empirical (`<1%` fits with structural `N`), not a
  derived loop calculation. Whether `f_g`'s `21/2` and the quark `N`'s are *forced*
  remains open.
- This is a re-organization of v59–v61, not new physics; its value is the map and
  the closed routes, not a new prediction.

---

## 6. Artifacts & verification

```
v62/  THESIS.md  README.md  verify_all.py(4/4)
      no_go/{NOGO.md, transcendence_nogo.py}            lean/PhaseNoGo.lean (1 cited sorry)
      flat_direction/{FLAT_DIRECTION.md, *_demo.py}     lean/EwVevHome.lean (0 sorry)
      residual_audit/{RESIDUAL_PARTITION.md, democracy_selection.py, phase_dynamical_home.py}
v63/  CEFF_TEST.md ceff_test.py                         lean/CeffNotStructural.lean (0 sorry)
      CONFORMAL_CHECK.md conformal_check.py
      INTEGRATE_BOTH.md integrate_density.py dh_localization.py  lean/PhaseOffLattice.lean (0 sorry)
      LOOP_RATIO.md loop_ratio.py
      G2_PRECHECK.md g2_loop_precheck.py
      LOOP_BOUND.md loop_bound.py                        lean/LoopBound.lean (0 sorry)
      CLOSEOUT.md (this file)
```

- **Python**: 7 scripts, all embedded self-checks pass.
- **Lean**: 5 modules build against the v59 Mathlib (`lake env lean`); only
  `PhaseNoGo.lean` carries a `sorry` (the cited Lindemann–Weierstrass input).

---

## 7. Recommendation (where v64 could go)

The phase sector is exhausted for parameter-free derivation. Productive directions
*not* about the phase:

1. **Confirm the reductions**: is `f_g`'s exponent `21/2 = dim Spin(7)/2`
   *derived*? Are the quark `N = 14, 10` *forced* by the G₂/Spin(5) branching?
   Each closed `−1` to the input count.
2. **The democracy term**: which `Sₙ`-invariant operator the EW potential actually
   contains (P5 shows democracy is selectable; v62 doesn't fix *which* term).
3. **Accept the residual**: treat `Q = dim G₂/dim Spin(7) = 2/3` as a single named
   coupling input (alongside `α` and the scale `a`) — the cleanest honest input
   ledger — and stop searching for a mechanism that the number-type forbids.

The number-type map is the deliverable: it converts a diffuse set of "value
conjectures" into a sharp ledger with **one** residual coupling, `Q = 2/3`.
