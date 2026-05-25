# φ = 2/9 program — progress (Tiers 0–4)

**Date**: 2026-05-24.  Plan: `2026-05-24-phase-2-9-task-scope.md`.  All modules build,
0 sorry / 0 native_decide, axioms = standard trio (or fewer).

## DONE

### Tier 0 — de-circularize (`BrannenPhase.lean`)
- `Q_phase_independent`: the Koide ratio `Q(a,t,φ)` is the **same for every φ** — Koide fixes
  the amplitude `t` (`t²=1/2 ⟺ Q=2/3`) but says *nothing* about φ.  So `φ = Q/3` is genuine
  extra content, not a corollary of `Q = 2/3`.  (This exposes that `LieDimensions.brannen_phi
  := Q/3` was a *definition*.)
- Power-sum moments: `sum_s = 3a`, `sum_s_sq = 3a²(1+2t²)` (both φ-free), and
  **`sum_s_cube = 3a³(1+6t²+2t³·cos 3φ)`** — the phase enters observables **only through
  `cos 3φ`** (the three generations / `Z₃` 120°-spacing force the third harmonic).  Bundled in
  `moments`.
- **Target sharpened to `3φ = Q`** (`target_three_phi_eq_Q`): the third-harmonic angle equals
  the Koide ratio.  Both sides dimensionless; the "factor 3" = #generations is now explicit.

### Tier 2 — gauge-sector 2/9 (`WeinbergPatiSalam.lean`)
- `sin2_thetaW_eq` / `sin2_from_coeffs`: **derive** `sin²θ_W = 2/9** from the Pati-Salam
  matching `1/g'² = 1/g_W² + 1/g_{B-L}²` with `(cW,cBL) = (5,2)·√α` (√α cancels):
  `g'² = (10/7)√α`, `sin²θ_W = (10/7)/(5+10/7) = 10/45 = 2/9`.  Previously
  `ScaleBridge.sin_sq_thW` just *defined* `2/9`.
- `cW_eq_dualCoxeter`: the `5` is `h∨(Spin(7)) = 7−2` (from `GaugePrefactorDualCoxeter`).

### Tier 3 — unification (`TwoNinthsUnification.lean`)
- `mass_eq_gauge`: the mass route `Q/3` and the gauge route `sin²θ_W` give the **same** `2/9`.
- Decompositions: `2/9 = 14/63 = dimG₂/(dimSpin7·3)` (mass) `= (9−dimImO)/9 = (9−7)/9` (gauge),
  bundled in `two_ninths_both_ways`.

### Tier 4 — quark check (numeric + `BrannenPhase.quark_phase_not_Q_div_three`)
- **`φ = Q/3` is lepton-specific.**  Fitted quark phases (`05_quark_sector.py`):
  `φ_d ≈ 0.1102` vs `Q_d/3 = 11/45 ≈ 0.244`; `φ_u ≈ −2.02` vs `Q_u/3 = 23/81 ≈ 0.284` — both
  far off (proved as rational gaps `>0.1` and `>2`).  Only the lepton matches.  This **constrains
  Tier 1**: the mechanism must be lepton / color-singlet specific — consistent with `J_c`
  living on the lepton singlet `{0,7}` (`ColorSU3`), not a universal Brannen relation.

## OPEN

### Tier 1 — *why* `3φ = Q` (the mechanism)
- **T1.2 (phase potential) ruled out in its simplest form** (`z3_potential_does_not_select_2_9`):
  a `Z₃`-symmetric potential is built from `cos 3φ`; its critical points are `3φ ∈ πℤ`
  (`φ ∈ {0, π/3, …}`), and `V'(2/9) ≠ 0` (since `sin(2/3) > 0`).  So naive vacuum alignment does
  **not** select `2/9`.
- **T1.1 (J∘Z₃)** remains the only live candidate: `φ = arg ξ` with `ξ` in the pinned
  `J_c = γ₀γ₅` plane (`ColorSU3`), the `Z₃` shift `S` on generations.  Needs the
  generation↔minimal-left-ideal (**Witt**) map — deferred in `SevenDAlgebra.lean`.  Not
  formalized; would have to tie an angle in the `G₂ ⊂ Spin(7)` geometry to `Q = 14/21`.
- **Honest status**: `3φ = Q` is numerically exact at fit precision and has the clean
  third-harmonic reading, but no derivation; the simplest dynamical mechanism fails.  It may be a
  coincidence at `m_τ` precision or await the `J∘Z₃` geometry.

### Tier 2.2 — the "2" in `g_{B-L}² = 2√α`
- The gauge `2/9` derivation rests on `cBL = 2` (analog of the derived `cW = 5`).  Candidates
  (B-L abelian factor / silent-direction count / Brannen numerator) not yet pinned.

### Tier 3.1 — deep identity vs coincidence
- Whether `φ_Brannen` (an angle) and `sin²θ_W` (a mixing ratio) being the *same* `2/9` is a deep
  identity (one shared origin) or numerical accident is unresolved.  Value-level equality and two
  independent structural decompositions are proven (`TwoNinths`); a *shared mechanism* is not.

## Negative results (`PhaseExclusions.lean`) — narrowing the design

Proving candidate mechanisms FALSE (the project's null-result method):

- **(E1) `phase_not_pi_rational`** — `2/9` (rational radians) is **not a rational multiple of
  π** (via `irrational_pi`).  So the phase is **not a geometric / holonomy / rotation / root
  angle** (those are π-rational).  ⇒ The whole "tie a `G₂⊂Spin(7)` geometric angle to `Q`" class
  of Tier-1 mechanisms is **ruled out**: the phase must be a *ratio* (like `Q/3`) inserted as a
  radian value, not an angle read off a rotation.  This is the strongest constraint — it
  redirects T1.1 away from holonomy/geometry toward a pure ratio identity.
- **(E2) `phase_ne_pi_div_nat`** — `φ ≠ π/n` for every `n` (corollary).
- **(E3) `cos6_potential_does_not_select_2_9`** (+ the `cos3φ` case in `BrannenPhase`) — neither
  of the two lowest Z₃-invariant harmonic potentials has a critical point at `2/9`; naive vacuum
  alignment (T1.2) fails.
- **(E4) `gauge_cBL_pinned`** — given `cW = 5`, the value `sin²θ_W = 2/9` forces **`cBL = 2`
  uniquely** (`cBL/(5+2cBL) = 2/9 ⟺ cBL = 2`); `1,3,4` give `1/7, 3/11, 4/13`.  So the open
  Tier-2.2 "2" is *not* a free parameter — it is pinned by `2/9` and the dual-Coxeter `5`.  The
  remaining work is to *interpret* the `2`, not to fit it.
- **(E5) `mass_gauge_distinct_reductions`** — the mass `2/9 = 14/63` reduces by `7 = dimImO`, the
  gauge `2/9 = 10/45` by `5 = h∨(Spin7)`: different pre-reduction integers, so a "single shared
  integer drives both" identity (T3.1) is constrained — evidence against a trivial deep identity.

## DECISIVE negative result (`PhaseAmbiguity.lean`) — the phase is not what it seemed

Computing the Brannen phase from PDG masses returns **`φ ≈ 2.3166 rad = 2/9 + 2π/3`** — a full
generation `Z₃` shift from "2/9".  Reason, now proven:

- **`s_shift` / `s_conj`**: `φ → φ+2π/3` (generation `Z₃`) and `φ → −φ` (conjugation) *permute*
  the Brannen amplitudes, so the **whole `S₃` leaves the masses invariant**.
- **`moments_Z3_invariant` / `phase_2_9_not_unique`**: hence `φ = 2/9` and `φ = 2/9 + 2π/3`
  give *identical* charged-lepton masses.  **"φ = 2/9" is a choice of representative, not a
  physical observable.**  (Reducing the fit mod `2π/3` returns `0.22223` — matching `2/9` to the
  `m_τ` precision; so Brannen's "2/9" is just the principal-branch representative.)
- **`invariant_is_cos_two_thirds`**: the *only* phase quantity the masses fix is
  `cos 3φ = cos(2/3) ≈ 0.7859`.  A structural derivation must therefore explain the (not nice)
  `cos(2/3)`, **not** the branch-/unit-dependent rational `2/9`.  Sharp contrast with the
  amplitude: `t² = 1/2` is `S₃`-invariant (it comes from Koide `Q`) and genuinely solid; the
  phase is convention-fixed.
- **`phase_invariant_ne_cos_sq_thetaW`** (rigorous, via `sin(1/3) < 1/3` ⇒ `cos(2/3) > 7/9`):
  the natural Tier-3.1 deep identity `cos 3φ = cos²θ_W = 7/9` is **FALSE** (gap ≈ 0.008, ~100×
  the data precision).  So the two `2/9`'s do *not* unify through their cosines — they coincide
  only as the *values* `Q/3` and `sin²θ_W`, reinforcing the coincidence reading (E5).

**Design impact**: the mass-phase target is no longer "derive 2/9" but "derive `cos(2/3)`,
lepton-only, non-geometric (E1), with the obvious gauge identity excluded".  This is a hard
constraint — it removes the niceness that motivated the conjecture and leaves a transcendental
target.

## POSITIVE knowledge (`LeptonPhaseEmpirical.lean`) — the agreement is Koide-tight

Computing from PDG-2024 masses (phase reduced to its principal `Z₃` branch):

| relation | value | |dev| | significance |
|---|---|---|---|
| Koide `Q = 2/3` | `0.66666051` | `6.2e-6` | `0.91σ` |
| phase `φ = 2/9` | `0.22222963` | `7.4e-6` | `0.89σ` |

- **`phase_as_tight_as_koide`**: `|2/9 − φ_meas| < 10⁻⁵` *and* `|2/3 − Q_meas| < 10⁻⁵` — the phase
  relation `φ = 2/9` holds at the **same `~10⁻⁵ / 0.9σ` precision as the Koide relation itself**
  (both `m_τ`-limited).  So this is NOT a loose `~0.22`-cluster coincidence; it is a tight,
  Koide-class numerical regularity.  (Corrects the earlier over-dismissal.)
- **`generation_count_pinned`**: among `φ = Q/n`, only `n = 3` matches (`Q/2, Q/4, Q` miss by
  `> 10⁻²`).  With `Q = dimG₂/dimSpin7 = 2/3`, the divisor is **pinned to the number of
  generations** → `φ = Q/(#generations) = 2/9`.

So the lepton sector's two dimensionless shape parameters are both, empirically, functions of
the single structural number `Q`: `t² = (3Q−1)/2` (Koide) and `φ = Q/3` — each verified at
`~10⁻⁵`.  `φ = Q/3` is a genuine **"second Koide relation."**

## Balanced verdict

**Positive**: `φ = Q/3 = 2/9` is a real, tight (Koide-class, `~10⁻⁵`) empirical regularity, with
the divisor uniquely fixed to the generation count.  **Negative**: it is *not* a kernel identity
(`Q_phase_independent`), *not* geometric (E1), *not* the Weinberg angle (`phase_invariant_ne…`),
lepton-specific (quarks fail), and its branch-invariant `cos(2/3)` matches no structural
candidate.  So it stands exactly where Koide stood for decades: an empirically precise relation
awaiting a mechanism — but now with the *easy* mechanisms rigorously excluded, so the eventual
explanation (if any) must be the `J∘Z₃` ideal structure, not geometry/potential/gauge-identity.

## Net
The bookkeeping/value tiers (0, 2, 3, 4) are done and rigorous; the gauge `2/9` is now derived,
not assumed, and the lepton target is sharpened and shown to be independent content.  The two
genuinely physical unknowns — the mass-phase mechanism (`3φ=Q`, T1.1) and the gauge "2" (T2.2) —
remain open, with the simplest mass mechanism (T1.2) rigorously excluded.
