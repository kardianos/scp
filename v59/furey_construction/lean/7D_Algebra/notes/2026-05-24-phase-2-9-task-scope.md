# Scope: the full task list for deriving φ = 2/9

**Date**: 2026-05-24
**Purpose**: enumerate every task touching the value `2/9`, with deliverables, dependencies,
and honest difficulty.  This is a planning document, not an implementation.

## The situation (why this is open)

`2/9` appears in **two structurally independent** places in v59:

1. **Mass sector** — the Brannen lepton **phase** `φ = 2/9 rad` in the cyclic mass kernel
   `M = a(I + ξS + ξ̄S²)`, `ξ = t·e^{iφ}` (`KernelEigenvalues.lean`).  Brannen's empirical
   fit is `δ ≈ 0.22222`; v59 elevates it to `φ = Q/3 = 2/9`.
2. **Gauge sector** — the weak mixing angle `sin²θ_W = 2/9` (`ScaleBridge.sin_sq_thW`),
   derived a *different* way via Pati-Salam (`11_u1_y_origin.py`).

**The Lean gap (crux).**  In `LieDimensions.lean`:
```
def brannen_phi : Rat := koide_Q / (n_generations : Rat)   -- φ := Q/3  (a DEFINITION)
theorem brannen_eq_koide_div_gens : brannen_phi = koide_Q / 3 := by unfold ...  -- vacuous
```
So `φ = Q/3` is *assumed by definition*; the `2/9` theorems are arithmetic on it.  What is
**proven**: the amplitude `t² = 1/2 ⟺ Q = 2/3 = dimG₂/dimSpin7` (`BrannenKernel`), and the
*matrix* eigenvalues realize the Brannen amplitudes (`KernelEigenvalues.lam_eq_brannen`).
What is **not proven**: that the mass-operator **phase** must equal `Q/3` (the physical
content `3φ = Q`, the factor 3 = three generations / Z₃).

---

## Tier 0 — Foundational (prerequisite for everything)

- **T0.1  De-circularize the phase.**  Introduce the physical phase `φ_phys` as an
  *independent* parameter of the mass kernel (the argument of `ξ` that fits the charged
  leptons), separate from the *defined* `brannen_phi := Q/3`.  Re-state the goal as a genuine
  theorem `φ_phys = Q/3` (currently false-by-absence).  Until this exists, all downstream work
  decorates a definition.
  *Deliverable*: a `def phi_phys` + a `theorem (to prove) phi_phys = brannen_phi`.  *Difficulty*: low (statement), but it exposes that the content is unproven.
- **T0.2  Audit empirical vs structural vs derived.**  One table: `φ_measured ≈ 0.22222`
  (empirical, `EmpiricalAgreement.brannen_phase_agrees`), `φ = Q/3` (conjectured-structural),
  derived = ∅.  *Difficulty*: trivial; clarifies the target.

## Tier 1 — Mass-sector derivation: *why* `3φ = Q` (the hard core)

Per project rule, try ≥3 mechanisms before declaring open.  Candidates:

- **T1.1  J∘Z₃ alignment (leverages the completed C-program).**  The phase is `arg(ξ)`, and
  `ξ` rotates in the plane of the now-pinned color complex structure `J_c = γ₀γ₅ ∈ Λ²`
  (`ColorSU3`).  The Z₃ generation shift `S` (`CyclicShift`, `ω = e^{2πi/3}`) acts on the
  three generations.  Conjecture: the relative alignment of the `J_c`-plane and the
  `S`-eigenbasis fixes `φ`.  *Needs*: the generation-space ↔ minimal-left-ideal embedding (the
  "exact Witt map" deferred in `SevenDAlgebra.lean` header).  *Difficulty*: high; most likely
  to be the *right* mechanism if any is.  *Depends on*: C complete (done), Witt-map refinement.
- **T1.2  Vacuum/extremization (analog of `XiVacuum`).**  `XiVacuum` fixed the amplitude
  `|ξ|²=1/2` by a Higgs-type potential minimum.  Is there a phase potential `V(φ)` (e.g. from a
  Z₃-symmetric term) extremized at `φ = Q/3`?  *Deliverable*: a potential whose critical point
  is `2/9`.  *Difficulty*: medium-high; risk that any such `V` is reverse-engineered (must be
  structurally motivated, not fitted — the rigor lesson).
- **T1.3  Second invariant / phase quantization.**  Find an algebraic condition on the 3
  charged-lepton masses, *independent* of Koide `Q`, that is equivalent to `φ = Q/3` and has a
  structural reading (a second symmetric-function identity à la Koide).  *Deliverable*: identify
  the invariant, check it numerically, then derive.  *Difficulty*: medium; good "Kepler" step
  even before a mechanism.
- **T1.4  Further methods if 1.1–1.3 stall**: phase from a spectral/η-invariant asymmetry of the
  kernel; phase from a G₂ holonomy/associator angle; phase as the Berry phase of the Z₃ orbit.
  *Difficulty*: speculative; enumerate, test cheaply, skip with notes.

## Tier 2 — Gauge-sector `2/9` (`sin²θ_W`): formalize + close the "2"

This half is *more tractable* — much is already structural in Python.

- **T2.1  Formalize the Pati-Salam derivation of `sin²θ_W = 2/9`.**  Chain (from
  `11_u1_y_origin.py`): `g_W² = 5√α` (have: `GaugePrefactorDualCoxeter`, 5 = h∨(Spin7)),
  `g_{B-L}² = 2√α`, `g'² = (10/7)√α`, then `1/g'² = 1/g_W² + 1/g_{B-L}²` and
  `sin²θ_W = g'²/(g_W²+g'²) = (1/5+1/2)⁻¹·… = 2/9`.  *Deliverable*: a Lean theorem deriving
  `2/9` from the three coupling ratios (currently `ScaleBridge.sin_sq_thW` just *defines* `2/9`).
  *Difficulty*: low–medium (rational arithmetic + the Pati-Salam relation as hypotheses).
- **T2.2  Derive the "2" in `g_{B-L}² = 2√α`** (the missing structural input, analogous to the
  `5`).  Candidates flagged in Python: `2` = the B-L abelian factor / the silent-direction
  count / the SU(2) rank / the "Brannen-phase numerator".  *Deliverable*: identify `2` as a
  specific invariant and prove the coupling relation.  *Difficulty*: medium-high; this is the
  gauge-sector analog of the open mass-phase mechanism.

## Tier 3 — The unification link (is the mass phase = the gauge angle?)

- **T3.1  Identity vs coincidence.**  Decide whether `φ_Brannen = sin²θ_W` (both `2/9`) is a
  genuine identity or numerical accident.  If genuine, find the **shared** structural object and
  prove both `Q/3` and `g'²/(g²+g'²)` reduce to it.  *Difficulty*: high; this is the real prize
  (would make one derivation serve both sectors).  *Depends on*: T1 and T2.
- **T3.2  Decompose `2/9` into `(numerator 2)/(denominator 9 = 3²)`.**  Show the same `2` and
  `9 = generations²` (or `= 63/7`, `= 10/45`, …) appear in both the mass `Q/3` and the gauge
  `sin²θ_W`.  Cross-link with `ScaleBridge` (already notes `2/9 = 10/45` from Pati-Salam and
  `2/9 = (dimG₂/dimSpin7)/3`).  *Difficulty*: medium; partly a bookkeeping/identity task.

## Tier 4 — Quark extension

- **T4.1  Quark phases.**  `05_quark_sector.py` fits up/down Brannen phases empirically
  (`φ_up ≈ −2.02`, `φ_down ≈ +0.110`).  With structural `Q_d = 11/15`, `Q_u = 23/27`
  (`Predictions`), test whether `φ_d = Q_d/3`, `φ_u = Q_u/3` (or another relation), then derive
  via the Tier-1 mechanism.  *Difficulty*: medium; first a numeric check (cheap), then derive.
- **T4.2  CKM/Cabibbo tie-in.**  `tan θ_C ~ √(m_d/m_s)` and `ScaleBridge.cabibbo_seven` — does
  the phase mechanism feed the mixing angles?  *Difficulty*: high; likely after T4.1.

## Tier 5 — Integration / verification (after any of the above lands)

- **T5.1**  Wire new theorems into `AxiomCheck`; update `EmpiricalAgreement` to point at the
  *derived* (not assumed) phase.  *Difficulty*: low.
- **T5.2**  Honest status update: after each task, restate what remains empirical (target: after
  Tier 1, the only lepton-sector input is the overall scale `a_lepton`).

---

## Dependency graph / suggested order

```
T0.1 ─┬─> T1.1 (needs Witt map)        ┐
      ├─> T1.2 / T1.3 / T1.4           ├─> T3.1 ─> T3.2
T2.1 ─┴─> T2.2                          ┘            │
T1.* + T4(numeric) ───────────────────> T4.1 ─> T4.2 │
                                                      └─> T5
```

**Most tractable first**: T0.1 (de-circularize) → T2.1 (formalize the existing Pati-Salam
`2/9`) → T1.3 (find the second invariant, a cheap "Kepler" check) → T4.1 (numeric quark check).
**Highest value / hardest**: T1.1 (J∘Z₃) and T3.1 (the unification), the genuine open physics.

**Honest expectation**: T0.1, T2.1, T4.1(numeric) are landable.  The *physical* derivation of
`3φ = Q` (T1.x) and the "2" (T2.2) are open research; they may yield a clean mechanism or may
join the documented null results — either outcome is a real result.
