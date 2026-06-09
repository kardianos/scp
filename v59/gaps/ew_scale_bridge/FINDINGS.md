# FINDINGS — Gap G1, the electroweak scale bridge `v_Higgs = 28²·a_ℓ²`

**Date:** 2026-05-25 · **Verdict in one line:** the bridge is **NOT yet a derivation** — it is a
*structural conjecture with one theorem-grade ingredient* (the 784-dim space is forced) and two
genuinely-open dimensionful residuals (R1, R2), now sharpened to a single crux: a **rank tension**
between the count (wants full rank 28) and the spectrum (wants rank 3).

---

## 1. Verdict: derivation or conjecture?

**Conjecture, with a real upgrade and a sharper gap.** Status by layer:

| layer | content | status |
|---|---|---|
| the algebraic reframing `v=28²a² ⟺ √v=28a ⟺ Σ√m/√v=3/28` | pure algebra | **[thm]** (Lean: `vHiggs_sqrt_form`, `ratio_form`, `ratio_value`) |
| the empirical match | `√v/a_ℓ = 28.0096` (0.034%); `v/a_ℓ² = 784.54` (0.068%) | **[emp≈3×10⁻⁴]** |
| **the space is 784-dim** | `𝒜_L = ⟨ad so(8)⟩ = End(L) = M₂₈(ℝ)`, dim 28² (Burnside) | **[thm]** — numerically confirmed and shown *generic* (every abs.-irred. so(n) adjoint → dim², so(4) reducible → caught) |
| lepton-specificity | `√v=D·a` works for D=28 only (d/u fail 1.4×/18×) | **[emp/thm]** (only color-singlet lepton's bilinear is wholly in L) |
| **R1: `v = ‖Y‖²_F`** | EW scale = Frobenius² of the L-bilinear (not trace/top-eig) | **[conj]** — *unproved physical identification* |
| **R2: component scale = `a_ℓ`** | so(8) forces all-784-equal, not the *value* | **[conj]** — *unproved dimensionful pin* |
| **25-dir / rank** | 28 light vs 3 generations | **[open]** — the crux (see §3) |

So Option 5 genuinely **upgraded** the headline: `784` is no longer "a noticed coincidence," it
is `dim End(L)`, *forced* once `L=so(8)` acts by its bracket. But "the EW scale is the Frobenius²
over that space, with per-component scale `a_ℓ`" (R1+R2) is **still posited**. Until R1+R2 are
derived, `v_Higgs` is rigorously a **second free scale** and the "1 input `a_ℓ`" headline is
*conditional* on this conjecture (consistent with `RIGOR_AUDIT.md`).

---

## 2. Most promising mechanism

**B + D combined** (see `ALTERNATIVES.md`):

- **B (spectral identity):** frame R1 as the *single* conjecture "**a physical scale = the
  Frobenius² of the relevant mass bilinear.**" This is not ad hoc to the EW sector — the
  **gravity charge** independently fits the *same* Frobenius² form (`Σm = Tr(M†M) = 9Qa²`, over
  the 3-dim generation space), and the two are tied by `Σm_ℓ = (6/784)·v` at 0.07%
  (`gravity_ew_coeff` `[thm]` for the coefficient). One conjecture now covers *both* the EW
  bridge and the gravity charge — a real economy and a passed consistency check.

- **D (give it a home):** identify the v59 physical Higgs (`XiVacuum`, the `ℍ≅ℝ⁴` Mexican hat
  with the proven `{λ,0,0,0}` spectrum = 1 Higgs + 3 eaten `W±,Z` Goldstones) with the **active
  block of the 28×28 bilinear** `Y ∈ End(L)`. This is the only route that connects to v59's
  *existing Lean-proven* Higgs mechanism. It requires building the explicit embedding
  `ℝ⁴ ↪ End(L)` (the deferred Witt/generation map, `SevenDAlgebra.lean`).

The so(8)-invariant potential (route A) is the most *concrete* and was tested directly: it shows
a Frobenius²-hat **can** pin `‖Y‖²_F`, and the sign of the quartic shape term `g` **dynamically
selects** a democratic (full-rank, all-784-comparable) vacuum for `g>0`. That is the closest a
potential gets to R1+R2 — but it leaves the *scale* (`2μ²/λ`) free and produces 28 light
directions, not 3.

---

## 3. The crux (the sharpest result of this session): a rank tension

The two pillars of the bridge demand **opposite ranks** of the vacuum bilinear `Y`:

- **The 784-count (R1, Option 5)** needs a **full-rank, democratic** `Y` (all 28 singular
  values comparable) so that `‖Y‖²_F = 28²·a² = 784 a²`. A full-rank `Y` has **28** light
  directions.
- **The 3-generation spectrum (Brannen kernel)** needs a **rank-3** `Y` (only √m_e,√m_μ,√m_τ
  nonzero). But a rank-3 `Y` has `‖Y‖²_F = m_e+m_μ+m_τ = Σm = 9Qa²` — the **gravity charge**,
  **NOT** `784 a²`.

These cannot both be read off the *same* `Y`. Confirmed numerically (`formalize_bridge.py` §4,
§6): the democratic potential vacuum has eff. rank ~23–28 (→ ~28 light), and **no** single-step
`so(8) → H` breaking leaves exactly 3 light directions (no 25-dim or 3-dim residual fits). So:

> **The EW bridge's `784 = dim(L)²` Frobenius² count and the 3-generation light spectrum are
> structurally in tension on the rank of the vacuum bilinear.** The "25 unaccounted directions"
> is the visible symptom of this tension, not a separate bookkeeping issue.

This is a *constructive* negative: it tells you the missing mechanism must explain how a *single*
`Y` simultaneously (a) has Frobenius² = 784 a² (full ambient weight) and (b) projects to a
rank-3 light spectrum. The natural resolution is **two objects**: an ambient 28-dim vacuum
condensate `⟨Y⟩` (full rank, sets the EW scale via `‖⟨Y⟩‖²_F`) and a *separate* 3×3 Brannen mass
operator (rank 3, sets the spectrum), with the 25 heavy/Goldstone directions living in `⟨Y⟩`'s
ambient. Route C (composite, two 28-legs) or D (v58-vacuum-as-bilinear) are the candidates to
realize this — but it must be **built**, not posited.

---

## 4. What a full derivation needs (checklist)

1. **A dynamical principle forcing the democratic (maximally-mixed) ambient vacuum**
   (`EQUIPARTITION_PRINCIPLE.md` argues this *one* principle would close Koide=2/3 *and*
   `v=784a²` together). Candidate: maximum-entropy / unique `SO(28)`-invariant fixed-norm state.
   Status: principle of *derivable type*, not yet derived from an energy/free-energy extremum.
2. **The identification `v = ‖Y‖²_F` (R1)** as a Lagrangian statement — e.g. the radial mode of
   route D, or the gap-equation vacuum energy of route C, coming out as a Frobenius². The
   spectral-identity framing (B) reduces this to *one* conjecture shared with gravity.
3. **The value pin `component = a_ℓ` (R2)** — tie the 28×28 ambient condensate's per-component
   scale to the 3×3 Brannen `a_ℓ = (Σ√m)/3`. This needs the explicit `3 generations ↪ 28-dim
   ambient` embedding map.
4. **Resolve the rank tension (§3)** — the decisive step. Show one consistent structure (likely
   ambient condensate + 3×3 spectrum, with 25 heavy/Goldstone directions) gives full-rank
   Frobenius² *and* rank-3 light masses. Falsifier: if every consistent embedding forces the
   *vector* (√28) scaling instead of the *bilinear* (28) scaling, the bridge is a coincidence and
   `v` is a second free scale.

---

## 5. Null results recorded (honest)

- **so(8)-invariant Higgs potential does NOT close R1+R2.** It pins `‖Y‖²_F` only to the free
  `2μ²/λ`, and a democratic minimum has ~28 light directions, not 3 (`formalize_bridge.py` §4).
- **Single-step `so(8)→H` breaking does NOT yield "3 light + 25 eaten."** No proper subalgebra
  of `so(8)` has dim 25; maximal proper subalgebras top out at dim 21 (`formalize_bridge.py` §6).
- **Trace and top-eigenvalue readings give the wrong magnitude** (only Frobenius² of a
  *democratic* `Y` gives 784; the singlet gives the excluded 28) — so R1 cannot be weakened to a
  more "obvious" invariant (`formalize_bridge.py` §5).
- **Vector/equipartition reading remains excluded** (`√28 ≠ 28`, `readings_differ` `[thm]`),
  reconfirming `03_higgs_bridge_result.md`.

---

## 6. Bottom line

Option 5 made the **784-dimensional space a theorem**; this session made the **remaining gap
precise**: it is not "why 784" (settled) but the *conjunction* R1 (`v = ‖Y‖²_F`) + R2
(`component = a_ℓ`), and these collide with the 3-generation spectrum in a **rank tension** that
no tested single-mechanism resolves. The most promising path is to read R1 as one
gravity-and-EW-shared conjecture (`scale = Frobenius²`) and give it a Lagrangian home by
identifying the v58 `XiVacuum` Higgs with the active block of `Y ∈ End(L)` — then *build* the
ambient-condensate + 3×3-spectrum structure that resolves the rank tension. Until that is done,
`v_Higgs` is rigorously a **second scale**, and the "1 input" headline stays conditional.
