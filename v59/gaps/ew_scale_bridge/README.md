# Gap G1 — the electroweak scale bridge  `v_Higgs = dim(L)²·a_ℓ²`

**Folder:** `v59/gaps/ew_scale_bridge/` · **Date:** 2026-05-25
**Status tags:** `[thm]` machine-checkable identity · `[emp≈X]` empirical match at precision X ·
`[conj]` ansatz with a structural integer, not derived · `[open]` unresolved.

> This is the single most load-bearing conjecture in the v59 theory. It ties the **electroweak
> scale** to the **lepton mass scale** by one pure integer, `28² = 784`. If that integer is a
> derivation, the charged-lepton+EW+Higgs block runs on **one dimensionful input** `a_ℓ`; if it
> is a coincidence, `v_Higgs` is a **second free scale**. `RIGOR_AUDIT.md` calls it "the most
> suspicious number." This README states the bridge precisely, records what is established, and
> defines the open residuals R1, R2, and the 25-direction problem.

---

## 1. The bridge, stated precisely

Let `L = Λ²⊕Λ⁶ ⊂ Cl(7)_even` be the 28-dim "lepton grade" (the skew = `so(8)` part; the proven
complex-structure/mass-bearing grade — `BladeSquareSign`, `LeptonRealityForcing` `[thm]`). Let
the charged-lepton √-mass quantum (Brannen scale) be

```
a_ℓ = (√m_e + √m_μ + √m_τ) / 3 = 0.56021525  GeV^(1/2)   [from PDG masses]
```

The conjectured bridge, in three algebraically equivalent forms (`vHiggs_sqrt_form`, `ratio_form`
`[thm]` for the equivalences):

```
   v_Higgs = dim(L)²·a_ℓ² = 28²·a_ℓ² = 784·a_ℓ²        [emp≈0.068%  (squared form)]
⟺  √v_Higgs = dim(L)·a_ℓ = 28·a_ℓ                       [emp≈0.034%  (sqrt form)]
⟺  (Σ√m_ℓ)/√v_Higgs = N_gen/dim(L) = 3/28              [emp≈0.034%]
```

Numerically (computed in `formalize_bridge.py`, PDG `v = 246.21965 GeV`):

| quantity | value | law wants | gap |
|---|---|---|---|
| `√v/a_ℓ` | 28.0096 | 28 = dim(L) | **0.034%** |
| `v/a_ℓ²` | 784.54 | 784 = dim(L)² | **0.068%** |
| `(Σ√m)/√v` | 0.1071062 | 3/28 = 0.1071428 | **0.034%** |

Equivalently a lepton **Yukawa sum rule** `Σ√y_ℓ = (3/28)·2^(1/4)` `[emp≈0.03%]`
(`HiggsVevReframe.yukawa_sum_form` `[thm]` for the algebra).

**Downstream (all hang on this `v`):** `m_W = ½g_W v` (0.04%), `m_Z = (3/√7)m_W` (0.02%),
`m_H = √(7/27)·v` (0.07%). And the gravity bonus `Σm_ℓ = (6/784)·v` (0.07%,
`gravity_ew_coeff` `[thm]` for the coefficient).

---

## 2. What Option 5 settled: the 784-space is *forced* `[thm]`

The earlier worry was "we noticed `v/a_ℓ² ≈ 784` and that `784 = 28²`." Option 5
(`synthesis/obe_options_2_5.py`, re-derived and stress-tested in `formalize_bridge.py` §2)
**upgrades 784 from a coincidence to a dimension theorem**:

- Octonion **non-associativity** forces `|Φ|²` onto the *associative left-action algebra*. For
  `L = so(8)` the product is the Lie bracket, so the relevant operator algebra is
  `𝒜_L = ⟨ad_X : X ∈ so(8)⟩`.
- `ad(so(8))` is the **absolutely irreducible** 28-dim adjoint. By **Burnside / Jacobson
  density**, the unital associative algebra it generates is the *full* matrix algebra:

  ```
  𝒜_L = End(L) = M₂₈(ℝ),     dim 𝒜_L = 28² = 784.          [thm]
  ```

- **Numerically confirmed and shown generic** (`formalize_bridge.py` §2): the defining rep
  (8×8) generates only `8²=64` (the wrong reading); the adjoint (28×28) generates exactly `784`.
  The same `dim² = full End` holds for *every* absolutely-irreducible `so(n)` adjoint tested
  (n=3,5,6,7,8 → 9,100,225,441,784); `so(4)` (whose adjoint `= su(2)⊕su(2)` is **reducible**)
  generates only 18 < 36 — a sanity check that the test detects reducibility. So `784 = dim
  End(L)` is the generic Burnside statement, *forced* once `L = so(8)` acts by the bracket.

- Lean: `endL_dim` / `endL_dim_eq_dimL_sq` (`EwScaleBridge.lean`) assert `dim_ℝ M₂₈(ℝ) = 784 =
  dim(L)²` via `Module.finrank_matrix` `[thm, written-not-built]`.

**The crucial caveat, stated honestly:** Option 5 derives that the *space the bridge lives in*
has dimension 784. It does **not**, by itself, say the EW scale equals a Frobenius² over that
space, nor that each component has magnitude `a_ℓ`. Those are R1 and R2.

---

## 3. The open residuals

### R1 — why `v_Higgs = ‖Y‖²_F` `[conj]`

The bridge needs `v` to be the **Frobenius² norm** of the L-grade mass bilinear `Y` (a 28×28
matrix), not its trace, top eigenvalue, or nuclear norm. `formalize_bridge.py` §5 shows that for
a *democratic* `Y` (every entry `= a`) **only** the Frobenius² gives 784:

| reading | democratic `Y` | so(8)-singlet `Y = a·δ` | bridge wants |
|---|---|---|---|
| Frobenius² `‖Y‖²_F` | **784 a²** | 28 a² (the excluded √28) | 784 a² |
| trace `Tr Y` | 28 a | 28 a | — |
| top singular value | 28 a | 1·a | — |

So R1 is really a **conjunction**: (i) the vacuum is *democratic on L* (all 784 components
comparable), and (ii) the physical `v` is read as the *Frobenius²* (= the L²/energy norm of a
bilinear). The Frobenius²-of-a-bilinear *is* the natural quadratic norm, and a Yukawa coupling
`ψ̄_L Y ψ_R H` *is* a bilinear — but the dimensionful identification "EW scale = `‖Y‖²_F`" has
**no v58/v59 dynamical home** (`03_higgs_bridge_result.md`): v59's actual Higgs (`XiVacuum`, a
Mexican hat on `ℍ≅ℝ⁴`, vacuum `|ξ|²=1/2`) is 4-dimensional and dimensionless; it does **not**
host the 28-dim bilinear.

### R2 — why each component magnitude `= a_ℓ` `[conj]`

`so(8)` democracy forces "all 784 components equal" but **not their value**. Tying the common
component scale of the 28×28 vacuum bilinear to the physical 3×3 Brannen scale `a_ℓ = (Σ√m)/3`
(`KernelEigenvalues.lean`) is a separate dimensionful identification, undriven.

### The "25 unaccounted directions" `[open]`

`L = so(8)` is 28-dim; a democratic full-rank `Y` has **28** singular directions (28 "light"
directions), but only **3** generations are light. The other `28 − 3 = 25` need a home — naively
the Goldstones of `so(8)` broken to a residual `H` with `dim(so(8)/H) = 25`, i.e. `dim H = 3`.

`formalize_bridge.py` §6 checks the classical subgroups: **no** proper subalgebra of `so(8)`
has dimension 25 (the maximal proper ones top out at `dim Spin(7) = 21`), and none leaves exactly
3 unbroken directions. So "3 light + 25 eaten" is **not** a clean single-step `so(8) → H`
breaking. This exposes the sharpest tension (see FINDINGS): the **full-rank democratic** picture
that gives the `784` count predicts **28** light directions, contradicting the **rank-3** lepton
spectrum. The count mechanism (Option 5) and the spectrum mechanism (Brannen 3×3 kernel) pull in
opposite directions on the rank of `Y`.

---

## 4. Lepton-specificity (a real structural constraint) `[emp/thm]`

`√v = D·a` is **lepton-only**. Using running quark scales (`formalize_bridge.py` §3; quark
numbers illustrative/RG-dependent):

| sector | `D` | `√v/a` | wants `D` | verdict |
|---|---|---|---|---|
| lepton | 28 | 28.01 | 28 | **match (0.03%)** |
| d-quark | 35 | 24.2 | 35 | fails ~1.4× |
| u-quark | 63 | 3.4 | 63 | fails ~18× |

Mechanism reason (`03_higgs_bridge_result.md`): only the **color-singlet** lepton's mass bilinear
lives *wholly* in `L`; quark bilinears carry color (F=35, L⊕F=63) and are not democratic-on-L.
The lepton 0.03% is the only *sharp* match; the quark failures are the point — the bridge is not
a generic "`√v = D·a` for any sector."

---

## 5. Files in this gap folder

- `README.md` — this background.
- `ALTERNATIVES.md` — candidate mechanisms for R1/R2/25-directions, each with a test + falsifier.
- `formalize_bridge.py` — self-contained (numpy) computation: empirical anchor, Burnside
  robustness, lepton-specificity, the so(8)-invariant Higgs-quartic R1 test, the norm-selection
  test, and the 25-direction count.
- `EwScaleBridge.lean` — `import Mathlib`, axiom-clean; proves the clean identities and `sorry`s
  exactly the open physics inputs (R1's `v=‖Y‖²_F` hypothesis, R2's value-pin). Written, not
  built this run.
- `FINDINGS.md` — verdict, most promising mechanism, what a full derivation needs.
