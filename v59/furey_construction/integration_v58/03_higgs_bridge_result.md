# Integration #1, the Cl(7) bridge: RESULT (computed)

*2026-05-24.  Per the steer "attempt the Cl(7) bridge first" — check, in the **actual**
Furey Cl(7)_even 8×8 representation, whether the lepton mass scale `a_l` relates to the
28-dim `L`-grade by the structural factor `28 = dim(L)` that `v_Higgs = dim(L)²·a_l²`
requires.  Computation: `bridge_vhiggs_cl7.py` (self-contained, numpy).  **Outcome: the
mechanism is now pinned — and the previously-stated "equipartition" mechanism is wrong.***

## What was computed

Built Cl(7) as seven anticommuting 8×8 matrices (`γ_a² = −I`), formed the 28 `L`-grade
blades (`Λ²`=21 ⊕ `Λ⁶`=7) and the 35 `F`-grade blades (`Λ⁴`), and verified:

- **grade ↔ square sign** (matches `LeptonComplexStructure.lean`): all **28/28** `L`-blades
  square to `−I` (complex structures), all **35/35** `F`-blades to `+I`.
- **the 28 `L`-blades are exactly Frobenius-orthogonal**: their Gram matrix
  `⟨A,B⟩ = Re Tr(A†B)` is `8·I₂₈` — diagonal, zero off-diagonal to machine precision.

That orthogonality is the decisive fact, because it fixes what "spread `√v` over the 28
`L`-directions" actually gives under the physical (energy / `L²`) norm.

## The result: 28 comes from the *bilinear*, not from equipartition

Two readings of "the EW vacuum lives on the 28-dim `L` with per-element scale `a_l`":

| reading | object | natural norm | `√v / a_l` | verdict |
|---|---|---|---|---|
| **(A) vacuum = vector in `L`** | `V = a_l Σ_{i=1}^{28} ê_i` (democratic) | `L²` (Frobenius / energy) | **`√28` = 5.29** | **WRONG** (off by `√28`) |
| (A′) same, but added **coherently** | `Σ` of amplitudes (`L¹`) | `L¹` (linear sum) | `28` | needs `L¹` = √-mass linearity (= Koide) |
| **(B) vacuum = bilinear on `L`** | Yukawa `Y` = a `28×28` matrix, each entry `~a_l` | `L²` (Frobenius², over the `dim(L)²` **entries**) | **`28`** ✓ | **CORRECT, natural norm** |

- **(A) is what `INTEGRATION_PLAN_v58.md §(1c)` asserted** — "equipartition of `√v` over those 28
  directions gives `a_l = √v/28`."  The rep shows that is **false**: 28 *orthogonal* directions
  each carrying `a_l` have total energy-norm `√28·a_l`, so equipartition predicts
  `√v = √28·a_l ≈ 5.3 a_l`, not `28 a_l`.  **The equipartition mechanism is excluded.**
- **(B) is correct and natural.**  A lepton mass / Yukawa coupling *is* a bilinear — a matrix
  `Y` in `ψ̄_L Y ψ_R H`.  On the 28-dim mass-bearing grade `L` it has `dim(L)² = 784`
  components.  With the **proven** `so(8)=L` symmetry making the 28 directions equivalent, a
  democratic `Y` has every component of common magnitude `a_l`, so

  > `v_Higgs = ‖Y‖²_Frobenius = (number of components) · a_l² = dim(L)²·a_l² = 28²·a_l²`,

  using the **ordinary `L²` norm**.  No `L¹` / coherence assumption.  **`28²` is literally the
  number of entries of the `L`-grade mass bilinear**, and the empirical `√v/a_l = 28.0096`
  (0.034%) is reproduced.

So the suspicious `784` is neither a fitted coincidence nor the (incorrect) vector
equipartition: it is **`dim(L)²` = #components of the so(8)-democratic lepton mass bilinear**,
under the natural Frobenius norm.

## Lepton-specificity (computed, confirms the picture)

`√v = D·a` is **lepton-only**.  Using running quark scales `a_q = (Σ√m_q)/3`:

| sector | `D` | `√v / a` (actual) | matches `D`? |
|---|---|---|---|
| lepton | 28 | **28.01** | **yes (0.03%)** |
| d-quark | 35 | 19.5 | no (1.8×) |
| u-quark | 63 | 3.3 | no (19×) |

Consistent with the construction: only the **color-singlet** lepton's mass bilinear lives
*entirely* in `L` (28-dim).  The d-quark sits in `F=Λ⁴` (35), the u-quark in `L⊕F` (63), both
with color structure — so their bilinears are not democratic-on-`L`, and `√v=D·a` correctly
fails.  (Same lepton-specificity flagged in `INTEGRATION_PLAN_v58.md`, now with the *reason*:
the Frobenius-bilinear count only applies to the singlet whose `Y` is wholly `L`-grade.)

## What this closes, and the irreducible residual

**Closed / sharpened:**
- The mechanism is the **Frobenius² of the so(8)-democratic lepton mass bilinear on `L`**,
  `v = dim(L)²·a_l²` — natural `L²` norm, `so(8)` democracy (proven `L=so(8)`) doing real work.
- The **equipartition reading is excluded** (gives `√28`, not 28) — a real correction to the
  prior note.
- **Lepton-specificity is derived**, not assumed: the count applies only to the singlet whose
  `Y` is entirely `L`-grade.

**Irreducible residual (the honest gap — neither v58 nor v59 currently supplies it):**
1. **(R1) why `v_Higgs = ‖Y‖²_F`** — why the EW scale is the *Frobenius²* of the lepton mass
   bilinear (not its trace, top eigenvalue, …), and why that equals 246 GeV.  This is the
   dimensionful "vacuum field = mass bilinear" identification; v58 has no lepton sector and
   v59's actual Higgs (`XiVacuum`, the Mexican hat on `ℍ≅ℝ⁴`, vacuum `|ξ|²=1/2`) is
   4-dimensional and dimensionless — **neither hosts the 28-dim bilinear.**
2. **(R2) why each component magnitude = the physical Brannen `a_l`** — `so(8)` forces "all 784
   equal" but not the *value*.  Tying the `28×28` vacuum bilinear's scale to the `3×3` physical
   Brannen kernel (`a_l = (Σ√m)/3`, `KernelEigenvalues.lean`) is a separate identification.
3. **Open consistency:** `L` is 28-dim but only **3** generations are light.  The democratic
   bilinear has 28 singular directions; the other 25 are unaccounted — plausibly the `so(8)`
   gauge/Goldstone content (cf. `XiVacuum`'s 3 Goldstones = `W±,Z`), but this is not worked out.

## Net

The Cl(7) bridge **does not "fully derive" `v_Higgs`** — (R1)+(R2) are genuine dimensionful
identifications with no dynamical home yet.  But it **converts the relation into a precise,
defensible structural statement** and **kills one wrong mechanism**:

> `v_Higgs = ‖Y_lepton‖²_Frobenius` on the 28-dim `L`-grade `= dim(L)²·a_l²`, with `so(8)`
> democracy fixing the component count `dim(L)²` and the equal magnitudes — **not** the
> (excluded) equipartition `√28`.

This is "as derived as Koide + the bilinear identification": the `28` is the same `dim(L)=so(8)`
that sets the lepton Koide ratio (`t²=1−14/28=1/2`), now appearing as the **component count** of
the mass bilinear.  Remaining is the dimensionful step `v = ‖Y‖²_F` with scale `a_l` — recommend
banking #1 here and moving to `α`.
