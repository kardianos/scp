# v59 Version Closeout — Positive and Negative Findings (detailed)

**Date**: 2026-05-25 · **Purpose**: the version-closing ledger for v59. Unlike
`v59/UNIFIED_THEORY.md` (the cohesive theory presentation), this document is an *audit*:
every positive and negative finding, with the critical equation and derivation step inline,
each traced to the file that backs it (path relative to the project root `scp/`), and a
verification status. It closes with the primary goal, the gaps to it, and suggested avenues.

**Status tags**: **[thm✓]** machine-checked Lean (built & verified this session) · **[emp≈X]**
empirical match at precision X · **[num✓]** numerically verified by a script that runs ·
**[conj]** conjecture, not derived · **[null]** tested and found false/absent · **[free]** input.

**Path convention**: paths are relative to the project root `scp/`. Bare core Lean module names
(e.g. `BladeSquareSign.lean`) live in `v59/furey_construction/lean/`; names prefixed `7D_Algebra/`
are in `v59/furey_construction/lean/7D_Algebra/`. Gap modules carry their full `v59/gaps/...` path.

**Verification basis** (see `v59/gaps/SYNTHESIS.md` §"Verification status"):
- Core Lean (`v59/furey_construction/lean/`, Lean 4.29.0 + Mathlib): `lake build` succeeds
  (8278 jobs); `AxiomCheck` shows every headline theorem uses only
  `{propext, Classical.choice, Quot.sound}`.
- Six gap modules `v59/gaps/*/*.lean`: all typecheck, 0 errors; `sorry` only on the flagged
  conjectures (α, R2, the α²¹ match, G4/G5/G6).
- All numerical scripts under `v59/synthesis/` and `v59/gaps/` run clean.

---

# PART I — POSITIVE FINDINGS (established)

## P1. Single source algebra and its grading **[thm✓]**

The fundamental object is `Φ ∈ 𝒜 = Cl(7)_even ≅ Cl(6) ≅ ℂ⊗ℍ⊗𝕆`, `dim_ℝ 𝒜 = 64`, graded
```
Λ⁰ ⊕ Λ² ⊕ Λ⁴ ⊕ Λ⁶  =  1 + 21 + 35 + 7  = 64 ,   C(7,2)=21=dim Spin(7),  C(7,6)=7=dim Im𝕆.
```
Bisection `L = Λ²⊕Λ⁶ (28) , F = Λ⁴ (35)`, with `28+35=63` and `28 = 2·dim G₂ = dim Spin(8)`.
*Files*: `v59/furey_construction/lean/SpinDimension.lean`, `Predictions.lean`
(`dim_cl7_even`, `lepton_ambient_decomp`, `v59_D_sum_identity`).

## P2. lepton = L is forced (not assumed) **[thm✓]**

A simple even k-blade squares to `(−1)^{k(k+1)/2}·𝟙`, so
```
B∈Λ²,Λ⁶ ⇒ B² = −𝟙  (complex structures)  ⇒ L = skew = so(8)
B∈Λ⁴    ⇒ B² = +𝟙  (real/involution)     ⇒ {Λ⁰}⊕F = symmetric.
```
A genuine complex structure `J` (carrying the Brannen phase, `J²=−𝟙`) must therefore lie in
`L`. The color `su(3)` and color complex structure `J_c=γ₀γ₅∈Λ²` split `ℝ⁸=ℂ⁴=𝟏⊕𝟑` (lepton
singlet ⊕ 3 colors) and pin the lepton `J` to ±canonical.
*Files*: `BladeSquareSign.lean` (`prod_sq`, `even_grade_complex_structure_dichotomy`),
`7D_Algebra/LeptonComplexStructure.lean`, `LeptonRealityForcing.lean`, `ColorSU3.lean`
(`colorSU3_closes`, `colorInvariant_pins_lepton`).

## P3. Koide `Q = 2/3` from the Brannen kernel **[thm✓ for the algebra; emp≈10⁻⁵ for the match]**

Mass kernel on the 3-generation space: `M = a(𝟙 + ξS + ξ̄S²)`, `S³=𝟙`, `|ξ|²=t²`. Eigenvalues
`λ_k = a(1+2t cos θ_k)`, `θ_k=φ+2πk/3`; identify `λ_k=√m_k`. Using `Σcos θ_k=0`,
`Σcos²θ_k=3/2`:
```
Σ√m = 3a ,  Σm = 3a²(1+2t²)  ⇒  Q ≡ Σm/(Σ√m)² = (1+2t²)/3 .
t² = 1/2  ⇒  Q = 2/3 ,  matched to  dim G₂/dim Spin(7) = 14/21 = 2/3 .
```
*Files*: `BrannenKernel.lean` (`Q_value`, `koide_iff_constraint`), `SpinDimension.lean`
(`koide_ratio_structural`), `CyclicShift.lean` (`sum_one_omega_omega_sq`). Masses reproduce to
3×10⁻¹⁴ given the one scale `a`.

## P4. Brannen phase `φ = Q/3 = 2/9`, and the "/3" is structural **[emp≈10⁻⁵; /3 thm-verified]**

`φ = Q/3 = 2/9` (|2/9 − φ_meas| ≈ 7×10⁻⁶, 0.9σ — Koide-class). The divisor `/3` is the
sedenion `S₃` automorphism `ψ` (order 3): `Aut(𝕊)=G₂×S₃`, `ψ³=𝟙`, an automorphism on all
16×16 basis pairs.
*Files*: `SpinDimension.lean` (`brannen_phase_structural`), `BrannenPhase.lean`
(`Q_phase_independent`), and the sedenion construction in `koide_phase_law/sedenion_s3.py`.

## P5. Gauge group, `sin²θ_W`, and the dual Coxeter 5 **[thm✓]**

```
Spin(7) = G₂ × SU(2)_L × SU(2)_R × U(1)_{B−L} ,   14+3+3+1 = 21 .
sin²θ_W = c_{B−L}/(c_W + 2c_{B−L}) = 2/(5+4) = 2/9 ,   cos²θ_W = 7/9 = t²_u .
5 = h∨(Spin 7) = dim Spin(7) − dim Cl(3,1) = 21 − 16 .
```
*Files*: `ScaleBridge.lean` (`sin_sq_thW_eq_brannen_phase`, `cos_sq_thW_eq_t_sq_u_quark`,
`spin7_pati_salam_decomp`), `GaugePrefactorDualCoxeter.lean`.

## P6. The 784-dim operator space is *forced* (the EW-bridge upgrade) **[thm✓ + num✓]**

The bridge `v_Higgs = dim(L)²·a_ℓ² = 784 a_ℓ²` lives on `End(L)`. Non-associativity forces
`|Φ|²` onto the associative left-action algebra; for `L=so(8)` the product *is* the bracket, so
the operator algebra is `𝒜_L = ⟨ad_X : X∈so(8)⟩`. `ad(so(8))` is the absolutely irreducible
28-dim adjoint, so by Burnside/Jacobson density
```
𝒜_L = End(L) = M₂₈(ℝ) ,   dim 𝒜_L = 28² = 784 .
```
Numerically confirmed and shown *generic* (defining rep → 64; adjoint → 784; every abs.-irred.
`so(n)` adjoint → `(dim)²`; reducible `so(4)` → 18 < 36, so the test detects reducibility).
*Files*: `v59/synthesis/obe_options_2_5.py`, `v59/gaps/ew_scale_bridge/formalize_bridge.py`
(num✓); `v59/gaps/ew_scale_bridge/EwScaleBridge.lean` (`endL_dim : finrank ℝ (Matrix (Fin 28)
(Fin 28) ℝ) = 784`, thm✓ via `Module.finrank_matrix`); `HiggsVevReframe.lean`
(`vHiggs_sqrt_form`, `ratio_form`: `√v=28a ⟺ Σ√m/√v = 3/28`, thm✓).

## P7. Gravity charge = second moment (equivalence-principle-exact) **[num✓]**

The live gravitational charge is the second moment of the mass kernel:
```
ρ_grav = Tr(M†M) = Σ_k m_k = 9 Q a²  (= 6a² for leptons) ,   charge/mass = 1 (all sectors) .
Σ m_lepton = (9Q/dim(L)²)·v = (6/784)·v  ≈  0.07% .
```
This is the same Frobenius² structure as the EW bridge (P6), over the 3-dim generation space.
*File*: `v59/synthesis/gravity_charge_test.py` (charge/mass = [1,1,1]; the constraint-deviation
candidate is dead — see N7).

## P8. Gravity is reachable as `1/r²` **[num✓]**

For a massless kernel `K` (`−∇²K=δ`), the OBE gravity term integrates by parts:
```
Ω_grav = f_g ∫ K ∇'ρ_grav  =  f_g ∇(K*ρ_grav) = f_g ∇Φ ,   −∇²Φ = ρ_grav .
```
Since `∫ρ_grav d³x = total mass ≠ 0` (nonzero monopole), `Φ ~ 1/r` (slope −1.000) and the force
`~ 1/r²` (slope −2.000). Recommended OBE form: `□ Ω_grav = f_g ρ_grav`.
*File*: `v59/synthesis/obe_radial_test.py` (massless+monopole → slope −2.000; zero-monopole or
massive → screened).

## P9. Gravity magnitude: the exponent is `dim Spin(7)` **[conj; exponent num✓ to 0.002%]**

```
α_G(electron) = (m_e/M_Pl)² = 1.752×10⁻⁴⁵ .   G_e = (21/16) α²¹ ,  α = α(0) (IR).
exponent fitting prefactor 21/16:  n = 21.0005 = dim Spin(7) .
hierarchy:  m_e/M_Pl = √(21/16) · α^(21/2) .
```
`f_g ~ α^(21/2)` — the historical V6 ~10⁴⁰ overshoot is *exactly* the missing `α²¹`. The
exponent is data-forced to 21 (neighbours miss by 10^±15); `α^(21/2)=det(√α·I₂₁)` = a Λ²¹
top-form = instanton `e^{−S}`, `S=21·ln(1/α)`.
*Files*: `v59/synthesis/gravity_magnitude_test.py`, `v59/gaps/gravity/g8_exponent_test.py`
(num✓); `v59/gaps/gravity/G8G9_Gravity.lean` (`gauge_index_additive : 5=21−16`,
`gravity_prefactor : 21/16`, thm✓; the empirical α²¹ match is the module's single flagged `sorry`).

---

# PART II — NEGATIVE FINDINGS (null results, demotions, dead ends)

*These are first-class results — the project's value is as much in what it killed as in what it kept.*

## N1. `g_W² = 5√α` is NOT an independent law **[null]** — it collapses into α

Given the SM-definitional `4πα = g_W² sin²θ_W` with `sin²θ_W = 2/9`:
```
4πα = (5√α)(2/9)  ⇒  √α = 5/(18π)  ⇒  α(M_Z) = 25/(324π²) .
```
So `g_W²=5√α ⟺ α(M_Z)=25/(324π²)` — a line meeting a √ at one point; it carries no information
beyond the single value of `α`. The "5"=h∨ survives [thm✓]; the `√α` form does not.
**Consequence**: G2 folds into G3; the rigorous input set is exactly `a_ℓ + α`.
*Files*: `v59/gaps/alpha_couplings/verify_forms.py`, `AlphaCouplingIdentities.lean`;
`v59/furey_construction/RIGOR_AUDIT.md`.

## N2. Quark Koide (11/15, 23/27) is a scale-convention artifact **[null as a prediction]**

`t²_N = 1 − 14/D_N` gives `Q_d=11/15, Q_u=23/27`, but Koide `Q` is not RG-invariant: at any
single physical scale the gap is **2–5%** (not 0.3%); the RG spread of `Q_d,Q_u` is **16–27× the
claimed gap**, and 12 simple rationals fall inside each RG band. The "0.3%" holds only at PDG's
mixed-scale convention. **Verdict: a soft pattern, not a prediction.**
*Files*: `v59/gaps/quark_flavour/quark_koide_rg.py`, `rg_invariance_test.py`, `QuarkKoide.lean`.

## N3. Gravity is a SCALAR — no `h=±2` tensor mode **[null — and fatal for LIGO]**

Rigorous helicity decomposition: the carrier `ρ_grav = Tr(M†M)` is a Lorentz scalar ⇒ `h=0`
only. The internal `so(8)/Λ²` index is **inert under spacetime rotations** (helicity is a
spacetime little-group label), so "gravity lives in the 28-dim bivector L" supplies **no** `h=±2`.
Only a symmetric *spacetime* rank-2 `h_μν` carries `h=±2`; the bivector can only *source* one
quadratically (`T_μν=Ω_μ^a Ω_ν^a`), which is `O(α)`/short-range, not a long-range graviton.
**This outranks the magnitude (P9): a perfectly-tuned `α²¹` scalar force still fails LIGO.**
*File*: `v59/gaps/gravity/g9_polarization_test.py`, `g9_spin2_route_test.py`.

## N4. The rank tension — the EW-bridge `Y` and the spectrum `Y` cannot be the same **[null/open]**

The 784 count (P6) needs a **full-rank** democratic `Y` (all 784 entries `~a`, `‖Y‖²_F=784a²`).
But the physical spectrum is **3 generations** = a **rank-3** `Y`, whose
```
‖Y_rank3‖²_F = Σm = 9Qa²   (= the gravity charge, P7)  ≠  784 a² .
```
The same mass bilinear is read two incompatible ways. The "25 unaccounted directions" is the
symptom, and **no single-step `so(8)→H` breaking gives 3 light directions** (no 25-dim subalgebra
of `so(8)`; maximal proper subalgebras top out at 21). This links the EW and gravity sectors and
is the sharpest open problem of the lepton/EW block.
*File*: `v59/gaps/ew_scale_bridge/formalize_bridge.py`; `v59/gaps/ew_scale_bridge/EwScaleBridge.lean`
(R1/R2 stated, the value-pin `sorry`-flagged).

## N5. `cos(2/3)` is not an algebra invariant; the phase is not geometric **[null]**

The phase invariant is `cos 3φ = cos(2/3)` — *transcendental*. No Lie/Casimir/√3 quantity equals
it (best near-miss 11/14, off by 1.7×10⁻⁴ ≫ the 10⁻⁵ data precision). Independently, `2/9` and
`2/3 rad` are **not π-rational**, so the phase is **not** a geometric/holonomy/rotation angle —
the geometric/Berry-phase mechanism was built and **ruled out** (v58-connection, color-Z₃, and
sedenion-ψ holonomies all fail to force the value).
*Files*: `v59/gaps/lepton_phase/search.py`, `second_invariant.py`; `PhaseExclusions.lean`,
`PhaseAmbiguity.lean` (`phase_not_pi_rational`, `koide_not_pi_rational`, thm✓).
*Positive remainder (N5→P)*: the **magnitude** `Q=2/3` does have non-geometric readings —
`t²=(D−dim G₂)/D` (max-mixing) and `t²=1/2` (equipartition) — both Lean-proved in
`LeptonPhaseMagnitude.lean` (`massSplit_lepton_half`, fully proved, 0 `sorry`); and `cos 3φ` is
exactly the standardized **skewness** of `√m` (`cos3phi_is_third_moment`, thm✓). What stays open
is the *radian-insert*: why `Q` re-enters as the phase argument.

## N6. The selection rule and quark phases are undriven **[open/free]**

`d-quark→F` is machine-checked, `u→L⊕F` thm-up-to-glue, but **lepton→L by availability is not
forced** (reduced to one decidable question: is the ℂ-unit `J∈Λ²`?). Quark phases `φ_d≈0.110,
φ_u≈−2.02` are **free** (`φ_q ≠ Q_q/3` — the lepton relation does not extend); the `n·α`
numerology only tests on the weak `cos3φ`; scale ratios (35, 72) are convention artifacts.
*Files*: `7D_Algebra/LeptonGradeForcing.lean`, `Z2Z2Forcing.lean`;
`v59/gaps/quark_flavour/quark_phases_scales.py`.

## N7. The constraint-deviation gravity charge is dead **[null]**

The earlier candidate `ρ_N = |Π_ℍΦ|² − (1−14/D) ~ a²(|ξ|²−½)` is **zero for leptons**
(`|ξ|²=½`) and scale-free (drops `a`), so it cannot be mass: charge/mass is non-universal
(2.4× spread) and the u/d ratio is 2.8 (factor) or 98 (×a²), not the empirical 40.6.
*File*: `v59/synthesis/gravity_charge_test.py` (superseded by P7).

## N8. `g² = h∨√α` does not generalize to the strong sector **[null]**

The EW ansatz fails for QCD: `g_s² = 3√α ≈ 0.26` vs measured `≈1.48` (off >80%). A `√α` form
cannot make an `O(1)` coupling. `α_s` is free.
*File*: `v59/gaps/frontier_sectors/pmns_strong_test.py`, `FrontierSectors.lean`
(machine-checkable falsifier `3√α < 1 < 1.48`, thm✓).

## N9. CKM/PMNS are mostly free; only `sin²θ_C` survives, coincidence-grade **[mostly null]**

`sin²θ_C = 7α(0)` (7=dim Im𝕆) → `sinθ_C = 0.226` vs 0.225 (0.45%), but coincidence-grade (needs
the IR α(0); no mechanism; ~10⁴× looser than the lepton Koide). `δ_CP` is hit by *two* unrelated
forms (the overfitting signature). `V_ub, V_td, δ_CP`, neutrino masses, and the PMNS matrix are
**free** (the neutrino Koide ambient `D_ν` fits no v59 ambient).
*Files*: `v59/gaps/frontier_sectors/ckm_overfitting_test.py`, `pmns_strong_test.py`.

## N10. The α value itself is not derived **[conj]**

Both `−lnα+2α=π²/2=8π²/16` (IR) and `α(M_Z)=25/(324π²)` (EW) are conjectures. Within a
form-family the v59 forms are uniquely tight, but the look-elsewhere across families is large (5
of 8 arbitrary templates already hit α to 0.03%); precision alone is not evidence.
*Files*: `v59/gaps/alpha_couplings/overfit_scan.py`, `rg_fixedpoint_test.py`, `AlphaZero.lean`.

---

# PART III — PARAMETER BUDGET (what v59 actually reduced)

| block | optimistic (accept conjectures) | rigorous (theorems + tight emp) |
|---|---|---|
| charged-lepton + EW + Higgs | **1 input: `a_ℓ`** | **2 inputs: `a_ℓ` and `α`** |
| full Standard Model | **~10–12 free** — quark scales ×2, quark phases ×2, CKM ×4, neutrinos ≥7, `θ_QCD`, `α_s` (all untouched) | same |

Turned dependent (machine-checked or tight-emp): the 2 non-scale lepton d.o.f. (→`Q`),
`v_Higgs`, `sin²θ_W`, `cos²θ_W`, `m_W`, `m_Z`, `m_H`, the gauge integers `(5,2)`. **Not**
reduced: the quark-flavour / CKM / neutrino / strong sectors. *Source*:
`v59/furey_construction/PARAMETER_BUDGET.md`, `RIGOR_AUDIT.md`.

---

# PART IV — PRIMARY GOAL, GAPS, AND SUGGESTED AVENUES

## Primary goal

> **A single algebraic field theory — the OBE on `Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆` — whose structure
> reproduces the Standard Model spectrum, gauge group, electroweak scale, and gravity from one
> dimensionful scale `a_ℓ` plus pure mathematics (no fitted dimensionless constants).**

Where v59 stands against it: the *algebraic skeleton* is theorem-grade (P1–P6); the lepton+EW
block reduces to `a_ℓ + α`; gravity has a live charge and radial law (P7–P8) and a striking
magnitude exponent (P9). But every *physical* reduction is an empirical match hanging on
conjectures, and three structural obstructions (N3, N4, N10) sit between the current state and the
goal.

## Gaps to the goal (ranked)

1. **G9 — no tensor gravity mode (N3).** Decisive: the OBE's gravity is scalar and fails LIGO.
   Gates the entire gravity sector; the magnitude (P9) is moot until this is solved.
2. **G1 — the rank tension (N4).** The EW-bridge `Y` (full-rank, 784) and the spectrum/gravity
   `Y` (rank-3, `Σm`) cannot be the same object; this is the real content of the bridge residuals
   R1/R2 and links the EW and gravity sectors.
3. **G3 — `α` is underived (N10);** subsumes the dead `g_W²=5√α` law (N1). The only dimensionless
   input.
4. **G7 — the radian-insert (N5).** Why `Q` re-enters as the phase argument (the magnitude is now
   live, the insert is not).
5. **Selection rule / quark flavour (N2, N6)** and **CKM/ν/strong (N8, N9)** — soft or free;
   patterns, not predictions.
6. **G13 — the dynamical Lagrangian.** Deliberately deferred: its field content depends on
   resolving G9 (does it need an `h_μν`?) and G1 (is `Y` full-rank or rank-3?).

## Suggested avenues of investigation

- **G9 (highest priority): an induced / emergent-metric recast.** Treat `Ω∈Λ²` as a fundamental
  2-form `B` (Plebański-style) with `h_μν` *derived* via a soldering/simplicity constraint to the
  spacetime Lorentz group; make-or-break tests = the soldering constraint exists, and a clean
  count of 2 transverse-traceless DOF. *(Start in `v59/gaps/gravity/`.)*
- **G1: resolve the rank tension via a two-piece `Y`.** Read R1 as one EW-and-gravity-shared
  conjecture ("scale = Frobenius² of the mass bilinear", consistent at 0.07% via
  `Σm_ℓ=(6/784)v`); identify the Lean-proven `XiVacuum` Higgs with the *active rank-3 block* of
  `Y∈End(L)`, and account for the other 25 directions as an `so(8)` Goldstone/gauge sector (a
  chain, not single-step breaking).
- **G3: an RG / β-function fixed point.** Run the theorem-grade `(c_W,c_R,c_{B−L})=(5,5,2)·√α`
  pattern and test whether a single scale jointly fixes `sin²θ_W` and `α(M_Z)` — this would
  *explain away* the `5√α` coincidence (N1) as the value standard running passes through at `M_Z`.
- **G7: a character/spectral origin for `cos(2/3)`.** Feed the live magnitude readings (max-mixing
  + equipartition) into the skewness reframing: seek a v59 character or weighted eigenvalue
  (e.g. `χ(ψ)` of the order-3 sedenion automorphism, weighted by the `G₂`-content Casimir) that
  yields `cos(2/3)` with `Q` entering as a *weight*, not a π-rational angle.
- **Strong CP (a clean "explain a zero"):** test whether `θ_QCD ≈ 0` is *forced* by the color
  complex structure `J_c` (theorem-grade) via the same reality/skewness argument that pins the
  lepton `J` — explaining an observed zero rather than fitting a number.
- **Only after G9 + G1 are resolved:** write the dynamical Lagrangian `ℒ_v59` whose
  Euler–Lagrange equations yield the OBE `Ω(x)` structure (the v60 deliverable).
