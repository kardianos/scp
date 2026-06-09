# Quark-Flavour Cluster — Background (Gaps G4, G5, G6)

**Date**: 2026-05-25 · **Cluster**: quark flavour sector of the v59 gap attack.
**Scope**: the Z₂×Z₂ selection rule (G4), the quark Koide ratios (G5), and the
quark Brannen phases & scales (G6). This document states each gap precisely, the
hypotheses already tried, and — the critical honesty point — quantifies how much
the quark Koide ratios move under RG running.

**Status tags**: **[thm]** machine-checked identity; **[emp≈X]** empirical match at
precision X; **[conj]** structural ansatz, not derived; **[free]** genuine input;
**[open]** no derivation exists.

Sources read: `UNIFIED_THEORY.md` §2/§3/§9-B; `furey_construction/FINDINGS_ckm_and_selection.md`;
`synthesis/FINDINGS_brannen_dynamics.md`; `cosserat_experiment/FINDINGS_selection.md`,
`FINDINGS_quarks.md`, `11_quark_sector.py`, `14_selection_rule_deep_think.py`;
`furey_construction/13_fock_mass_forcing_report.md`;
Lean `Predictions.lean`, `7D_Algebra/{Z2Z2Forcing,LeptonGradeForcing}.lean`.

---

## The single source and its bisection (shared context) **[thm]**

The parent algebra is `Cl(7)_even ≅ Cl(6) ≅ ℂ⊗𝕆` (dim 64). Its Clifford-grade
decomposition on ℝ⁷ is

```
  Λ⁰ = 1     Λ² = 21 (= dim Spin(7))     Λ⁴ = 35     Λ⁶ = 7 (= dim Im𝕆 = dim S⁷)
  1 + 21 + 35 + 7 = 64
```

The **L ⊕ F bisection** is the grade-mod-4 eigenspace split of the involution
`μ|_{Λᵏ} = (−1)^{k/2}`:

```
  L = Λ² ⊕ Λ⁶   (μ = −1, dim 28 = dim Spin(8))   "Lie-algebra content", NO G₂-singlet
  F = Λ⁴        (μ = +1, dim 35)                  "G₂-form content", contains the
                                                   coassociative 4-form *φ (the unique
                                                   even-grade G₂ singlet)
```

The bisection itself, the dimensions, and the additive identity `28 + 35 = 63`
are **[thm]** (`QuarkKoide.lean`, and `Predictions.lean` Tier 7/8). The three
sector ambients are graded sums of the same parent:

| Furey N | sector | (Bit-L, Bit-F) | ambient | D_N |
|---|---|---|---|---|
| 0 | lepton (color singlet) | (1, 0) → L | Λ²⊕Λ⁶ | **28** |
| 1 | d-quark (color triplet) | (0, 1) → F | Λ⁴ | **35** |
| 2 | u-quark (color antitriplet) | (1, 1) → L⊕F | Λ²⊕Λ⁴⊕Λ⁶ | **63** |
| 3 | ν? | ? | ? | ? |

---

## G4 — the Z₂×Z₂ selection rule

**Claim.** Each fermion sector N picks a grade-bit pair `(Bit-L, Bit-F)`, fixing its
mass-space ambient dimension D_N. The bisection L⊕F is **[thm]**; lepton=L is argued
**[thm]/[conj]** (see below); the full map `N ↦ (Bit-L, Bit-F)` is **[open]**.

### What IS settled (more than the gap brief credits)

A genuine, machine-checked advance exists in `7D_Algebra/Z2Z2Forcing.lean` and
`LeptonGradeForcing.lean` (post sign-fix to a genuine Cl(7); `decide`-checked on the
8×8 generators):

- **`color_splitting_requires_F`** **[thm]**: every L-grade generator (all 28 of
  Λ²∪Λ⁶) is diagonal-free, so L carries **no** color Cartan; the F-grade supplies a
  **rank-≥2** color-block Cartan (two explicit 4-forms give `[1,1,−1]` and `[−1,1,−1]`).
  ⇒ a color triplet (d/u quark) *cannot* obtain its color structure from L and **must
  use F**. This pins **d-quark = F**.
- **`composite_color_requires_LF`** **[thm-ingredients]**: L is not closed —
  `γ_iγ_j · γ_kγ_l` (disjoint) lands in F (`disjoint_bivectors_mul_isF`, a universal
  grade law). A *composite* (N=2, two raisings in the Fock/Witt construction) is a
  product of L-excitations and so spills into F; combined with the F color Cartan,
  this pushes **u-quark = L⊕F**. The grade arithmetic is **[thm]**; the "N=2 is a
  product of two raisings" identification is interpretive ([conj]).
- **`lepton_L_not_forced_by_availability`** **[thm]** — an honest *negative* result:
  lepton = L is **NOT** forced by mass-channel availability. In the real 8×8 picture,
  F offers leptons a *richer* (diagonal/symmetric) mass channel that L lacks; on the
  lepton singlet pair {0,7}, **every L-generator is antisymmetric** and **every
  F-generator is symmetric**. So the *real* metric mass actually sits in F — opposite
  to the assignment. The forcing reduces to a single missing input: **which grade
  carries the ℂ-imaginary unit J of the Furey ideal**. If `J ∈ Λ² (⊂ L)`, then the
  L antisymmetric block is the genuine Hermitian lepton mass and lepton=L is forced;
  if `J ∈ F`, the assignment flips. Establishing where J lives is a statement about
  the ℂ-structure of the minimal left ideal, **outside** the real 8×8 matrices — and
  is **[open]**.

Note: `ColorSU3.lean` separately pins the lepton complex structure to ±canonical via
an explicit su(3) annihilating the singlet (UNIFIED_THEORY §2 calls lepton=L
"forced"). That forcing is real at the *complex-structure* level; the
`LeptonGradeForcing` result is the more careful statement that the *grade* assignment
of the mass term still needs J's home pinned. These are consistent: both say the open
residue is "locate J".

### Selection-rule hypotheses tried and why they failed

From `14_selection_rule_deep_think.py` (six hypotheses) + the 7-angle attack in
`13_fock_mass_forcing_report.md`:

| # | hypothesis | result | why it failed |
|---|---|---|---|
| H1 | Witt bidegree (number-preserving ops, p=q) | **null** | gives 1+9+9+1 = 20 for all N; not 28/35/63, and N-independent |
| H2 | Hodge / co-grade involution within Cl(7)_even | **incomplete** | Λ⁰↔Λ⁶, Λ²↔Λ⁴ pairing doesn't yield a clean per-N rule |
| H3 | Z₃ triality × Furey N | **underdetermined** | triality already broken by the Spin(7) decomposition; no Lagrangian |
| H4 | topological winding W = 2N selecting grades | **null** | W↔grade map contradicts observed selections |
| H5 | centralizer of the Furey-N projector in End(ℂ⁸) | **close, wrong** | gives 50/34/34/50; 34 vs 35 is ~3% off, not exact |
| H6 | SU(3)_c-content filtering of the Cl(7) grades | **most physical, partially done** | the G₂→SU(3) branching (`15_su3_branching.py`) is the route that EVOLVED into the `Z2Z2Forcing` color-Cartan argument — i.e. H6 is the one that *partly worked* (d-quark=F), but does not by itself fix lepton=L vs F |
| H7 | Yukawa-Lagrangian projectors / WZW windings | **not derived** | no executable Lagrangian; framework only |

**Verdict on G4**: the d-quark→F half is now **[thm]** (color Cartan); the
u-quark→L⊕F half is **[thm]**-ingredients + interpretive glue; the lepton→L half is
**[open]**, reduced to the well-posed sub-problem "where does J live". The full
dynamical map `N ↦ (Bit-L,Bit-F)` is **not yet derived from a Lagrangian**.

---

## G5 — quark Koide ratios

**Claim** **[emp≈0.3%, RG/scheme-dependent — weak]**:
`t²_N = 1 − dimG₂/D_N = 1 − 14/D_N` with the same Brannen closed form
`Q_N = (1 + 2t²_N)/3` gives

```
  Q_d = 11/15 = 0.73333   (from t²_d = 3/5, D_d = 35)
  Q_u = 23/27 = 0.85185   (from t²_u = 7/9, D_u = 63)
```

The closed forms and the cross-sector invariant `(1 − t²_N)·D_N = 14` are **[thm]**
(`QuarkKoide.lean`). The PHYSICS claim — that the *measured* quark Koide ratio equals
these — is the gap, and it is **weak for two independent reasons**:

1. **It rests on G4.** The denominators D_N are exactly the ambient dimensions whose
   assignment G4 leaves undriven (lepton half open). So G5 inherits G4's weakness.
2. **Quark masses run; Koide is NOT RG-invariant.** This is the decisive honesty
   point, quantified in `quark_koide_rg.py` (see below).

### The RG-invariance problem — quantified

`quark_koide_rg.py` computes Q_d, Q_u from documented MS-bar masses at several scales:

| dataset | Q_d | gap to 11/15 | Q_u | gap to 23/27 |
|---|---|---|---|---|
| **PDG native (mixed scales)** | 0.7314 | **+0.26%** | 0.8490 | **+0.34%** |
| MS-bar @ 2 GeV | 0.7831 | −6.36% | 0.8942 | −4.74% |
| MS-bar @ M_Z | 0.7487 | −2.05% | 0.8879 | −4.06% |
| MS-bar @ m_t | 0.7476 | −1.91% | 0.8932 | −4.63% |
| MS-bar @ 1 TeV | 0.7479 | −1.94% | 0.8959 | −4.92% |

**The 0.3% match appears ONLY at the PDG mixed-scale convention** — where each quark
is quoted at a *different* reference scale (u,d,s at 2 GeV; c at m_c; b at m_b; t at
its pole). That is **not a physical scale**; it is a bookkeeping convention. At any
single, physically meaningful renormalization scale the match degrades to **2–5%**.

- **RG spread of Q_d** across scales = 0.052, which is **27×** the native 0.3% gap.
- **RG spread of Q_u** = 0.047, which is **16×** the native gap.
- Restated on the "clean" invariant: `(1−t²_d)·35` ranges over **[11.4, 14.1]** and
  `(1−t²_u)·63` over **[9.8, 14.3]** — the target 14 is hit only at the mixed
  convention.

**Overfitting control**: within the RG band, **12 distinct simple rationals** p/q
(q≤27) fall inside for *each* of Q_d and Q_u. 11/15 and 23/27 are among them, but so
are 7/9, 3/4, 13/17, 6/7, 8/9, … — being "a simple rational near the band" is not
special.

**Sharper RG diagnosis** (`rg_invariance_test.py`): Koide Q is **exactly invariant
under a common (flavor-universal) rescaling** `m_i → η·m_i`, and one-loop QCD running
*is* flavor-universal (`γ_m` is the same for all flavors). So pure LO-QCD running does
**not** move Q at all. The observed drift therefore comes from two sources: (i)
mismatched *input* scales/threshold conventions between datasets, and (ii) the
**flavor-non-universal top Yukawa**, which gives m_t a flavor-specific anomalous
dimension. Since Q_u is utterly dominated by √m_t, Q_u is **intrinsically unprotected**
(varying m_t from 150→384 GeV swings Q_u's gap from −3.3% to −7.6%). The 0.34% Q_u
match is thus a coincidence of one particular m_t convention (the pole-ish 172.69), not
an RG-invariant statement. Q_d, with no large Yukawa, is better-protected against
running but is dominated instead by the **s-quark mass uncertainty** (PDG ±9%).

**Verdict on G5**: a **pattern, not a prediction**. The structural closed form is a
theorem; its agreement with data is real only at an unphysical scale convention and
inside an RG band wide enough to host a dozen competing rationals. Honest precision
is ~2–5% (the common-scale gap), not 0.3%.

---

## G6 — quark Brannen phases & scales

**Facts** (from `FINDINGS_brannen_dynamics.md`, re-tested in `quark_phases_scales.py`):

- Fundamental-domain phases: **φ_d ≈ +0.1086 rad**, **φ_u ≈ −0.0725 rad**
  (the literature also quotes φ_u ≈ −2.02 rad; same orbit under the generation S₃).
- **φ_q ≠ Q_q/3** **[thm-negative]**: the *exact* lepton relation `φ_l = Q_l/3 = 2/9`
  does **not** extend. Q_d/3 = 11/45 ≈ 0.244 vs φ_d ≈ 0.109 (off by >2×); Q_u/3 has
  the wrong sign entirely.
- Prior numerology: φ_d ≈ 14·α(M_Z) [dim G₂·α], φ_u ≈ −10·α(0) [dim Spin(5)·α] at
  ~0.6–0.8%; and 3φ_d ≈ 1/9·3 = 1/3, 3φ_u ≈ −3/14.
- Scale ratios: a_u²/a_d² ≈ 35, a_u²/a_l² ≈ 72.

### Status (with overfitting caveats)

- **The phase is observable only through cos(3φ)** (φ is fixed only mod S₃). Because
  the quark phases are small, cos(3φ) ≈ 0.95–0.98, so *every* candidate matches the
  observable cos(3φ) to <0.7% — cos(3φ) is a **weak discriminator** here.
- The `n·α` forms: only one integer n lands within 1% of each φ (grid spacing of n·α
  near the phase is 7–10% of φ). So `14·α` / `10·α` are *not* obviously overfit on φ
  itself — but with the wrong-α ambiguity (α(0) vs α(M_Z)) and the choice of two
  integers, the "structural integer × α" reading is **suggestive at best**, not a law.
- **Scale ratios are NOT scale-stable** **[emp/artifact]**: a_u²/a_d² = 35 and
  a_u²/a_l² = 72 hold only at the PDG mixed-scale convention (same artifact as G5). At
  M_Z, a_u²/a_d² ≈ 49; at 2 GeV, a_u²/a_l² ≈ 153. (a_l² is scale-stable since QED
  running of charged leptons is negligible; the quark a's move with m_t, m_b.)

**Verdict on G6**: φ_d, φ_u and a_u, a_d are **[free]**. No structural law derives
them; the lepton φ=Q/3 relation explicitly fails; the α-numerology survives only the
weak cos(3φ) observable; the scale ratios are convention artifacts.

---

## Cross-gap dependency

```
   G4 (assignment N→D_N)  ──undriven (lepton half open)──►  feeds D_N into
        │
        ▼
   G5 (Q_N = (1+2(1−14/D_N))/3)  ──RG/scheme-dependent──►  pattern, not prediction
        │
        ▼
   G6 (φ_q, a_q)  ──free, no law──►  numerology only on weak observables
```

The whole quark block is the weakest tier of v59: G4's lepton half is open, G5 is a
soft RG-dependent pattern, G6 is free. See `ALTERNATIVES.md` for proposed avenues and
`FINDINGS.md` for the consolidated verdict.
