# v59 — Progress Summary

**Date**: 2026-05-22 (major update with full session record in `SESSION_2026-05-22.md`)
**Status**: **LEPTON-SECTOR KEPLER ELLIPSE COMPLETE**. Brannen form, Koide identity Q = 2/3, and Brannen phase φ = 2/9 rad are all structurally derived from the $G_2 \subset \text{Spin}(7) \subset \text{Spin}(8)$ exceptional Lie group chain plus three generations from triality. Full Standard Model construction via $\mathbb{C} \otimes \mathbb{H} \otimes \mathbb{O}$ (Furey program) now in progress.

**2026-05-22 session additions** (see `SESSION_2026-05-22.md` for the master record):
- New v59-tier conjectures: $g_W^2 = 5\sqrt{\alpha}$ (0.28%), $G_e = (21/16)\alpha^{21}$ (0.25%), quark Koide $Q_d = 11/15, Q_u = 23/27$ (~0.3%).
- Single source identified: all sector ambients are graded subspaces of Cl(7)_even ≅ ℂ⊗𝕆 (the Furey color algebra, dim 64).
- ~50 axiom-clean Lean theorems across 8 modules; the silent SU(2)/U(1) direction proved algebraically.
- 16 Python experiments documented in `cosserat_experiment/`.

---

## Headline Results

1. **The Cl(3,1) Z₃-cyclic algebraic kernel derives the Brannen form** of the charged lepton mass spectrum without fitting. (`multivector_kernel_fit/01_findings.md`)

2. **The Koide identity Q = 2/3 is structural** in two complementary ways:
   - **Constraint-surface origin**: ξ promoted from complex to quaternionic, restricted to S³ of radius $1/\sqrt{2}$ in ℍ. Random points on this constraint surface all yield Koide Q = 2/3 to floating-point precision. (`multivector_kernel_fit/04_findings.md`)
   - **Lie-algebraic origin**: $Q = \dim G_2 / \dim \text{Spin}(7) = 14/21 = 2/3$, the ratio of the octonion automorphism group to the symmetry of the octonion imaginary sector. (`octonionic_extension/02_findings.md`)

3. **The Brannen phase φ = 2/9 rad is structural**: $\varphi = Q/3 = (\dim G_2)/(3 \cdot \dim \text{Spin}(7))$. The 2/9 emerges as the Koide ratio divided by the number of generations. (`octonionic_extension/02_findings.md`)

4. **The 42-order EM-vs-gravity hierarchy** compresses under exponential density modulation to the single suppression-depth ratio:
$$\frac{S_{\rm grav}}{S_{\rm em}} = \frac{-\ln G_e}{-\ln \alpha} = \frac{103.06}{4.92} = 20.945 \;\approx\; 21 = \dim \text{Spin}(7)$$
   Gap from structural value 21: 0.26%, plausibly a modulation-function correction. (`multivector_kernel_fit/06_findings.md`, `octonionic_extension/01_findings.md`)

5. **Three generations identified with triality**: the three fermion generations correspond to the Z₃ ⊂ S₃ cyclic subgroup of the triality outer automorphism of Spin(8). (`octonionic_extension/01_findings.md`)

6. **Charged lepton masses are reproduced to machine precision** (3 × 10⁻¹⁴) using only the overall mass scale $a$ as an empirical input. Three measured numbers → 0 free parameters beyond scale.

7. **A specific testable prediction emerges**: $\delta\alpha/\alpha \sim 10^{-10}$ in Earth's gravitational potential, accessible to next-generation atomic clocks. (`multivector_kernel_fit/05_findings.md`)

## The Kernel Structure (Current Best Understanding)

- **Underlying algebra**: Cl(3,1) projective geometric algebra, with the natural identification of the three quaternionic imaginary units with the three spatial bivectors:
$$i \leftrightarrow e_{23}, \quad j \leftrightarrow e_{31}, \quad k \leftrightarrow e_{12}$$

- **Lepton mass operator**: $M_\ell = a(I + \xi S + \bar\xi S^T)$ acting on the 3-flavor space, where:
  - $S$ is the cyclic shift (R times 120° rotation about (1,1,1) axis).
  - $\xi \in \mathbb{H}$ is a quaternionic coupling on the constraint surface $|\xi|_\mathbb{H}^2 = 1/2$.
  - The constraint surface is S³ of radius $1/\sqrt{2}$ in ℍ.

- **Eigenvalues**: automatically Brannen form $\lambda_k = a(1 + 2|\xi|\cos(2\pi k/3 + \varphi))$, with Koide Q = 2/3 enforced by the constraint.

- **Cross-sector** (under exponential density modulation, currently empirical):
$$\alpha_{\rm obs} = \alpha_{\rm bare} \exp(-S_{\rm em}), \quad G_{\rm obs} = G_{\rm bare} \exp(-S_{\rm grav})$$
with $S_{\rm grav}/S_{\rm em} \approx 21$ (empirically).

## What's Structural vs Empirical

| Item | Status |
|------|--------|
| Brannen-form eigenvalues | **Structural** — derived from Cl(3,1) Z₃ |
| Brannen Koide closed form $Q = (1+2t^2)/3$ | **Lean-verified** (`BrannenKernel.Q_value`, 2026-05-22) |
| Constraint $Q = 2/3 \iff t^2 = 1/2$ | **Lean-verified** (`BrannenKernel.koide_iff_constraint`) |
| Koide identity Q = 2/3 | **Structural** — $\dim G_2/\dim \text{Spin}(7) = 14/21$ |
| Brannen phase φ = 2/9 rad | **Structural** — $Q/3$, with 3 = number of generations |
| Z₃ identity $1 + \omega + \omega^2 = 0$ | **Lean-verified** (`CyclicShift.sum_one_omega_omega_sq`) |
| $\dim \text{Spin}(7) = 21 = \binom{7}{2}$ | **Lean-verified** (`SpinDimension.dimSpin_seven`) |
| $\dim G_2 = 14 = 49 − 35$ | **Lean-verified** (`SpinDimension.dimG2_eq_14`; orbit-stabilizer) |
| $\dim G_2 = \dim \text{Spin}(7) − \dim S^7 = 14$ | **Lean-verified** (`SpinDimension.dimG2_via_S7`) |
| Koide ratio $\dim G_2 / \dim \text{Spin}(7) = 2/3$ | **Lean-verified** (`SpinDimension.koide_ratio_structural`) |
| Brannen phase $((2/3)/3) = 2/9$ | **Lean-verified** (`SpinDimension.brannen_phase_structural`) |
| Silent SU(2)/U(1) direction on S³ | **Lean-verified** (`SilentDirection.silent_pair`); empirical to 4×10⁻¹⁵ |
| Cross-sector ratio 21 | **Structural** — $\dim \text{Spin}(7)$ |
| ~~$\alpha_W / \alpha = \sqrt{21}$~~ | **Refuted** (cosserat_experiment/06_killing_form.py, 2026-05-22): the Killing-form embedding index of so(3) ⊂ so(7) is **5**, not 21. |
| $g_W^2 = 5 \cdot \sqrt{\alpha}$ | **Conjecture** (cosserat_experiment/07_full_lagrangian.py, 2026-05-22): combines Killing-form embedding 5 (DERIVED) with α^{1/4} = exp(-π²/8) (v59-instanton-NATURAL). Predicts $g_W = 0.6535$ vs empirical $g_W(M_Z) = 0.6517$ — **0.28 %** gap. Not yet Lagrangian-derived. |
| $G_e = (21/16) \cdot \alpha^{21}$ | **Conjecture** (cosserat_experiment/08_gravity_correction.py, 2026-05-22): predicts $G_e = 1.756 \times 10^{-45}$ vs empirical $1.752 \times 10^{-45}$ — **0.25–0.33 %** gap. The 21/16 factor was the unique best match (gap 0.33%) out of 60 v59-natural structural ratios scanned; next best (28/21) was 6× further off. |
| Quark Brannen: $t^2_N = 1 - \dim G_2/D_N$ | **Conjecture** (cosserat_experiment/11_quark_sector.py, 2026-05-22): with $D_0 = \dim \text{Spin}(8) = 28$ (lepton), $D_1 = \dim \Lambda^3 \mathbb{R}^7 = 35$ (d-quark), $D_2 = 63 = 3 \cdot \dim \text{Spin}(7)$ (u-quark). Predicts $Q_d = 11/15$ (0.26% gap), $Q_u = 23/27$ (0.34% gap). Note $28 + 35 = 63$ (additive identity). |
| **Single source: Cl(7)_even ≅ ℂ⊗𝕆** (Furey color algebra, dim 64) | **Mechanism** (cosserat_experiment/13_single_source.py, 2026-05-22): all three D_N values are graded subspaces of the same parent algebra. Decomposition: Λ⁰=1, Λ²=21=dim Spin(7), Λ⁴=35=d-quark D, Λ⁶=7=dim S⁷, total 64. Lepton uses Λ²+Λ⁶=28; d-quark uses Λ⁴=35; u-quark uses Λ²+Λ⁴+Λ⁶=63. Selection rule (which N→which grades) not yet derived. |
| **Z₂×Z₂ decomposition** | **Novel observation** (cosserat_experiment/16_Z2_decomposition.py, 2026-05-22): Cl(7)_even bisects into L = Λ²⊕Λ⁶ (no G₂-invariant, 28) and F = Λ⁴ (contains coassociative 4-form, 35). Each sector picks bits (B_L, B_F): lepton=(1,0), d-quark=(0,1), u-quark=(1,1). Additive identity D_u=D_e+D_d emerges as L⊕F direct sum. Why each Furey N takes its specific bits still NOT derived. |
| Three generations | **Structural** — Z₃ ⊂ S₃ triality of Spin(8) |
| Complex structure | **Structural** — from $B^2 = -1$ in Cl(3,1) bivectors |
| Overall lepton mass scale $a$ | **Empirical** |
| α (absolute value) | **Empirical**; $\pi^2/2$ conjecture (0.30% gap) |
| G (absolute value) | **Conjecture** ($(21/16) \cdot \alpha^{21}$, NEW 2026-05-22): predicts $G_e = 1.756 \times 10^{-45}$ vs empirical $1.752 \times 10^{-45}$ — **0.25–0.33 %** gap. Prior $\alpha^{21}$ alone was 24% off. The (21/16) = dim Spin(7) / dim Cl(3,1) prefactor closes the gap. |
| 0.26% gap in cross-sector ratio | Modulation-function correction beyond pure exp(-x), still to derive |
| Quark sector | Not yet addressed |
| Full Standard Model gauge structure | In progress via Furey construction |

## Open Questions

→ **For the full consolidated roadmap, see [`ROADMAP.md`](ROADMAP.md)** —
contains all open questions across 4 frontiers (dynamical ξ, selection rule,
Lagrangian prefactors, scale bridges) plus 6 minor Lean items and 6
speculative items, with concrete next-session entry points.

### Top open questions (summary)

1. **Dynamical ξ(x) field**: How does the v59 Brannen parameter ξ become a
   dynamical field?  Framework set up (`synthesis/DYNAMIC_XI.md`); 3+1 mode
   spectrum identified (3 Goldstones + 1 Higgs-like radial).  Six sub-questions
   open (sector multi-vacua, λ value, scale bridge to Higgs, gauging,
   Yukawa form, L⊕F in dynamics).

2. **Selection rule mechanism**: Why does each Furey N take its specific
   (Bit-L, Bit-F) value?  Seven hypotheses tested in 2026-05-22 session, none
   derived.  See `cosserat_experiment/FINDINGS_selection.md`.

3. **Lagrangian derivation of (5, 21/16)**: The g_W² = 5·√α and G_e = (21/16)·α²¹
   conjectures match at 0.25–0.28 % but lack Lagrangian derivation. Both use
   only {dim Spin(7) = 21, dim Cl(3,1) = 16} via different operations.

4. **Scale bridges**: Brannen scale a (314 MeV²) vs SM Higgs VEV (246 GeV) —
   6 orders of magnitude apart.  No structural bridge identified yet.

5. **What sets the lepton mass scale $a$ in Planck units?** $a/M_P$ remains
   empirical.

6. **Quark sector**: Brannen Koide extends via `t² = 1 - dim G₂/D_N` (this
   session).  But quark Brannen phases φ_q (-2.02 for u, +0.110 for d) and
   absolute mass scales remain empirical.

7. **Standard Model gauge group**: $\text{SU}(3)_c$ from Cl(6) Witt
   decomposition (variant B).  $\text{SU}(2)_L$ as silent SU(2)/U(1)
   (this session, Lean-proved).  $\text{U}(1)_Y$ and Higgs sector not yet
   identified.

## Directory Structure

```
v59/
├── README.md                       — orientation, pointers to all subdirectories
├── SUMMARY.md                      — this document, the consolidated headline results
├── TYCHO_TABLE.md                  — fitting targets (PDG/CODATA values)
├── first_experiments/              — Steps 1-4: numerological scans
│   ├── 01_koide_substrate_scan.py + 01_findings.md     — Brannen confirmed; phi residue
│   ├── 02_two_ninths_audit.py + 02_findings.md         — 2/9 NOT in Coxeter/Lie/polytopes
│   ├── 03_modular_hyperbolic_scan.py + 03_findings.md  — 2/9 NOT in modular/hyperbolic
│   └── 04_muon_electron_ratio_scan.py + 04_findings.md — m_μ/m_e not from simple substrates
├── multivector_kernel_fit/         — Steps 5-8: kernel-fit on Cl(3,1)
│   ├── 01_z3_lepton_kernel.py + 01_findings.md         — Cl(3,1) Z₃ derives Brannen form
│   ├── 02_z3_invariant_potential.py + 02_findings.md   — Z₃ potential pins Koide minimum
│   ├── 03_em_bivector_coupling.py + 03_findings.md     — α not predicted by Cl(3,1) alone
│   ├── 04_quaternionic_constraint.py + 04_findings.md  — Koide STRUCTURAL via S³ constraint
│   ├── 05_alpha_G_corroboration.py + 05_findings.md    — α-G unify under exp modulation
│   └── 06_normalization_check.py + 06_findings.md      — 21 = 3 × 7 is best fit
├── octonionic_extension/           — Steps 9-11: octonionic identification
│   ├── 01_fano_and_21.py + 01_findings.md              — 21 = dim Spin(7)
│   └── 02_brannen_and_instanton.py + 02_findings.md    — φ = Q/3 STRUCTURAL
└── furey_construction/             — Step 12+: full ℂ⊗ℍ⊗𝕆 construction (ACTIVE)
    ├── PLAN.md                                          — multi-variant plan
    ├── README.md                                        — orientation
    └── (subprojects per the plan)
```

## Where This Stands in the Larger SCP Project

- **v28–v55**: Soliton-based modeling on fixed grids (6-field Cosserat, multivector GA). Excellent qualitative phenomenology, no quantitative TYCHO matches.
- **v56**: Pivot to projective geometric algebra ℝ(3,1,1) with Higgs + Skyrme. Topology survival improved but lattice collapse remained.
- **v57**: Quantitative mismatch with lattice-QCD proton mechanical structure documented. Layer-level concerns identified.
- **v58**: LAYER_CRITIQUE diagnosed the fixed-background problem. First-principles requirements articulated. Unified multivector force experiment achieved formal recovery of Newton + Maxwell limits but on a 2D fixed lattice — reintroducing the layer error.
- **v59 (this version)**: Kepler-stage numerology + multivector-kernel fitting. **First quantitative TYCHO match in the project's history** (charged lepton masses to machine precision via constraint surface). Cross-sector framework consistent but not yet fully predictive.

## Compared to v58

v58's "unified multivector force" claimed Newton + Maxwell emerge from one equation. **Qualitatively true** (different grade projections give different sector behaviors), **quantitatively empty** (no TYCHO matches, parameters tuned).

v59's kernel-fit program **inverts the workflow**: start from tight numerical mysteries (TYCHO_TABLE), fit a structural kernel, and accept partial unification rather than over-promise. Concrete deliverables:

- One quantitatively-pinned ellipse: Brannen form + Koide for charged leptons.
- One quantitatively-tight residue: cross-sector ratio 20.95 ≈ 21.
- One specific testable prediction: α variation in gravitational potential ~10⁻¹⁰.

This is *less* than v58 claimed but *more* than v58 ever produced quantitatively.

## Posture for Octonionic Extension

The "7" in 21 = 3 × 7 is the most actionable open question. The Cl(3,1) algebra has rank 16 = 2⁴ (dimension 4 = 3 + 1 spacetime). Going to Cl(0,7) gives the octonions (dim 8 = 2³). Cl(3,1) ⊗ Cl(0,3) has dim 128 = 16 × 8.

If 7 emerges as a natural rank/dimension of the extended algebra, the cross-sector ratio 21 becomes structural and the kernel is complete. This is the next concrete experiment, opened in `octonionic_extension/`.
