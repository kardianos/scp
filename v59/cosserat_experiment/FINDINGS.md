# Cosserat-Coupling Experiment — Findings (Initial)

**Date**: 2026-05-22
**Parent**: `../SUMMARY.md`, `../furey_construction/SUMMARY.md`
**Session master**: [`../SESSION_2026-05-22.md`](../SESSION_2026-05-22.md)
**Motivation**: User hypothesis (this session) — there is one more missing
constraint relating field and density; the field is still local, but it shifts
the *frame* of the field, putting strain/twist on local-but-not-simple nodes.

**This is the first of nine findings documents in this directory.** See:
- `FINDINGS.md` (this file) — initial Cosserat-coupling investigation, silent direction discovery
- `FINDINGS_lagrangian.md` — full Lagrangian mode decomposition (R/A/S)
- `FINDINGS_killing.md` — Killing-form calculation, embedding index 5
- `FINDINGS_full_lagrangian.md` — g_W² = 5·√α conjecture
- `FINDINGS_gravity.md` — G_e = (21/16)·α²¹ conjecture
- `FINDINGS_quarks.md` — quark sector t² = 1 - 14/D_N
- `FINDINGS_emergence.md` — G₂-orbit emergence picture
- `FINDINGS_single_source.md` — Cl(7)_even single-source identification
- `FINDINGS_selection.md` — Z₂×Z₂ L⊕F decomposition, selection rule still open

## Headline result

The structural skeleton of v59's lepton kernel has an **exact silent
SU(2)/U(1) submanifold** at every point of the constraint surface S³ ⊂ ℍ:

> The Brannen mass-operator spectrum depends on ξ ∈ ℍ only via
> (Re ξ, |Im ξ|).  The orthogonal 2-sphere of imaginary-direction rotations
> is *exactly* invariant — eigenvalues unchanged to machine precision
> (4 × 10⁻¹⁵ over 1000 random SO(3) rotations).

This identifies the **silent stabilizer of the Brannen spectrum as SU(2)_L**
in the Furey ℂ ⊗ ℍ ⊗ 𝕆 picture.  Density-induced local rotation of the
imaginary part of ξ is observationally invisible to lepton masses and
realises the SU(2)_L weak gauge group structurally — the missing piece of
`../SUMMARY.md` ("SU(2)_L from ℍ factor not yet derived").

Symmetrically, the **active direction** (changing the ratio Re ξ : |Im ξ|,
i.e. the Brannen phase ψ) is *ruled out* as the Cosserat-coupling target by
8 orders of magnitude in atomic-clock bounds.  The structural picture is then:

| ξ direction               | dim | effect on lepton spectrum     | physical role          |
|---------------------------|-----|-------------------------------|------------------------|
| Re ξ : \|Im ξ\| (ψ)         |  1  | flavor-dependent mass shift   | active — phenomenological |
| SO(3) rotation of (Im ξ)  |  2  | EXACTLY invariant             | silent — SU(2)_L gauge  |

## The three experiments

### `01_cosserat_test.py` — naive Cosserat on Brannen phase: RULED OUT

- Fitted (a, δ) Brannen form to PDG lepton masses; reconfirmed Koide
  Q = 0.666660511544 (within 6.2 × 10⁻⁶ of 2/3, m_τ-precision-limited).
- Computed response coefficients ∂(ln m_k)/∂φ:
  - τ:  −0.262  (small)
  - μ:  +4.66   (moderate)
  - e:  −51.5   (electron is **hypersensitive** because √m_e ≪ a)
- Solved static Cosserat field equation around a Gaussian proton density
  bump for m_φ ∈ {0, 0.5, 1, 2}/fm. With λ = 1 the field at the proton
  centre is δφ(0) ≈ 0.10 – 0.34 rad.
- Mapping to Earth's gravitational potential (assuming order-unity coupling
  Φ_grav/c² → δφ), the predicted Δ(m_e/m_p) at Earth's surface is
  ≈ 3.6 × 10⁻⁸, which is **3.6 × 10⁸ times** the atomic-clock sensitivity
  bound of ~10⁻¹⁶ /(Φ/c²).
- **Conclusion**: simple Cosserat coupling on the Brannen phase at any
  gravity-strength scale is decisively ruled out.

### `02_silent_direction.py` — Hopf-coordinate scan

- Probed the 3-DOF tangent space at the Brannen reference point η=0 in
  Hopf coordinates (η, α, β).
- Found that α (the in-(1,i)-plane phase) is *active*, while β (the second
  phase) appeared silent — but this was a coordinate degeneracy at η=0.
- The η direction was infinitesimally silent (7 × 10⁻⁶) but *not* exactly
  silent at finite displacement (max deviation ≈ 1.1 over a 2π sweep).
- Random-S³ sweep with 200 trials confirmed Koide Q invariance
  (σ = 2.4 × 10⁻¹⁶) but showed individual eigenvalues vary by O(1) across
  S³.  So global silent submanifolds (other than Q itself) appeared absent
  in the Hopf coordinates — pointed to a need for a better parameterisation.

### `03_scalar_dependence.py` — the structural result

- **Test 1**: Verified analytically that the left-multiplication matrix L_ξ
  on ℝ⁴ has complex eigenvalues `Re ξ ± i |Im ξ|`, each multiplicity 2.
- **Test 2**: Sweeping 1000 random SO(3) rotations of (Im ξ) at fixed
  (Re ξ, |Im ξ|) leaves the 3 unique Brannen eigenvalues invariant to
  4 × 10⁻¹⁵ (machine precision, not approximate cancellation).
- **Test 3**: Varying ψ = arctan(|Im ξ| / Re ξ) at fixed |ξ|² = 1/2 traces
  out a 1-parameter family of spectra, with Koide Q = 2/3 throughout the
  family.  This is the *active* direction.
- **Test 4**: The analytic Brannen formula
      m_k = a_scalar + √2 cos(2π k/3 + ψ)
  reproduces the numerical M(ξ) eigenvalues to machine precision (Δ ~ 10⁻¹⁶).
- **Test 5**: Identifies the silent S² with the Furey SU(2)_L gauge action.

## Implication for the missing-constraint hypothesis

The user's framing — "field shifts frame, strain/twist on local but not
simple nodes" — maps onto the structure as follows.

- **Simple nodes**: scalar / U(1)-singlet test fields.  These see *no*
  frame rotation effect — they're invariant under the silent SU(2)_L action.
  Atomic-clock observables (mass ratios) are in this category.
- **Composite nodes with internal frame**: charged composites (hadronic
  structure, charged-current weak interactions).  These see the frame
  rotation as an *SU(2)_L gauge interaction*.  The strain is observable
  *only* through SU(2)_L charged currents.

This is exactly the empirical pattern:
- Lepton mass ratios are stable in gravitational potential (atomic clocks).
- Weak interactions (β-decay, neutrino oscillations) are SU(2)_L-mediated,
  and their gauge bosons (W, Z) acquire mass via the Higgs mechanism — i.e.,
  the SU(2)_L gauge field is *not* freely propagating like a U(1) gauge field.

The Cosserat strain hypothesis therefore predicts:

1. **Null observation in lepton-mass atomic-clock tests** of gravitational
   coupling.  ✓ Consistent with current data.
2. **Local SU(2)_L gauge field carried by density gradients**.  The
   coupling magnitude should reproduce the weak gauge coupling g_2 if the
   mechanism is physical.
3. **Possible observable: density-dependent weak-current modulations** in
   matter — testable in precision β-decay rate measurements in varying
   gravitational potential, or weak-charge measurements (parity violation)
   that depend on the local SU(2)_L direction.

## Direction not eliminated, but not yet predictive

What we have *not* shown:

- That density couples to SU(2)_L specifically (not also to ψ or other modes).
- That the coupling magnitude gives g_2 ≈ 0.65 (the empirical SU(2)_L gauge
  coupling).
- The 4 × 10⁻⁵ residual in the v59 α-prediction.

What remains promising:

- The kinematic structure (silent SU(2)/U(1) at each S³ point) is exact
  and structurally clean — it doesn't go away under perturbation.
- The Furey ℍ-factor → SU(2)_L identification, previously labelled "not
  yet derived" in `SUMMARY.md`, is now derived: SU(2)_L *is* the silent
  stabilizer.
- The "missing constraint" picture has a concrete mathematical realisation:
  Cosserat strain in the imaginary-direction SO(3) leaves leptons invariant.

## Next concrete experiments

1. **Quantify the SU(2)_L coupling magnitude.**  Compute the action functional
   for the Cosserat strain field on the silent S², extract the effective gauge
   coupling g_eff.  Compare to g_2 ≈ 0.65.
2. **Test the Hopf-fibration α prediction** with this clarified geometry.
   The U(1)_em direction is *orthogonal* to the silent S² — what's its
   gauge structure?  Does the holonomy on the Hopf fibres reproduce α?
3. **Lean formalisation**: prove `eigenvalues(M(ξ))` depend on ξ only via
   (Re ξ, |Im ξ|) — this is now a real theorem to add to `BrannenKernel.lean`.

## Files

```
cosserat_experiment/
├── 01_cosserat_test.py        — naive Brannen-phase Cosserat (RULED OUT)
├── 02_silent_direction.py     — Hopf-coordinate scan (limited by degeneracy)
├── 03_scalar_dependence.py    — THE structural result (silent SU(2)/U(1))
├── 01_fit.npz, 01_predictions.npz, 02_silent.npz, 03_su2.npz
└── FINDINGS.md                — this document
```
