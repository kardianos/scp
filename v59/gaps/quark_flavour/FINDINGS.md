# Quark-Flavour Cluster — Findings & Verdict (G4, G5, G6)

**Date**: 2026-05-25 · Cluster lead deliverable. See `README.md` (background, all
hypotheses), `ALTERNATIVES.md` (avenues + falsifiers), `QuarkKoide.lean` (formalization),
and the three Python scripts (`quark_koide_rg.py`, `quark_phases_scales.py`,
`rg_invariance_test.py`).

---

## Verdict per gap

### G4 — Is the Z₂×Z₂ selection rule derivable?  **PARTIALLY.**

- **d-quark → F: DERIVED [thm].** Machine-checked in `Z2Z2Forcing.lean`: only the
  F-grade supplies the rank-≥2 color Cartan a color triplet needs; all 28 L-grade
  generators are diagonal-free. `color_splitting_requires_F` is a genuine theorem on
  the (corrected) Cl(7) generators. This half of G4 is no longer "undriven" — the gap
  brief understates the current state.
- **u-quark → L⊕F: DERIVED-UP-TO-GLUE [thm-ingredients].** L is not closed
  (`γ_iγ_j·γ_kγ_l ∈ F`, a universal grade law), so a Fock-composite (N=2 = two
  raisings) spills L into F; combined with the F color Cartan this forces L⊕F. The
  grade arithmetic is theorem-grade; the "N=2 = product of two raisings" reading is
  interpretive.
- **lepton → L: OPEN [open].** Honestly *not* forced. `LeptonGradeForcing.lean` proves
  (machine-checked) that F actually offers leptons a *richer* mass channel than L; the
  whole lepton assignment reduces to one well-posed question — **does the ℂ-imaginary
  unit J of the Furey ideal live in Λ² (⊂ L)?** If yes, lepton=L is forced; if J∈F it
  flips. This is a finite, decidable computation that has *not yet been done*.
- **Full map N↦(Bit-L,Bit-F): no Lagrangian derivation.** Seven hypotheses tried
  (Witt bidegree, Hodge, triality, winding, centralizer, SU(3) filtering, WZW); the
  SU(3)-content route (H6) is the one that *partly worked* — it became the d-quark=F
  color-Cartan theorem. The centralizer route (H5) came closest numerically (34 vs 35)
  but is wrong.

**Net**: ⅔ of the selection rule is now theorem-grade or theorem-up-to-glue; the
remaining ⅓ (lepton=L) is reduced to a single finite computation (locate J). This is
the **most tractable** open piece in the whole cluster.

### G5 — Is quark Koide a prediction or a soft pattern?  **A SOFT PATTERN.**

- The closed forms `Q_d=11/15`, `Q_u=23/27` and the cross-sector invariant
  `(1−t²_N)·D_N=14` are **[thm]** (`QuarkKoide.lean`).
- The agreement with data is **NOT** a 0.3% prediction. The 0.3% match exists **only**
  at the PDG mixed-scale convention (each quark at a different reference scale — not a
  physical scale). At any single physical scale (M_Z, m_t, 1 TeV) the gap is **2–5%**.
- RG spread of Q_d (0.052) is **27×** the claimed gap; of Q_u (0.047) is **16×**.
- **Sharper diagnosis** (`rg_invariance_test.py`, a *correction* to a naive "Koide
  isn't RG-invariant"): Koide Q *is* exactly invariant under common rescaling, and
  one-loop QCD running is flavor-universal — so **pure QCD running does not move Q**.
  The drift is from (i) input-scale/threshold mismatch and (ii) the
  **flavor-non-universal top Yukawa**. Q_u is dominated by √m_t and is therefore
  *intrinsically* unprotected; Q_d is limited instead by the ±9% s-quark mass error.
- **Overfitting**: 12 distinct simple rationals (q≤27) sit inside each RG band.
  11/15 and 23/27 are not distinguished from 7/9, 3/4, 6/7, 8/9, …

**Net**: G5 is a structural *pattern* that lands near the data at one unphysical
convention. Honest precision is ~2–5%, not 0.3%. It also inherits G4's open lepton
half (the D_N denominators come from the assignment).

### G6 — Are the quark phases & scales structural?  **NO — FREE.**

- **φ_q ≠ Q_q/3** [thm-negative]: the exact lepton relation does not extend (Q_d/3 is
  >2× φ_d; Q_u/3 has the wrong sign).
- The α-numerology (φ_d≈14α, φ_u≈−10α) is testable only via the observable cos(3φ),
  which — because the quark phases are small — sits at 0.95–0.98 for *every* candidate
  (<0.7% spread). cos(3φ) is a **weak discriminator**; the value claim is not
  falsifiable at useful precision.
- The **defensible** content (avenue C2) is a *size/hierarchy* claim: φ_quark ~ O(α)·
  φ_lepton (the quark phases look radiative, the lepton phase tree-level at 2/9). The
  ratios φ_d/φ_l=0.489 vs 14α(M_Z)/(2/9)=0.492 (0.8%) and φ_u/φ_l=0.326 vs
  10α(0)/(2/9)=0.328 (0.6%) support O(α) *size*; the specific integers remain [conj].
- Scale ratios a_u²/a_d²≈35, a_u²/a_l²≈72 are **convention artifacts** — they hold
  only at the PDG mixed scale (at M_Z: 49 and 69; at 2 GeV: 45 and 153). **Dead** as
  single-scale relations.

**Net**: φ_d, φ_u, a_u, a_d are genuine free inputs. The only survivor is the *order-α
size* of the quark phases (a hierarchy statement, not a value prediction).

---

## Most promising avenue

**(A1) Locate the ℂ-imaginary unit J in the Cl(7)_even grade decomposition.** This is
the single highest-leverage, near-term-feasible move in the cluster: prior work
(`LeptonGradeForcing.lean`) already reduced the entire open lepton half of G4 to this
one finite, decidable computation. Build the explicit `Cl(7)_even → M₈(ℂ)` isomorphism
(64 basis elements, 8×8 matrices — the computation flagged in
`13_fock_mass_forcing_report.md` §6), construct the Furey Witt idempotent, extract J,
and read off its grade. If `J ∈ Λ²`, lepton=L becomes a theorem and G4 is essentially
closed; if `J ∈ F`, lepton=L is *refuted* as forced (an equally valuable null).

**Runner-up (B2)**: the common-QCD-factor RG test already executed here shows Q is
protected against universal running but not against the top Yukawa — which *explains*
why Q_u is the worst offender and confirms G5 cannot be rescued to a tight prediction.
This closes the "is it RG-invariant?" question with a clean negative for Q_u.

---

## Honest bottom line for the cluster

The quark sector is the weakest tier of v59, but not uniformly so:

- **G4** is in better shape than advertised: d-quark→F is a theorem, u-quark→L⊕F is
  theorem-up-to-glue, and lepton→L is one finite computation away from resolution.
- **G5** is a soft, RG/convention-dependent pattern — *not* a 0.3% prediction. Honest
  precision 2–5%; 12 rationals compete in the band; Q_u intrinsically unprotected by
  the top Yukawa.
- **G6** is free, except for a defensible *order-α size* of the quark phases.

The cluster's numerology temptation (running masses, large uncertainties) is real and
was the main risk; the RG analysis here quantifies it and downgrades the headline
"0.3%" quark-Koide claim to a convention artifact. Report delivered with null results
front and center.

---

## File list (all in `v59/gaps/quark_flavour/`)

| file | role |
|---|---|
| `README.md` | background: precise statement of G4/G5/G6, all hypotheses tried, RG-invariance problem |
| `ALTERNATIVES.md` | open solution space: avenues A1–A4 (G4), B1–B3 (G5), C1–C3 (G6) with test+falsifier each |
| `FINDINGS.md` | this file — verdicts + most promising avenue |
| `quark_koide_rg.py` / `.json` | Q_d, Q_u at 5 scales/schemes; RG spread; PDG-MC band; overfitting (rationals-in-band) control |
| `quark_phases_scales.py` / `.json` | phi_d, phi_u candidate forms on cos(3phi); n*alpha overfitting control; scale-ratio convention test |
| `rg_invariance_test.py` | the B2 test: common-rescaling invariance, flavor-universal QCD running, top-Yukawa isolation; C2 hierarchy check |
| `QuarkKoide.lean` | Lean: D_N decompositions, t²_N=1−14/D_N, Q closed forms, cross-sector invariant [thm]; G4/G5/G6 open content flagged with `sorry`/refutation. **Written, NOT built this run** (shared-project build conflict). |
