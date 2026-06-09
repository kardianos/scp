# v60 PLAN — Multi-Variant Attack on G9 (Induced Metric) and G1 (Rank Tension)

**Date**: 2026-05-25 (kickoff)
**Parents**: [`README.md`](README.md), [`ROADMAP.md`](ROADMAP.md), [`../v59/CLOSEOUT.md`](../v59/CLOSEOUT.md), [`../v59/gaps/SYNTHESIS.md`](../v59/gaps/SYNTHESIS.md)
**Goal**: Execute bounded, exhaustive, documented attacks on the two highest-priority gaps. Produce either (a) working resolutions that enable the dynamical Lagrangian deliverable, or (b) the strongest possible negative results with clear restructuring implications for the program.

**Discipline**: No giving up on a variant until it has been run to completion or honest documented failure. Every ansatz gets a `0X_findings.md` (or Lean theorem) that states what was tried, what the outcome was, and why it succeeded or failed. Clean falsifications are first-class outputs.

---

## Overall Attack Structure

Two parallel but communicating tracks:

1. **G9 Track** (`v60/gravity_recast/`) — Induced-metric / emergent-gravity recast.
2. **G1 Track** (`v60/gaps/rank_tension/` or `v60/synthesis/rank_tension/`) — Two-piece Y + breaking-chain resolution.

The tracks share the synthesis / consistency layer (`v60/synthesis/`) and the Lean corpus (`v60/lean/` or extensions of the v59 tree).

A third gated track (Lagrangian) begins only when at least one of G9 or G1 has a live candidate.

---

## G9 Track — Induced Metric Recast (Plebański-Style and Variants)

**Reference implementations / literature anchors**:
- Plebański 2-form gravity (self-dual 2-forms, simplicity constraints).
- Ashtekar variables / loop quantum gravity simplicity constraints (for inspiration only — we are not doing LQG).
- Existing v59 polarization analysis (`g9_polarization_test.py`, `G8G9_Gravity.lean`).

**Success criterion**: A concrete ansatz in which a fundamental 2-form (derived from or projected from the v59 \(\Omega \in \Lambda^2\)) plus a soldering / simplicity constraint yields a symmetric spacetime tensor with *exactly* two transverse-traceless physical degrees of freedom in the linearized spectrum around flat space, while the long-range force still couples to the second-moment charge \(\rho_{\rm grav}\) at the observed strength (or a clean, bounded search that demonstrates no such natural ansatz exists within a stated family).

### Variant Family A — Direct Internal 2-Form

- **A1**: \(\Omega\) itself (or its self-dual projection) is the fundamental Plebański 2-form \(B^{ab}\). Simplicity constraint of the form \(B^{ab} \wedge B^{cd} \propto \epsilon^{abcd}\) or the v59-analog using the internal metric on \(\Lambda^2\).
- **A2**: Same as A1 but with an additional "soldering" map that identifies a preferred Lorentz subalgebra inside the internal \(\mathfrak{so}(8)\).
- **A3**: Non-chiral (full) 2-form version; allow both self-dual and anti-self-dual sectors.

**Expected output per sub-variant**: `gravity_recast/A1_plebanski_internal.py`, `A1_findings.md`, and (where possible) a Lean statement of the constraint and the helicity content.

### Variant Family B — 2-Form Valued in the Full Algebra

- **B1**: The fundamental object is a 2-form valued in the whole \(\mathrm{Cl}(7)_{\rm even}\) (or the L⊕F bisection), projected onto the \(\Lambda^2\) sector that sources gravity.
- **B2**: Use the coassociative 4-form or other structural elements of the algebra to define the simplicity constraint (exploiting the 35 = F grade that was "left out" of leptons).
- **B3**: Triality-rotated versions (different embeddings of the 2-form under the Z₃ action).

### Variant Family C — Auxiliary Soldering Field (Tetrad or Frame)

- **C1**: Introduce an auxiliary soldering field (tetrad \(e^a_\mu\)) that couples the internal 2-form to spacetime indices; integrate it out or constrain it to be non-dynamical.
- **C2**: Dynamical tetrad + 2-form system (inspired by first-order formulations of GR); check whether the v59 scalar charge emerges as the source for the emergent metric.
- **C3**: Purely algebraic soldering (no new fields) using the existing grade structure and the complex structures \(J \in L\).

### Variant Family D — Limits and No-Go Theorems (Defensive)

- **D1**: Linearized analysis around flat space for the broadest natural class of 2-form + constraint systems compatible with the v59 algebra; prove a general helicity content theorem (reuse/extend the little-group SO(2) decomposition from `g9_polarization_test.py`).
- **D2**: "Why 2 TT DOF are impossible" — exhaustive scan of low-dimensional ansatz families with symbolic or numeric mode counting; document the representation-theoretic obstruction if one exists.
- **D3**: Short-range tensor mode only — if the tensor content is massive or confined, show the long-range force remains scalar (and therefore still fails LIGO at observable distances).

**Order of attack**: A1 → A2 → B1 → C1 (the most "natural" first). Defensive D variants run in parallel from the start to bound the search.

**Lean targets for G9**:
- `InducedMetric.lean`: statements of "2-form + simplicity constraint → helicity content", "soldering exists", "DOF count = 2".
- Reuse and extend `G8G9_Gravity.lean` (gauge_index_additive, gravity_prefactor, helicity lemmas).
- New: formalization of the Plebański action and the constraint in a geometric algebra setting (may require Mathlib developments or a lightweight custom theory).

**Numerical / symbolic targets**:
- `01_plebanski_setup.py`: SymPy linearization of the 2-form fluctuations, extraction of the physical mode spectrum.
- Polarization test harness (fork of `g9_polarization_test.py`) that ingests a candidate constraint and returns the helicity content + coupling to \(\rho_{\rm grav}\).
- Scale and magnitude checks: does the emergent \(h_{\mu\nu}\) reproduce the \(\alpha^{21}\) strength in the weak-field limit?

**Documentation per variant**:
- `0X_<variant>_setup.py` + `0X_findings.md` (what was tried, exact equations, numerical or symbolic outcome, why it succeeded/failed or was abandoned).
- If a no-go is proved: a short Lean theorem or a clear "representation obstruction" note.

**Exit criteria for the G9 track**:
- A live candidate with 2 TT DOF + correct long-range charge coupling, or
- A bounded search (e.g., 4–6 major families, 2–3 weeks of effort) that yields a strong negative result with a recommended pivot (e.g., "gravity in this algebra is short-range only", "tensor modes require a different carrier", "the program must modify GR or drop long-range gravity").

---

## G1 Track — Rank Tension and Two-Piece Y (Parallel Track)

**Reference material**:
- `v59/gaps/ew_scale_bridge/formalize_bridge.py`, `EwScaleBridge.lean` (the 784 = dim End(L) theorem, R1/R2 statements).
- `v59/synthesis/gravity_charge_test.py` (the 6/784 ratio and Frobenius² unification).
- Representation theory of \(\mathfrak{so}(8)\) (already partially formalized in the v59 Lean corpus: `SpinDimension.lean`, `LieDimensions.lean`).

**Success criterion**: An explicit decomposition or breaking chain inside (or acting on) the 28-dimensional space such that (a) the full 784-dimensional End(L) structure is preserved for the bridge scale, (b) exactly three (or a small number of) directions acquire the physical Brannen eigenvalues and source \(\rho_{\rm grav}\), and (c) the 0.068% numerical match is either preserved or explained as a derived consequence rather than an input.

### Variant Family A — Explicit Block Decomposition of End(L)

- **A1**: Direct 28 × 28 matrix split into a 3 × 3 physical block + 25 × 25 (or 3+25 mixed) Goldstone block, with the off-diagonal blocks constrained by the algebra.
- **A2**: Use the known complex structure \(J \in L\) and the color \(\mathfrak{su}(3)\) action to define invariant subspaces.
- **A3**: Generation-space democracy vs. sector democracy — different splits for leptons vs. quarks that still respect the additive identity \(D_u = D_e + D_d\).

### Variant Family B — Subalgebra Chain (the "chain, not single step")

- **B1**: Exhaustive classification of maximal proper subalgebras of \(\mathfrak{so}(8)\) (dim 21, 16, …) and all chains down to a 3-dimensional (or rank-3) stabilizer. Check which chains can leave exactly three light directions.
- **B2**: Map each step in the chain onto the already-proved Spin(7) → Pati-Salam → SM breaking pattern.
- **B3**: Count Goldstone bosons at each step and check consistency with the silent SU(2)_L and the would-be Goldstones from the Higgs mechanism.

### Variant Family C — Two-Piece Y with Shared Scale Conjecture

- **C1**: Treat the 784 count as applying to the *total* Frobenius norm across both pieces; the active rank-3 piece sources gravity while the Goldstone piece contributes only to the bridge "background".
- **C2**: Derive the 6/784 ratio from the representation theory of the split rather than imposing it.
- **C3**: Extend the split to the quark sectors (35 and 63) and check consistency with the observed mass hierarchies and the universal Koide deviation \((1 - Q_N) D_N = 28/3\).

### Variant Family D — No-Go / Alternative Interpretations (Defensive)

- **D1**: Prove (or numerically demonstrate to high confidence) that no 25-dimensional subalgebra of \(\mathfrak{so}(8)\) can be the "heavy" complement to a 3-dimensional light sector.
- **D2**: "The 784 bridge is not End(L)" — explore whether the operator algebra that resolves non-associativity is smaller than the full End(L) once the two-piece structure is taken into account.
- **D3**: Composite or effective interpretation of the bridge (the 784 is an *emergent* count after integrating out heavy modes).

**Order of attack**: B1 (chain classification) + A1 (explicit split) in parallel, then C1 once a candidate split exists. Defensive D variants from day one.

**Lean targets for G1**:
- Extend `EwScaleBridge.lean` (or new `RankTension.lean`): formal statement of the tension (R1/R2), proof that 784 = dim End(L) is generic for the adjoint, proof of absence of suitable 25-dim subalgebras (or the positive statement of a working chain).
- Reuse `HiggsVevReframe.lean`, `XiVacuum.lean` (the vacuum that enforces the Koide constraint surface).
- New: theorems about the Frobenius norm on the split space and the resulting gravity charge.

**Numerical targets**:
- `01_two_piece_y.py`: scan of splits, computation of \(\|Y_{\rm rank3}\|_F^2 / \|Y_{\rm total}\|_F^2\), check against the empirical 6/784 for leptons and the corresponding numbers for quarks.
- Update to `gravity_charge_test.py` that ingests a candidate split and verifies EP universality + u/d ratios.
- Bridge precision audit: does the split preserve (or derive) the 0.068% match?

**Documentation per variant**: Same pattern as G9 — setup script + findings file with explicit equations, outcome, and rationale for continuation or abandonment.

**Exit criteria for the G1 track**:
- A concrete, representation-theoretically natural split + chain that accommodates both the bridge scale and the 3-gen spectrum at the observed precision, or
- A strong negative result (e.g., "no chain of subalgebras of \(\mathfrak{so}(8)\) leaves a 3-dimensional light sector while preserving the 784 count") with a clear recommendation (reinterpret the bridge, change the algebra, accept two scales, etc.).

---

## Gated: Dynamical Lagrangian Construction (`v60/lagrangian/`)

**Gate condition**: At least one live candidate from G9 *or* G1 (ideally both) that survives the initial LIGO polarization and rank-tension numerical filters.

**Target specification** (living document `TARGET_SPEC.md`):
- The action \(\mathcal{L}\) must be written on the algebra (or the appropriate bundle after the G9 recast).
- Its EL equations must admit the static algebraic \(\Omega(x)\) structures (Brannen kernels on the grade-selected sectors) as exact or long-lived approximate solutions.
- Linearization around the vacua must reproduce the observed gauge content, the Brannen mass spectra, and (post-G9) a viable tensor gravity sector.
- All structural integers (5, 21/16, 2/9, 28, 784, …) must appear as derived consequences or as the unique values that satisfy the EL equations + boundary conditions.
- The theory must pass the existing v59 precision tests (0.02–0.07% bridges) and the new LIGO / EP filters.

**Early (pre-gate) work**:
- Port the v59 effective Lagrangian into `v60/lagrangian/v59_target.md` with explicit "derive this" annotations for every structural feature.
- Sketch the minimal kinetic + potential terms that could be consistent with a 2-form carrier (once G9 has a candidate).
- Identify which pieces of the existing Lean corpus (e.g., the Brannen kernel eigenvalues, the constraint surface) can be re-used as "on-shell conditions" derived from the EL equations.

**Post-gate variants** (to be expanded once the gate opens):
- Variant L1: Minimal 2-form + scalar + gauge Lagrangian whose EL equation is the curved-space generalization of \(\square \Omega = f_g \rho_{\rm grav}\).
- Variant L2: Full multivector \(\Phi(x) \in \mathrm{Cl}(7)_{\rm even}\) with potential engineered from the grade projectors and the Z₂×Z₂ selection (once G1 clarifies the selection mechanism).
- Variant L3: First-order / BF-style formulation that makes the soldering constraint part of the dynamics.

**Lean target (post-gate)**: `LagrangianEL.lean` — statement of the action and derivation of (at least the quadratic + free) EL equations; proof that the algebraic solutions satisfy them.

---

## Synthesis, Consistency, and Cross-Checks (`v60/synthesis/`)

- Reusable test harnesses: polarization / helicity, gravity charge / EP, bridge precision, Koide reproduction, scale factors.
- End-to-end script: take a candidate (induced metric + two-piece Y), compute the implied Lagrangian, verify it reproduces the v59 algebraic structures + new geometric requirements.
- All runs archived as timestamped JSON + figures (small files stay in the repo; large data to `/space/scp/`).
- Figures and summary tables for the eventual CLOSEOUT.

---

## Lean Engineering Notes

- New files under `v60/lean/` (or as a subdirectory that can be added to the v59 lakefile if desired).
- Every new headline theorem must eventually be proved or explicitly `sorry`-flagged with a clear conjecture statement.
- Run `lake build` + `AxiomCheck` after every significant addition (or at least at each session close).
- Import strategy: prefer importing the proved theorems from the v59 corpus rather than copying; only duplicate when a v60-specific re-statement is required for clarity.

---

## Documentation Rhythm

- Raw daily notes: `v60/notes/2026-05-25-*.md` (or per agent).
- After each variant completion or major insight: `0X_findings.md` in the appropriate subdir + entry in `v60/synthesis/FINDINGS_*.md`.
- Weekly (or per major milestone): update `ROADMAP.md` and this `PLAN.md` with status, demotions, and new priorities.
- At version close (or major checkpoint): full `CLOSEOUT.md` + delta to `UNIFIED_THEORY.md`.

---

## Resource and Scope Guardrails

- **Time box on search**: If after 3–4 variants per track (or 2–3 weeks of focused effort) neither G9 nor G1 has a live candidate, escalate to a "restructuring review" rather than indefinite continuation.
- **No kernel modification**: The unified simulation kernel remains off-limits unless the user explicitly authorizes changes for a specific dynamical test.
- **SFA only**: Any future simulation output (e.g., to test the new Lagrangian on a lattice) must use the SFA format and the existing `scp-runner` tooling.
- **Honesty**: A strong negative result on G9 ("no natural induced metric yields 2 TT DOF within the algebra") is a valid and valuable v60 outcome. It does not constitute failure of the project; it constitutes progress in mapping the space of viable theories.

---

**This plan is a living document. It will be revised in-session as the first results arrive. The current authoritative snapshot is always this file + the per-variant findings.**