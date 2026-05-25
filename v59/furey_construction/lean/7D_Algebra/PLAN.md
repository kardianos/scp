# Implementation Plan — 7D_Algebra Lean Realization (Cl(7)_even Spinor Action)

**Date**: 2026-05-23 (Phase 1 — living document; will be updated after each milestone)
**Latest update**: 2026-05-23 after Angle 1 execution — Phase 1 core (7D generators + matrices on |Ω_N⟩) delivered and verified (lean checks, concrete numbers obtained via mirror). See notes/ for details. PLAN remains the contract; deviations documented.
**Folder**: `v59/furey_construction/lean/7D_Algebra/`
**Owner**: This subagent (and any follow-on Lean specialists)
**Related**: `README.md` (this dir), `../README.md`, `../../13_fock_mass_forcing.py`, `../../../algebra/cl7_even.py`, `../../../density_algebra/lean/OctonionAlgebra.lean`

## 1. Strategic Choices & Justifications

### Location Decision
Chose `v59/furey_construction/lean/7D_Algebra/` (subdirectory of the existing Furey Lean package) rather than a new top-level `v59/lean_7d_algebra/`.

**Rationale**:
- The work is **Furey-construction-specific**: it exists to close the exact gap called out in the Furey Fock-space agent's report (13_fock_mass_forcing.py:382 "explicit Cl(7)→Cl(6) grade correspondence in the matrix rep") and the synthesis notes.
- Co-location with `../Octonions.lean`, `../SpinDimension.lean` etc. allows direct reuse/extension (fano_triples, dimSO7=21, L/F mentions in Predictions) and a single `lake build` for the whole furey/lean/ tree.
- The Python drivers (02/07/13, cl7_even.py, brannen_kernel.py) all live under `v59/furey_construction/`, so the Lean artifacts (exported matrices, theorems) can be cross-checked against the same directory's .json / .npz outputs.
- A central `v59/lean_7d_algebra/` would duplicate lake manifests/toolchains and split the "Furey algebraic physics" concern from the "density protection stability" concern (the latter already has its own `v59/density_algebra/lean/` with separate OctonionAlgebra + 4-phase roadmap).
- Future integration (e.g. using the new matrices inside a BrannenKernel or Predictions extension) is simpler when everything Furey-Lean is under one Lake package root.

If the parent lakefile becomes unwieldy, a future refactor can turn 7D_Algebra into its own Lake package with `require` on the parent libs; for now, we add it as an additional `lean_lib` entry.

### Scope Boundaries (Strict — Do Not Creep)
- **In**: 7D vector algebra (Fano/octonion table or full Cl(7) generators), even subalgebra (or at least its action via products), 8D spinor module (Fock/Witt or octonion regular rep), explicit matrix representations (at least for a generating set of L-grade bivectors and F-grade 4-forms / *φ), grade projectors (L vs F), action on the labeled |Ω_N⟩ basis (N=0..3 with color subsets), verification of diagonality or selection rules for a few concrete operators, documentation of how it plugs into Z₂×Z₂ forcing + Hessian.
- **Out (for this dedicated effort)**: Full C⊗H⊗O 64D tensor product, full SM gauge embedding proofs, numerical mass spectra, gravity sector, new physics predictions. Those stay in Python or other Lean modules. Also out: porting the entire 64D regular representation matrices (too heavy for first pass; the 8D irrep is the target).
- **Pragmatism rule**: Use `axiom` / `sorry` / `noncomputable` only when the mathematics is known correct but the Lean proof is long (e.g., openness of G₂ orbit, full associativity of Clifford product for all 64×64). Document every one with a precise "TODO: replace by proof" + reference to the Python/Maxima that already checked it numerically. Prefer `#eval` and decidable instances for concrete matrices.

### Naming & Terminology (Project Convention)
- Use `OmegaN`, `|Ω_N⟩`, `N0_lepton`, `N1_dquark`, `N2_uquark`, `N3_lepton` for the sectors.
- `L_content`, `F_content`, `L_mask`, `F_mask`, `mu_projector` (following cl7_even.py and brannen_kernel.py).
- `FanoTriple`, `octMultTable` (match density/OctonionAlgebra.lean).
- `Cl7EvenBasis`, `Grade`, `CrossingSign`.
- "braid" only for z-aligned sub-component (quark analog) if it appears; here mostly "Fock state" or "spinor basis vector".
- "particle" for the physical interpretation of the |Ω_N⟩.

## 2. High-Level Architecture

### Layer 1 — Algebra Core (Cl7Even or 7D GA)
- Either:
  (A) Full even-subalgebra: type `Cl7Even` ≅ `Fin 64 → ℝ` (or `Vector ℝ 64`), basis indexed by even-cardinality subsets of `Fin 7` (exactly as `BASIS` + `BASIS_IDX` + `GRADES` in cl7_even.py).
  Multiplication via a precomputed table or a function `mul (S T : EvenSubset) : (Sign × EvenSubset)` using the `crossing_sign` / bubble-sort logic ported from Python.
- Or (B) (lighter start) the 7D imaginary algebra only: `ImagOct := Fin 7`, with multiplication table from the 7 Fano triples (extend `fano_triples` in Octonions.lean). This gives the generators; the even algebra is then the algebra they generate inside End(Spinor).
- **Decision for Phase 1/2**: Start with (B) + explicit 8×8 matrices for the 7 left-multiplication operators (using the exact `octMultTable` from density/OctonionAlgebra.lean). This immediately gives "7-dimensional algebra realization" + computable action. Later lift to full multivector if needed for grade projectors inside the algebra itself.
- Grade decomposition can be axiomatized at first ("L generated by bivectors + 6-vectors that land in the μ=-1 eigenspace") and later made computable by tagging basis elements.

### Layer 2 — Spinor Module (8D Irrep)
- `SpinorBasis := Fin 8` (or better: a structure `FockState` carrying the subset `List (Fin 3)` or bitmask `Fin 8` + `N : Nat` computed as popcount, to match Python `basis_subsets`, `N_values`, `subset_to_idx` exactly).
- Two dual views (both needed):
  1. Octonion regular view: identify basis 0..7 with 1, e1..e7; action of the 7 generators = left multiplication matrices extracted from octMultTable. This gives faithful 7D algebra action.
  2. Fock/Witt view: the exterior algebra ∧* ℂ³ (or complex 3D Fock space of Cl(6)). The 8 states are exactly the same vector space; the identification map `oct_basis ↔ fock_subset` can be chosen (or left as a `def equiv` with axioms for now) so that the color weights (SU(3) action on the 3 indices) line up with the G₂ action on the octonion imaginaries.
- The target output: for any algebra element g (e.g. a specific bivector e_i e_j or a 4-form), the 8×8 matrix `actionMatrix g : Matrix (Fin 8) (Fin 8) ℂ` (or ℝ) w.r.t. the ordered Fock basis. Columns = coefficients of g · |Ω_j⟩ expressed in the |Ω_i⟩ basis.

### Layer 3 — Grade Mapping & L/F Selectors
- Define `L_generators : List (Matrix ...)` (the 28 basis elements of Λ²⊕Λ⁶ mapped to 8×8 operators).
- `F_generators : List (Matrix ...)` (35 of Λ⁴).
- Concrete sample operators:
  - L-sample: the three so(3) generators on the first three directions (bivectors 01,02,12) — these should correspond to the "gauge" content used by leptons.
  - F-sample: the G₂-invariant 4-form *φ (coassociative) and a few others that carry the color structure for quarks.
- Projectors: `P_L`, `P_F` as 8×8 (or on the algebra side first).
- The key theorem target: "for lepton_indices, the matrix elements of any pure-F generator between N=0 or N=3 states satisfy <i| F_op |j> = 0 when i,j both lepton (or the diagonal mass term vanishes / is inconsistent with Brannen embedding)" — or the rep-theoretic version that only the observed bit choice preserves both the irrep and the |ξ|² constraint.

### Layer 4 — Integration Hooks
- Export functions or `#eval` that print the matrices in a format importable by Python (JSON of arrays) for cross-validation against 13_fock_mass_forcing.py numeric output.
- Theorems that the existing Lean facts (dim 21 for L bivectors, 35 for F, 14/21=2/3) are compatible with the dimension of the image of the representation map.
- Later: feed into `BrannenKernel.lean` extension ("the ξ embedding for leptons lives in the image of L_generators under the rep").

## 3. Detailed Implementation Phases (Execution Order)

### Phase 0 — Infrastructure (this session, immediate)
- Create folder + PLAN.md + ROADMAP.md + README.md (done).
- Add the new lib to the parent `lakefile.lean` (or create a minimal self-contained one for the subdir).
- Create `notes/2026-05-23-*.md` with this plan.
- Verify that `lake build` still works from parent after adding stub.

### Phase 1 — Core 7D Algebra + 8×8 Action (primary goal of this pass)
1. Port/copy the exact `octMultTable : List (List (Int × Nat))` (8 rows) from `density_algebra/lean/OctonionAlgebra.lean` into a new `FanoTable.lean` or inside the main file. Prove `length = 8`, each row length 8, and the 7 Fano incidences.
2. Define `Imag7 := Fin 7`, `OctBasis := Fin 8` (0 = 1, 1..7 = e_i).
3. `def octLeftMult (k : Fin 7) (m : OctBasis) : (Int × OctBasis)` using table lookup (or direct function).
4. `def gamma (k : Fin 7) : Matrix (Fin 8) (Fin 8) Int` — the 8×8 matrix whose (p, m) entry is the coeff from octLeftMult k m.
5. Prove (or `axiom` + test `#eval`) the Clifford relations for a generating set: `gamma k * gamma k = -1` (or +1 per project convention), `gamma i * gamma j + gamma j * gamma i = 0` for i≠j in Fano triples.
6. Define the Fock labeling exactly as in 13_fock_mass_forcing.py:
   ```lean
   def basisSubsets : List (List (Fin 3)) := ... -- or better inductive / Fin 8 -> Subset
   def NValue (i : Fin 8) : Nat := ...
   def leptonIndices : List (Fin 8) := [0,7]  -- N=0 and N=3 (the full set)
   def dQuarkIndices, uQuarkIndices : ...
   ```
   (Use `List.range 3`, `combinations` ported or hardcoded 8 entries for determinism.)
7. Implement matrix-vector action and matrix-matrix for the gammas.
8. For a concrete L-grade operator, e.g. bivector generator corresponding to rotation in plane 0-1 (or the first Fano line), compute `M := gamma 0 * gamma 1` (adjusted for pure bivector), print via `#eval` the 8×8 matrix in the Fock ordering, and the 2×2 subblock on lepton indices.
9. Repeat for one or two F-grade proxies (e.g. triple products or the volume element projected).
10. Document in a note: "what we see for diagonals", "how the 3D Fock bivectors already give the universal vanishing on pure rotations", "what the 7D adds (the extra 4 directions allow *φ singlet in F that commutes with G2)".

**Success for Phase 1**: A single `.lean` file that `lake build`s cleanly, `#eval`s at least one 8×8 action matrix for an L-generator and one for F, and has the Fock labeling + sector indices as defs. At least 3 dated notes.

### Phase 2 — Grade Projectors + More Generators (next session or continuation)
- Implement (or axiom) the full 64-basis even multivector type + mul (port crossing_sign logic; this is pure arithmetic, fully decidable).
- Define embedding of the 7 gammas into the even algebra (as products of two vectors).
- Define the 28 L and 35 F basis operators as explicit lists of 8×8 matrices (generated by taking products of the gamma's that land in the correct grade under the iso).
- `def isDiagonalOnLeptons (M : Matrix ...) : Prop` (computable predicate on the submatrix).
- First theorems: "all pure so(3) bivectors from the 3 color directions act with zero diagonal on the lepton singlet states" (replicates the 3D numeric result formally).
- Axiom or proof that the G2 singlet 4-form *φ (unique in F) has non-zero action on the color singlets but is "skipped" for consistency with the Brannen ξ placement (which is hardcoded in L bivectors per brannen_kernel.py).

### Phase 3 — Forcing Argument & Cross-Checks
- Encode the rep-theoretic forcing as a theorem (using the G2 branching facts already in SpinDimension.lean + new content of L vs F under G2).
- Export the matrices to JSON (via a simple printer) and verify byte-for-byte or numerically against a new Python driver that calls the Lean? (or just manual diff).
- Add to `Predictions.lean` or new `FockMassForcing.lean` the statements "the observed Z₂×Z₂ bit pattern is the unique one compatible with N-color irrep + G2 content of grades + Brannen D constraint".
- Remove `sorry`s where possible; replace axioms with proofs using Mathlib's `Matrix` and `Finset` tactics.

### Phase 4 — Polish & Hand-off
- Full documentation in all files.
- `AxiomCheck.lean` analogue listing what is still axiomatic.
- Update parent furey/lean/README.md and the main PLAN.md / 13_fock... report.
- Produce a "matrices for the 8 states" summary table in a note.
- Hand-off list for the next agent: "take the gamma matrices, implement the iso map from abstract Cl7 grades to these concrete 8x8, prove the selection rule as `theorem lepton_skips_F : ∀ op ∈ F_generators, diagonal_on lepton_indices op = 0` (or the precise version)".

## 4. Risks & Mitigations
- **Build bloat / Mathlib compile time**: Keep the new lib minimal; import only `Mathlib.Data.Matrix`, `Mathlib.Algebra.CliffordAlgebra.Basic` (if it exists and is lightweight) or avoid and roll our own `Matrix` ops with `Array`/`Fin` for `#eval` speed. Use `Float` or `Int` first, `ℚ`/`ℂ` later.
- **Signature mismatch ( +1 vs -1 squares)**: The project Python uses +1 for Cl(7,0) vectors; oct table gives -1 for left mult. Mitigate by defining `gamma k := i * L_{e_k}` or documenting the convention difference and using the project Fano for structure only. The combinatorial grade counts and L/F dims are signature-independent.
- **Choosing the "right" 7 directions / iso to the 3 Witt**: The 7 octonion directions contain the 3 color + 4 "broken" G2 directions. We can start with the first 3 for the Fock α's (matching 02/07), and treat the extra 4 as the ones supplying the F-grade *φ. Exact map can be axiomatized initially.
- **Non-associativity**: Octonions non-assoc, but the Clifford product **is** associative. The left-mult matrices **are** associative (matrix mult). Use the matrix representation as the source of truth for the action; the abstract table is for the algebra structure constants only.
- **Proof length for full 64D mul**: Enormous. Do not attempt in one pass; the 8D rep + a few generators is sufficient for the "explicit matrices for the action of generators on the 8 states".

## 5. Tooling & Verification
- Every new `.lean` must appear in the lakefile roots.
- After every significant chunk: run `lake build`, capture output, record in a dated note whether it succeeded and what `#eval` printed.
- Cross-check key numeric matrices by adding a small Python snippet that hardcodes the same Fano table and reproduces the Lean `#eval` output (store in `notes/`).
- Use `lean --version` and the exact v4.29.0 + mathlib v4.29.0 to match parent.

## 6. Success Criteria (Minimum for "Complete Phase 1")
- Dedicated folder with the three required docs + at least one dated note.
- One primary `.lean` module that builds cleanly.
- Fock basis labeling + sector indices as computable defs matching Python 13 exactly (up to ordering).
- At least the 7 gamma matrices defined and one derived L-bivector and one F-proxy operator with their 8×8 action matrices `#eval`'able.
- Clear comment / example showing how a future `FockMassForcing.lean` would import these matrices and state the Z₂×Z₂ theorem.
- This PLAN.md and ROADMAP.md updated with "Phase 1 status: X% — concrete next 3 tasks".

**This plan is the contract.** Any deviation must be justified in an update note and the document itself.

---

*Written at the start of the dedicated high-effort Lean 4 + algebraic physics pass.*