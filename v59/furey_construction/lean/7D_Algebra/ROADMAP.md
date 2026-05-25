# Roadmap — 7D Algebra Lean Realization (Cl(7)_even / Spinor 8 Action)

**Status**: Living. Updated after every phase or major discovery.
**Current overall progress**: Phase 0 COMPLETE. Phase 1 Angle 1 (table + gammas + exact Fock labeling + first L/F 8×8 matrices) COMPLETE (see notes/2026-05-23-Phase1-Angle1-*.md). `lean SevenDAlgebra.lean` type-checks cleanly. Concrete matrices delivered. 60% of Phase 1 minimum achieved in first execution pass.
**Parent docs**: `PLAN.md` (detailed tactics), `README.md` (orientation + references), `../README.md`
**Overall Goal**: Deliver a working, buildable Lean module (or small set) that supplies explicit, computable 8×8 matrices for the action of L-grade and F-grade generators of the 7D algebra on the Fock basis |Ω_N⟩ (N=0,1,2,3). This removes the last technical obstruction identified by the Furey Fock-space agent and multiple synthesis notes.

---

## Milestone Phases

### Phase 0: Setup & Documentation (COMPLETE — this session start)
**Dependencies**: None (pure creation).
**Deliverables**:
- `v59/furey_construction/lean/7D_Algebra/` directory + `notes/` subdir.
- `README.md`, `PLAN.md`, `ROADMAP.md` (initial versions) referencing all required Python/Lean/Maxima files.
- First dated note `notes/2026-05-23-initial-setup-and-plan.md`.
- Decision recorded: location under furey/lean/ (justified in PLAN).
**Success Criteria**: All three docs exist, are non-empty, and the justification for folder choice is explicit. `lake build` from parent still succeeds (no breakage yet).
**Owner**: This agent.
**Estimated effort**: 1-2 hours (high-quality writing + cross-reference verification via tools).

### Phase 1: Core 7D Algebra + Explicit 8×8 Spinor Action Matrices (PRIMARY — aggressive push this session)
**Dependencies**: Phase 0 docs; access to `octMultTable` (from density/OctonionAlgebra.lean and Maxima), `fano_triples` (from Octonions.lean), the Fock construction in 02/07/13 Python (for exact labeling).
**Sub-milestones** (must all be attacked; use 7 angles/strategies where possible per the spirit of prior Lean agents):
1. **Fano / Multiplication Table internalization**
   - Port exact table (8×8 of (sgn, target) pairs) into Lean.
   - Prove basic theorems: 7 lines, 21 incidences, reproduction of the triples from Octonions.lean.
   - Define `octLeftMult` and `gamma k` (Fin7 → 8×8 Int matrix).
2. **Spinor / Fock Basis Labeling**
   - Exact replica of Python `basis_subsets = [combinations(range(3), k) for k in 0..3]`, `N_values`, `lepton_indices = [0,7]` (assuming ordering ∅, singles, doubles, full), `dQuarkIndices`, `uQuarkIndices`.
   - `def OmegaN (n : Nat) (color : Fin 3) : Fin 8` or indexed accessors.
   - Ordering comment explaining any permutation vs Python.
3. **Action Implementation**
   - Matrix multiplication helpers (or use Mathlib.Matrix if imported cleanly).
   - Compute at least one concrete L-grade operator matrix (e.g. bivector 0-1 or first Fano line product gamma0 * gamma1) w.r.t. the Fock-ordered basis.
   - Compute at least one F-grade proxy (e.g. a 4-form realized as product of 4 gammas, or the *φ singlet identified via G2 invariants).
   - `#eval` the full 8×8 (or sparse nonzero entries) and the 2×2 lepton subblocks, 3×3 d-quark blocks etc.
4. **Basic Verification**
   - Confirm Clifford relations on the gammas (at least squares and a few anticommutators) via `#eval` or decidable `Prop`.
   - Show that the L-bivector action on the 3D Fock sub-system already produces the universal diagonal-zero on pure rotations (replicate 13.py Part 3).
   - Document the structural meaning of the extra 4 octonion directions (they supply the F-grade content that only quarks see).
5. **Documentation & Notes**
   - At least 3 dated progress notes (one after table+labeling, one after first matrices, one after verification).
   - Update PLAN.md + this ROADMAP with "Phase 1 status" + concrete "next 3 tasks".
**Success Criteria (minimum for "Phase 1 done")**:
- Single primary `.lean` file (e.g. `7D_Algebra.lean` or split `FanoTable.lean` + `SpinorRep.lean`) that `lake build` accepts with no errors.
- Fock labeling + sector lists are **identical** (modulo documented permutation) to the Python in 13_fock_mass_forcing.py.
- At least two distinct 8×8 action matrices (one L, one F) are printed by `#eval` and recorded in a note.
- Clear comment block: "How a future FockMassForcing.lean imports these to state the Z₂×Z₂ theorem".
- No scope creep into full 64D regular rep or C⊗H⊗O.
**Risk / Angle strategy**: 7+ approaches tried in parallel where blocked (e.g. pure matrix vs inductive CliffordAlgebra, Fin 8 vs custom FockState type, Int vs Rat coeffs, hard-coded 8 states vs generating function, etc.). Record every angle in notes even if abandoned.
**Owner**: This agent (high-effort pass).
**Exit to Phase 2**: When the above minimum is met and a summary "what compiles / what still needs work / concrete next steps" is written.

### Phase 2: Grade Projectors, More Generators, μ-Bisection (Next Dedicated Pass)
**Dependencies**: Phase 1 matrices + labeling; cl7_even.py grade logic; existing Lean dim facts (L=28, F=35).
**Milestones**:
- Full (or schematic) 64-basis even multivector type + `grade : Cl7Even → Nat` and `P_L`, `P_F` as algebra-level projectors (port from cl7_even.py).
- Map from abstract L/F basis elements to their 8×8 images under the spinor rep (at least a generating subset of 28+35).
- Concrete *φ (the G₂ singlet 4-form in F) realized as a specific linear combination of gamma products; prove it commutes with the G2 action (link to existing G2 facts in SpinDimension.lean).
- First machine-checked statements:
  - "Bivectors from the color-3 directions (L content used by leptons) act diagonally or with selection rules compatible with singlet N=0,3."
  - "4-forms carrying the *φ component (F) have non-trivial action precisely on the colored sectors."
- Update AxiomCheck-style list.
**Success Criteria**: At least 5-10 additional generators have explicit matrices; one theorem (even if with `sorry` on a geometric input) linking the observed bit pattern to the G2 content difference between L and F.
**Notes**: 7+ dated notes expected (one per major angle/strategy).

### Phase 3: Computability, Forcing Theorems, Cross-Validation
**Dependencies**: Phase 2 + Python numeric matrices from 13.py + brannen_kernel.py embeddings.
**Milestones**:
- Decidable predicates `DiagonalOnSector (M : Matrix 8 8) (sector : List Fin8) : Bool`.
- `#eval` that reproduces (within rounding) the numeric diagonals from 13_fock_mass_forcing.py Part 3 for the sample L-op.
- JSON exporter (simple `IO.println` of arrays) so that `lean --run` or lake script produces a `.json` consumable by Python for bit-exact comparison.
- Theorem (or high-confidence axiom + numeric cert) "the only bit pattern (L for leptons, F for d, both for u) that simultaneously (a) gives non-zero mass diagonal, (b) preserves color irrep of the N sector, (c) matches the Brannen |ξ|² = 1-14/D values is the observed Z₂×Z₂".
- Integration point: the new matrices appear in an extended `Predictions.lean` or new `FockMassForcing.lean` inside the parent lean/.
**Success Criteria**: End-to-end "from abstract grade → concrete 8×8 matrix → verified selection rule on |Ω_N⟩" for at least the lepton and d-quark cases. At least one fully `sorry`-free structural theorem that future agents can cite.

### Phase 4: Polish, Hand-off, Project Integration
**Dependencies**: Phase 3 + review by synthesis / Furey agent.
**Milestones**:
- All code heavily commented with direct pointers back to the exact line numbers in 13_fock_mass_forcing.py, cl7_even.py, 02_sm_idempotent.py etc.
- `notes/` contains a "final status + hand-off list" (what is buildable today, what the next 5 concrete edits are, which axioms are the highest priority to discharge).
- Parent `furey_construction/lean/README.md` updated with a new row in the table: "7D_Algebra.lean — explicit Cl(7)_even action on Fock 8 spinors; matrices for L/F generators".
- Parent PLAN.md / 13_fock report updated: "the explicit matrix rep obstruction has been removed; the rep-theoretic forcing is now formalizable".
- Optional: a small `test/` or `examples/` showing usage (e.g. "compute action of the lepton ξ bivector on |Ω_0⟩ and |Ω_3⟩").
**Success Criteria**: The folder is "future-agent ready" — a new Lean specialist can `lake build`, read the top comments, and immediately start writing the vanishing or consistency theorems without re-exploring the Python.

---

## Cross-Cutting Concerns (All Phases)
- **Build hygiene**: Every edit must leave `cd v59/furey_construction/lean && lake build` succeeding (or a documented temporary breakage with recovery steps). Prefer adding to the existing package over a second lakefile.
- **Dated notes discipline**: After every 1-2 hours of work (or after each sub-milestone), write `notes/YYYY-MM-DD-HHMM-PhaseX-AngleY.md` or similar. Even "Angle 3 failed because ..." counts. This mirrors the successful 7-angle attack pattern from the density Lean roadmap agents.
- **Living docs**: After every note, do a search-replace on PLAN.md and ROADMAP.md to mark completed sub-steps and add "learned X, therefore next is Y".
- **Citation & traceability**: Every definition that comes from a Python file must have a comment `-- port of cl7_even.py:157 (P_L) and 13_fock...py:73 (lepton_indices)`.
- **Lightweight vs Mathlib**: Start with `import Mathlib.Data.Fin.Basic` + manual matrix (Array or `Matrix (Fin 8) (Fin 8) Int` from Mathlib if the import is cheap). Only pull in heavy linear algebra when the simple version is insufficient for proofs.
- **Signature / convention tracking**: Maintain a single source-of-truth comment block "Project convention: Cl(7,0) with e_i² = +1 (per cl7_even.py). Octonion left-mult gives squares -1 → we will use i*gamma or document the Wick rotation / dual picture for the spinor rep."

## Overall Success Metric for the Whole Effort
When a future agent (or the original Furey Fock agent) can open `7D_Algebra.lean`, write 20 lines that say

```lean
import V59.Furey.7D_Algebra.SpinorRep
-- ...
theorem lepton_skips_F_for_mass (op : FGenerator) :
    ∀ i j ∈ leptonIndices, (actionMatrix op)[i, j] = 0 ∨ (inconsistentWithBrannenEmbedding i j op) := by
  -- now provable because we have the matrices
  sorry  -- or real proof
```

and `lake build` + the theorem statement type-checks, the mission is complete.

The current session (this subagent) is responsible for **Phase 0 + substantial, documented progress on Phase 1** (at minimum the table, labeling, first two matrices, and the hand-off summary).

---

*Roadmap written at the opening of the dedicated high-effort session. Will be updated in place.*