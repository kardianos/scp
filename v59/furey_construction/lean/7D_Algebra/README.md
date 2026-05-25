# 7D Algebra Realization for Furey Fock-Spinor Action (Cl(7)_even / Octonion GA)

**Location**: `v59/furey_construction/lean/7D_Algebra/`
**Parent**: `../README.md` (Furey Lean), `../../PLAN.md`, `../../13_fock_mass_forcing.py`, `../../../algebra/cl7_even.py`
**Date**: 2026-05-23 (initial setup)
**Status**: Phase 1 docs + early implementation in progress. Dedicated sub-project to deliver the missing explicit Lean matrices for generator action on |Ω_N⟩ 8D spinors.

## Purpose
Provide a buildable, documented Lean 4 formalization of the 7-dimensional geometric algebra (primarily the 64-dimensional even subalgebra of Cl(7,0) ≅ ℂ ⊗ 𝕆 ≅ Cl(6)) together with its irreducible 8-dimensional spinor representation. This supplies the **computable matrix representations** of L-grade (Λ²⊕Λ⁶, 28-dim) and F-grade (Λ⁴, 35-dim) operators with respect to the Fock basis |Ω_N⟩ (N=0,1,2,3 labeling the lepton/d/u sectors).

This directly enables the machine-checked part of the Z₂×Z₂ (L/F bit) forcing argument for mass-term selection rules that was identified as the critical missing piece in:
- `13_fock_mass_forcing.py` (and its report)
- The Furey Fock-space synthesis notes
- Prior Lean stability-bounds and density agents

Once the explicit action matrices exist (even for a few key generators), future agents can:
- Prove (or refute) diagonal vs off-diagonal patterns per sector
- Internalize the grade-mapped mass operators for the Hessian/Brannen kernel formalization
- Close the loop on "why leptons skip F, d-quarks skip L"

## Key References (must be consulted in all extensions)
- **Python Furey infrastructure**:
  - `v59/furey_construction/02_sm_idempotent.py`: Witt basis α_i / ᾱ_i for Cl(6) ≅ M₈(ℂ), Fock vacuum |0⟩, primitive idempotents, SM particle assignment to N=0..3 states.
  - `v59/furey_construction/07_full_generation.py`: explicit construction of the 8 right-handed states |Ω_N⟩ via α products on vacuum; charge formula Q = (2N-3)/3 etc.; generation of full SM content.
  - `v59/furey_construction/13_fock_mass_forcing.py`: the 8 basis_subsets via combinations(range(3)), N_values, lepton_indices (N∈{0,3}), d_quark (N=1), u (N=2); clifford_vector_action, bivector_action; sample L-op matrix; conclusion that "explicit 7D Cl(7) matrix rep + iso to Cl(6) + grade-mapped ops" is required for the strict vanishing theorem.
  - `v59/algebra/cl7_even.py`: authoritative concrete realization — 64 even-subset basis of Cl(7)_even, crossing_sign MULT_TABLE, grade projectors P_L / P_F / P_LF, μ-eigenspace bisection (L = μ=-1), left_mult_matrix on the 64D regular rep.
  - `v59/algebra/brannen_kernel.py`: explicit embedding of ξ (ℍ) into L-slice bivectors for leptons; *φ coassociative 4-form in F for quarks.
- **Maxima (Fano / octonion table source)**: `v59/density_algebra/maxima/octonion_sensitivity_analysis.mac` (and full_octonion_perturbation.mac) — supplies the exact 7-line Fano multiplication table used for all project octonion work (copied into `OctonionAlgebra.lean`).
- **Existing Lean**:
  - `v59/furey_construction/lean/Octonions.lean`: fano_triples, incidence count 21 = dim Spin(7).
  - `v59/furey_construction/lean/SpinDimension.lean`, `LieDimensions.lean`: dimG2=14, dimSO7=21, structural ratios.
  - `v59/density_algebra/lean/OctonionAlgebra.lean` + `StabilityBounds.lean`: the 8×8 octMultTable (from Maxima), L/F masks, effective Hessian models.
  - `v59/furey_construction/lean/README.md`: documents the current (pre-7D) state; explicitly flags "The Furey ℂ ⊗ ℍ ⊗ 𝕆 representation theory (Cl(6) Witt decomposition...)" as not yet verified in Lean.
- **Theory docs**: `v59/furey_construction/13_fock_mass_forcing_report.md`, `09_ckm_and_selection.py` (μ bisection), `v59/INTEGRATION.md`, `v59/SUMMARY.md`, CONCEPT.md (L/F single-source).

## Folder Layout (planned)
```
7D_Algebra/
├── README.md          # this file
├── PLAN.md            # living implementation plan (this session + future)
├── ROADMAP.md         # phased milestones + success criteria
├── 7D_Algebra.lean    # main module (or split: Cl7Even.lean, SpinorRep.lean, LFActions.lean, FockBasis.lean)
├── lakefile.lean      # optional: self-contained package or point to parent
├── lean-toolchain
├── notes/
│   ├── 2026-05-23-initial-setup.md
│   └── YYYY-MM-DD-*.md (dated progress after each chunk)
└── (future) data/     # exported matrices as JSON or Lean literals for cross-check with Python
```

## How This Plugs Into the Larger Project
- Supplies the "explicit iso Cl(7)_even → 8×8 matrices with grade basis" demanded by 13_fock_mass_forcing.py:379.
- Enables formal statement inside Lean of the rep-theoretic forcing (G2 content of L vs F + N→color irrep map selects the observed bit pattern).
- Feeds the Furey Fock agent and the stability-bounds Lean (internal 8×8 Hessian from octonion table can now use grade-projected operators).
- Once matrices for sample L-bivectors and F-4-forms are computable, we obtain machine-checked analogues of the numeric diagonals computed in 13.py Part 3/4.
- Long-term: axiom-free theorems "the only G2-singlet mass channel compatible with lepton singlets lives in L", "d-quark color weights require F", etc.

## Build / Usage (once implemented)
```bash
cd v59/furey_construction/lean/7D_Algebra
lake build
# or (preferred for integration): from parent lean/ dir, after updating lakefile.lean to include the new lib
lake build
# Then #eval or theorem statements that compute/print the 8×8 action matrices for key generators
# and prove e.g. "diagonal entries on lepton_indices for wrong-grade ops satisfy the selection rule"
```

All code must respect project conventions (SFA not relevant here; use existing naming L_content, F_content, Brannen kernel slices, μ projector, |Ω_N⟩ or OmegaN states, Fano triples).

This is the concrete engineering step that turns the high-level "representation-theoretic forcing" narrative into verifiable Lean objects.

---

*Initial setup by Grok Build subagent — Lean 4 + algebraic physics specialist. High-effort dedicated pass on the repeatedly flagged missing piece.*