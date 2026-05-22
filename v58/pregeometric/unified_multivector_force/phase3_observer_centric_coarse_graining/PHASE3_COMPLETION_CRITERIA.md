# Phase 3 Completion Criteria — Observer-Centric Coarse-Graining

**Phase 3 is considered complete only when ALL of the following conditions are simultaneously satisfied:**

## 1. Internal Observer Model Implemented
- At least one explicit "internal observer" model has been implemented in the Python relational framework.
- The observer must be realized as a **stable, high-density, protected-chirality lump** (or small cluster of such lumps) whose internal structure is built from the same multivector degrees of freedom as the rest of the graph.
- The observer must maintain internal "clocks" and "rulers" defined purely through its own causal interactions (e.g., sequences of protected bivector oscillations or repeated causal round-trips within the lump).

## 2. Observer Reconstructs Its Own Effective Geometry
- The observer must be able to reconstruct a local effective geometry (distances, local light-cone structure, local `d_eff`, or local metric proxies) **solely from the causal relations it can access internally**.
- This reconstruction must be performed without global knowledge of the full graph — only using the observer's local causal past and internal structure.

## 3. Comparison with Global Map
- The geometry reconstructed by the observer(s) must be quantitatively compared to the global emergence map produced in Phase 2 (or an updated version of it).
- Concrete error / agreement numbers must be produced (e.g., average deviation in local distances, agreement on local `d_eff`, isotropy consistency, etc.).
- The comparison must be performed on real evolved graphs (≥ 500 nodes) under the exact living candidate.

## 4. Lean Certification
- At least one non-trivial property of the observer reconstruction process or the agreement between observer geometry and global map must be machine-checked in Lean on real exported data.
- Examples of acceptable properties:
  - The observer's local distance function is consistent with (a monotonic transformation of) the global retarded distance.
  - Local light-cone isotropy as seen by the observer lies within a stated bound of the global isotropy.
  - Some invariance or error-bound holds for the observer-reconstructed quantities on the exported living-candidate data.

## 5. Python ↔ Lean Alternation
- At least **two full alternation cycles** (Python implementation / observer evolution / reconstruction + export → Lean ingestion / proof → feedback) must be completed and logged.

## 6. Logging of Changes to Previous Phases
- Any modifications to Phase 1 or Phase 2 code/artifacts must be logged in `PHASE3_CHANGES_TO_PREVIOUS_PHASES.md` with clear descriptions.

---

**Notes for the agent**
- The spirit of Phase 3 is to close the loop between "global map" and "what an internal observer actually experiences."
- The observer should be as minimal as possible while still demonstrating genuine internal reconstruction.
- All work must remain strictly background-free and use the living candidate exactly.
- Reuse and extend the existing Python and Lean infrastructure from Phases 1 and 2 wherever possible.

This file is the authoritative definition of when Phase 3 is complete. Do not declare completion until every condition above is met with evidence.