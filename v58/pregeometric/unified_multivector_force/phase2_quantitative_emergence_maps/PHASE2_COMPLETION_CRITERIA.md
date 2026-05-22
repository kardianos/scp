# Phase 2 Completion Criteria — Quantitative Emergence Maps

**Phase 2 is considered complete only when ALL of the following conditions are simultaneously satisfied:**

## 1. Quantitative Emergence Map Implemented
- At least one concrete, quantitative "emergence map" has been implemented in Python that takes a background-free relational graph (produced under the living candidate) and outputs derived geometric quantities.
  Examples of acceptable first maps (any one is sufficient to start):
  - Effective retarded distance between regions (shortest causal path length, or weighted by edge strength / integrated |Ω|).
  - Local volume element or growth-based metric proxy (ball size vs. depth, local `d_eff` fit, curvature proxy via second differences).
  - Effective light-cone / isotropy measure derived from branching or null-direction statistics.

- The map must be applied to ensembles of graphs with ≥ 500 nodes from Phase 1 (or newly generated equivalent relational graphs under the exact living candidate).

## 2. Error Quantification Against Phase 1 Criteria
- The outputs of the map have been measured for error / consistency against the §3 criteria in `EMERGENCE_OF_EFFECTIVE_3D_GEOMETRY.md` (local isotropy, effective dimensionality in the 3 ± δ range, bounded effective curvature, stability under compensating dynamics).
- Quantitative error numbers or bounds are produced and logged (not just qualitative statements).

## 3. Lean Certification on Real Data
- At least one non-trivial property of the emergence map (or of the quantities it produces) has been machine-checked in Lean on **real exported snapshots** from the Python implementation.
  Examples:
  - The map preserves monotonicity of causal balls.
  - Under the safe-band living candidate + protection fraction ≥ X, the extracted local `d_eff` lies in [2.7, 3.3] for all tested exported graphs of size ≥ N.
  - Some invariance or error-bound theorem holds on concrete exported configurations.

## 4. Python ↔ Lean Alternation Demonstrated
- At least **two full alternation cycles** of the form:
  Python (implement / refine map + apply to graphs + export map outputs) → Lean (ingest + prove property of the map) → feedback into Python (refinement or new map variant).

## 5. Logging of Any Phase 1 Adjustments
- If any code or data in the Phase 1 folder (`../phase1_minimal_relational_models/`) is modified to support Phase 2 work, all changes must be logged in a file called `PHASE2_CHANGES_TO_PHASE1.md` (or equivalent) inside this Phase 2 folder, with clear descriptions of what was changed and why.

---

**Notes for the agent:**
- The goal of early Phase 2 is to produce the *first* quantitative, certifiable bridge from the relational field configurations to effective geometric quantities, not to produce a perfect or final map.
- Reuse as much Phase 1 infrastructure as possible (`RelationalNode`, export format, living-candidate evolution, certified extractors).
- Keep all work strictly background-free.
- The living candidate must continue to be used exactly (no parameter changes unless explicitly logged and justified).

This file is the authoritative definition of "Phase 2 done" for this sub-folder. Do not declare completion until every condition above is met with evidence.