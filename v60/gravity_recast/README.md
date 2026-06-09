# v60 Gravity Recast — G9 Induced / Emergent Metric Attack

**Priority**: #1 (decisive blocker per v59 CLOSEOUT and gaps/SYNTHESIS).

**Goal**: Recast the v59 scalar connection \(\Omega \in \Lambda^2\) (or its OBE dynamics) as a fundamental 2-form whose soldering or simplicity constraint to the Lorentz group derives a symmetric spacetime tensor carrying exactly two transverse-traceless degrees of freedom, while preserving the second-moment gravity charge and the long-range force law.

**Primary literature / code anchors**:
- `../v59/gaps/gravity/g9_polarization_test.py` and `g9_spin2_route_test.py` (helicity decomposition baseline).
- `../v59/gaps/gravity/G8G9_Gravity.lean` (Lean statements of the scalar nature and the 21/16 prefactor).
- Plebański 2-form gravity (self-dual 2-forms + simplicity constraints) — for inspiration and technical toolkit.

**Subdirectory layout** (to be populated):
- `01_plebanski_setup.py` + `01_findings.md` (first ansatz)
- `0X_<variant>_*.py` + findings for each major family (A = direct internal 2-form, B = full-algebra valued, C = auxiliary soldering field, D = defensive no-go scans).
- `Lean/` or pointers to new modules under `../lean/` (e.g., `InducedMetric.lean`).
- `figures/`, `data/` for archived runs.

**Success / exit criteria**: See `../PLAN.md` (G9 Track) and `../ROADMAP.md`.

**Guardrails**:
- Bounded search (do not chase indefinitely).
- Every variant gets a documented outcome (success, clean failure with reason, or "abandoned after X because...").
- All candidates must pass the polarization test harness (2 TT DOF + correct coupling to \(\rho_{\rm grav}\)).
- No modification of the unified simulation kernel without explicit authorization.

**Status (2026-05-25)** — authoritative artifacts are `04_*` and `../lean/G9Soldering.lean`:
- `01_8space_to_spacetime_bindings.py` + `01_findings.md`: binding *enumeration* (useful context).
- ~~`02_constrained_helicity_count.py`~~: **CIRCULAR — superseded.** It inserted the
  spin-2 generator by hand (`Jz_on_sym2_transverse ⊗ Id`) and reported back ±2; the
  "constraint" did no work. The Lean it cited did not compile.
- **`04_soldering_helicity_honest.py`**: the honest, non-tautological computation.
  One machinery, four scenarios; the only thing that reaches ±2 is **soldering**
  (the internal Lorentz index co-rotating). Runs clean; all assertions pass.
- **`05_dof_and_weakfield.py`**: exact massless-spin-2 DOF count = **2** (rank
  computation, helicities ±2) and the weak-field bridge showing v59's `□Ω=f_g ρ`
  is the trace sector of linearized GR. See `05_findings.md`.
- **`06_lorentz_commutant.py`**: the Schur no-go — internal `Spin(7)` commutant in
  `End(ℝ⁸)=Cl(7)_even` is exactly 1 ⇒ no `so(3,1)` fits inside ⇒ the carrier is
  *forced* into the spacetime factor. (C3 derivation.)
- **`07_unified_algebra.py`**: explicit `Cl(3,1) ⊗ Cl(7)_even` (v59's own stated
  factorization). Spacetime Lorentz = the 3φ+3θ Cosserat bivectors; commutes
  exactly with the full internal `Spin(7)` (incl. `SU(2)_L`, color, triality).
  **Closes C3.** See `07_findings.md`.
- **Lean** (`../lean/`, all build clean against v59 Mathlib, axiom-only):
  `G9Soldering` (soldering = Minkowski sum; no-go vs ±2; charpoly X²+4; DOF=2;
  Schur obstruction `no_internal_lorentz`), `G9Unification` (general tensor-factor
  commutation `spacetime_internal_commute` + dim identities), `G9InducedMetric`,
  `G9ToyHelicity` (both rewritten from the broken originals).

See `04_findings.md` → `05_findings.md` → `07_findings.md` for the full arc:
soldering mechanism → exact 2 DOF + weak-field bridge → C3 closed in
`Cl(3,1)⊗Cl(7)_even`. The remaining open item is the full Plebański action + EOM.

See the G9 section of `../ROADMAP.md` for the detailed open questions (Q-G9-1 through Q-G9-5) and session priorities.