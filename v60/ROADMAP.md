# v60 ROADMAP — G9 + G1 Attack and Dynamical Lagrangian

**Date**: 2026-05-25 (kickoff)
**Status**: Master open-questions list and session priorities for v60. Built directly on the v59 closeout audit.
**Parents**: [`../v59/CLOSEOUT.md`](../v59/CLOSEOUT.md), [`../v59/gaps/SYNTHESIS.md`](../v59/gaps/SYNTHESIS.md), [`../v59/ROADMAP.md`](../v59/ROADMAP.md) (legacy frontiers now subordinated)

This document supersedes the v59 four-frontier roadmap for the purpose of v60 work. The algebraic skeleton is trusted; the blocking problems are now **geometric/dynamical consistency** (G9) and **spectrum-bridge compatibility** (G1).

---

## Big Picture: the v60 mandate (from CLOSEOUT §IV)

> Only after G9 + G1 are resolved: write the dynamical Lagrangian \(\mathcal{L}_{v60}\) whose Euler–Lagrange equations yield the OBE \(\Omega(x)\) structure.

**Ranked gaps (unchanged order from gaps/SYNTHESIS + CLOSEOUT)**

1. **G9 — Tensor gravity mode (induced-metric recast)** — *Decisive blocker*
2. **G1 — Rank tension (two-piece Y / Goldstone chain)**
3. **G3 — α value (RG fixed-point or other mechanism)**
4. **G7 — Radian-insert (why Q re-enters as phase argument; magnitude now live)**
5. Selection rule / quark flavour / CKM / neutrinos / strong sector (mostly soft or free)
6. Full dynamical Lagrangian (gated on 1+2)

---

## Frontier G9: Induced / emergent metric (highest priority)

**Core question**: Can we recast the v59 scalar connection \(\Omega \in \Lambda^2\) (or its OBE dynamics) as a fundamental 2-form \(B\) on spacetime such that a soldering / simplicity constraint derives a symmetric tensor \(h_{\mu\nu}\) (or tetrad) carrying exactly two transverse-traceless degrees of freedom, while the long-range force still reproduces the second-moment charge \(\rho_{\rm grav} = {\rm Tr}(M^\dagger M)\) and the \(\alpha^{21}\) magnitude?

**Why this is #1**:
- A purely scalar long-range force (N3) is incompatible with LIGO observations of \(h_\pm\).
- The internal \(\mathfrak{so}(8)/\Lambda^2\) index is inert under spacetime little-group rotations.
- Magnitude (P9) and radial law (P8) are moot until the carrier has the right polarization content.
- Suggested avenue (CLOSEOUT): Plebański-style formulation — treat \(\Omega\) (or a projection) as the fundamental 2-form \(B\), impose a simplicity constraint that "solders" it to the Lorentz group, derive \(h_{\mu\nu}\) as an emergent composite.

### Open sub-questions (G9)

**Q-G9-1 (DOF count)**: Does a concrete soldering/simple constraint on the 2-form exist that yields precisely 2 TT DOF (and no more) while preserving the equivalence-principle-exact charge \(\rho_{\rm grav}\)?

**Q-G9-2 (Action principle)**: What is the natural action for the 2-form that reproduces \(\square \Omega = f_g \rho_{\rm grav}\) (or its curved-space generalization) at long range, and also sources the emergent \(h_{\mu\nu}\)?

**Q-G9-3 (Lean's formalization)**: Can the helicity decomposition and the soldering constraint be stated and (partially) machine-checked in Lean (extending `G8G9_Gravity.lean`)?

**Q-G9-4 (Limits of the recast)**: In the flat-space / weak-field limit, does the theory reduce to the v59 scalar OBE + small tensor corrections, or does it modify the inverse-square law at observable scales?

**Q-G9-5 (LIGO falsification)**: If no clean 2-TT-DOF story exists after exhaustive ansatz search, what is the principled exit (modified gravity, composite graviton, short-range-only tensor mode, abandonment of long-range gravity from the algebra, etc.)?

### v60 Roadmap for G9 (session priorities)

**Session 1–2 (setup + first ansatz)**:
- Formalize the Plebański 2-form + simplicity constraint in differential geometry language.
- Write prototype Python/SymPy code for the soldering map and DOF count (linearized).
- Extend `G8G9_Gravity.lean` with a `InducedMetric` module skeleton (helicity + constraint statements).
- **Done (2026-05-25)**:
- `01_8space_to_spacetime_bindings.py` + `01_findings.md` (binding *enumeration* — useful context).
- ~~`02_constrained_helicity_count.py`~~ **CIRCULAR, superseded**: it inserted the spin-2 generator by hand and tensored with an inert identity; the "constraint" did no work. The Lean it cited (`G9ToyHelicity.lean`, "no sorry on spectrum") **did not compile** (bad import path; `decide` on `Multiset ℂ`).
- **`04_soldering_helicity_honest.py` + `04_findings.md`**: the honest result. Soldering = Minkowski sum of co-rotating helicity charges is necessary *and* sufficient for ±2; reproduced as a clean four-scenario contrast (only the co-rotation term differs). The carrier belongs in the **spacetime `Cl(3,1)` bivectors** (already present, = the 6-field Cosserat `so(3,1)`), with `ρ_grav` as the scalar source — not the internal `Λ²(V^8)`. Routes A2/B1 demoted (they attempt the impossible: spacetime helicity from internal indices).
- **`../lean/G9Soldering.lean`** (machine-checked, axiom-only, no `sorry`): `solder_reaches_two`, `soldering_is_the_difference`, `twoform_no_spin2`, `JzTT_charpoly = X²+4`, `JzTT_has_eigenvalue_2i`, and `G9_resolved_if_soldering` (open claim isolated as a hypothesis). `G9InducedMetric.lean`/`G9ToyHelicity.lean` rewritten so they build too.

**Geometric origin synthesis** (new `03_geometric_origin_octo_constraint.md`): traces the required soldering to the concrete G₂ objects already in the v59 octo-space (associative 3-form φ + coassociative 4-form *φ in F-grade, J in L-grade). Shows how these same forms can define the Plebański constraint. Proposes the lift of the current integrated OBE to a local constrained 2-form dynamics. Direct tie-back to octo-space + fundamental equations.
- Lean side: `v60/lean/G9InducedMetric.lean` (all proved structural theorems from v59 restated + new open claims for the constrained 2-form case with explicit `sorry` discipline + minimal lake scaffolding). Full lake build not possible in the current shell (no `lake` in PATH), but the proved parts are verbatim from the already-built v59 module.

**Session 3–5 (attack variants)**:
- Variant A: Direct 2-form on the internal \(\Lambda^2\) bundle.
- Variant B: 2-form valued in the full Cl(7)_even, projected.
- Variant C: Introduce an auxiliary soldering field (tetrad) and integrate it out.
- Test each for 2 TT DOF, EP preservation, and recovery of the v59 radial law.
- Document clean failures with "why it cannot work" theorems (Lean where possible).

**Multi-session (completion or pivot)**:
- If a working recast exists: derive the full curved-space OBE equation and match to P8/P9.
- If none exists after bounded search: write the falsification paper / restructuring proposal for the program.

**Make-or-break test**: Clean count of exactly two transverse-traceless modes + non-zero coupling to \(\rho_{\rm grav}\) at the observed strength.

---

## Frontier G1: Rank tension resolution (second priority, runs in parallel)

**Core question**: The EW bridge requires a democratic full-rank \(Y \in {\rm End}(L)\) with \(\|Y\|_F^2 = 784 a_\ell^2\). The physical spectrum + gravity charge is a rank-3 object with \(\|Y\|_F^2 = 9 Q a^2 = (6/784) v\). These cannot be the same matrix. How is the tension resolved without destroying either the 0.07% bridge match or the 3-generation spectrum?

**Status from gaps/SYNTHESIS**: "No single-step \(\mathfrak{so}(8) \to H\) breaking gives 3 light directions (maximal proper subalgebras top out at 21)."

**Suggested avenue (CLOSEOUT)**: Two-piece \(Y\).
- Active rank-3 block = the physical mass matrix (identified with the Lean-proven `XiVacuum` Higgs).
- The remaining 25 directions are an \(\mathfrak{so}(8)\) Goldstone / gauge sector (a breaking *chain*, not a single step).
- Read R1 (the 0.07% bridge) as a shared conjecture: "the scale is the Frobenius² of the mass bilinear" — consistent across sectors.

### Open sub-questions (G1)

**Q-G1-1 (Two-piece split)**: Can we exhibit an explicit decomposition \(Y = Y_{\rm rank3} \oplus Y_{\rm Goldstone}\) (or block form) inside \({\rm End}(L) \cong M_{28}(\mathbb{R})\) that preserves the 784 count for the bridge while letting only 3 directions acquire vevs and masses?

**Q-G1-2 (Goldstone counting)**: What is the precise 25-dimensional representation? Is it a coset \(\mathfrak{so}(8)/H\) for some natural \(H\)? Does it produce exactly the right number of would-be Goldstones to be eaten or to become heavy gauge bosons?

**Q-G1-3 (Lean formalization)**: State the rank-tension conjecture (R1/R2) cleanly in `EwScaleBridge.lean` (or a v60 successor) and prove the generic 784 result while isolating the single `sorry` on the value pin.

**Q-G1-4 (Gravity link)**: Show that the rank-3 block alone sources the correct \(\rho_{\rm grav}\) while the Goldstone block is either inert or contributes only to short-range forces.

**Q-G1-5 (Breaking chain)**: Map the chain \(\mathfrak{so}(8) \supset \cdots \supset H\) onto the SM gauge group + the silent directions, and check consistency with the Pati-Salam / Spin(7) decomposition already proved.

### v60 Roadmap for G1

**Session 1–2**:
- Re-express the tension in representation-theoretic language (28 = 3 + 25? or other split).
- Scan maximal subalgebras of \(\mathfrak{so}(8)\) (already partially done) and classify possible chains.
- Write `v60/gaps/rank_tension/01_two_piece_y.py` + findings.

**Session 3+**:
- Lean module `RankTension.lean` (or extension of EwScaleBridge) proving the generic End(L) result and the absence of a 25-dim subalgebra that could be the light sector.
- Numerical test: does the two-piece split preserve the 6/784 ratio for leptons while allowing quark sectors their own rank-3 blocks?
- If a clean chain exists: integrate with the G9 recast (does the Goldstone sector source tensor modes?).

**Exit criteria**: Either a concrete two-piece + chain model that keeps the 0.068% bridge match, or a sharp statement of why the 784 construction itself must be reinterpreted (e.g., the bridge is not End(L) but something smaller).

---

## Subordinate Frontiers (G3, G7, selection, etc.)

These remain live but are **deprioritized** until G9 + G1 have a plausible path. Work on them continues opportunistically or in parallel by separate agents, but the main v60 thread is G9/G1/Lagrangian.

- **G3 (α derivation)**: RG fixed-point test of the theorem-grade \((c_W, c_R, c_{B-L}) = (5,5,2)\sqrt{\alpha}\) pattern. See `v59/gaps/alpha_couplings/`.
- **G7 (radian-insert)**: Character / spectral origin for \(\cos(2/3)\) using the now-live magnitude readings (max-mixing + equipartition). Skewness reframing already Lean-proved in `LeptonPhaseMagnitude.lean`.
- Selection rule, quark phases, CKM, neutrinos, strong CP: continue the "explain a zero" and pattern searches; strongest lead remains \(\theta_{\rm QCD} \approx 0\) possibly forced by \(J_c\) (color complex structure, theorem-grade).

---

## The Dynamical Lagrangian (gated deliverable)

Once a viable carrier (G9) and a viable mass bilinear (G1) are in hand, the v60 deliverable is:

**Write \(\mathcal{L}\) on \(\mathrm{Cl}(7)_{\rm even}\) (or its appropriate sub-bundle) such that the EL equations are satisfied by the algebraic \(\Omega(x)\) structures previously derived, and the linearized spectrum around the vacua reproduces the Brannen kernels, the gauge content, and the emergent gravity.**

This is *not* the tree-level effective form already written in `v59/LAGRANGIAN.md`. That document is now the **target specification**.

**Status update (2026-05-26, dynamical-Lagrangian loop GEN1)**: the gravity-sector
gate is now partially open. `lagrangian/10_findings.md` shows (SymPy residual 0 +
axiom-free Lean) that the OBE is the **connection-eliminated trace sector of a
first-order parent action** with an independent `Cl(3,1)` connection — resolving
`09`'s open item and fixing the direction `PARENT ⟹ OBE` (and `PARENT ⟹ Plebański`).
See `lagrangian/LOOP_LOG.md` for the per-generation aspect ledger.

**DELIVERABLE MET (2026-05-26, loop complete, GEN1–GEN8)**: `lagrangian/LAGRANGIAN_v60.md`
is the assembled dynamical Lagrangian `ℒ_v60` whose EL equations yield the OBE
`Ω(x)` structure, with a LIGO-viable 2-TT-DOF gravity sector, the Koide-cone
vacuum (Q=2/3 derived), EP-exact coupling, a stable ghost/tachyon-free spectrum,
and genuine nonlinear time evolution (GEN7). Both gates addressed: **G9** by
GEN1–5, **G1** by GEN6. Residual = 4 value-conjectures (α, `v=784a²`, `φ=2/9`,
`f_g`), all isolated as inputs, not dynamical gaps. Full regression 13/13 + 7 Lean
modules clean. See `lagrangian/CLOSEOUT.md`. Next program (v61): nonlinear/curved
backreaction + the EW-vev home (R1).

### Required ingredients (post G9+G1)
- Kinetic term for the multivector field (or the 2-form \(B\)).
- Potential whose minima enforce the Brannen constraint surface(s) and the sector equilibria.
- Gauge interactions (including the silent SU(2)_L gauging of the Goldstones).
- Yukawa / mass terms that realize the two-piece Y.
- Gravity sector consistent with the induced metric.
- Proof (or numerical evidence) that the EL equations admit the static algebraic solutions as exact or approximate backgrounds.

**Lean / formalization target**: a `LagrangianDerivation.lean` module that states the action and derives (parts of) the EL equation, at minimum for the free + quadratic sectors.

---

## Session Entry Points (2026-05-25 onward)

1. **Immediate (this session)**: Read `v59/CLOSEOUT.md` §"Suggested avenues", `v59/gaps/gravity/ALTERNATIVES.md`, and `v59/gaps/ew_scale_bridge/ALTERNATIVES.md`. Choose first concrete ansatz for G9 (Plebański 2-form) and G1 (representation split of the 28).
2. **Lean**: Fork or extend the existing gap modules under `v60/lean/` (or keep importing from the v59 tree while adding new files). Maintain the `AxiomCheck` discipline.
3. **Geometry**: Stand up `v60/gravity_recast/01_plebanski_setup/` with SymPy + numerical linearization for DOF counting.
4. **Documentation**: Keep `v60/notes/2026-05-25-*.md` for raw thinking; promote clean results to `v60/synthesis/FINDINGS_*.md` or dedicated `gravity_recast/` and `lagrangian/` FINDINGS.
5. **Cross-check**: Every candidate must be run through the LIGO polarization test (reuse/extend `g9_polarization_test.py`) and the rank-tension numerical audit (reuse `formalize_bridge.py`).

---

## Success Metrics for v60 Closeout

- **G9 resolved** at the level of a working (or cleanly falsified) induced-metric ansatz with documented 2-TT-DOF count.
- **G1 resolved** at the level of a concrete two-piece Y + breaking chain that preserves the existing bridge precision while accommodating the 3-gen spectrum.
- **Dynamical Lagrangian** written and shown (analytically or numerically) to have the algebraic structures as solutions.
- Full re-audit in `v60/CLOSEOUT.md` with updated status tags, new theorems, and either a clear path to a unified theory or a principled restructuring recommendation.
- All new Lean modules built and axiom-checked; all numerical scripts run clean and archived.

**If G9 or G1 have no viable path after bounded exhaustive search, the closeout will state the strongest possible negative result and the recommended pivot for v61+ (new algebra, composite gravity, pregeometric rethink, etc.).**

---

*This roadmap will be updated in-session as concrete results arrive. The authoritative source of truth for any claim remains the per-experiment `FINDINGS_*.md` and the Lean sources.*