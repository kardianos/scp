# Unified Multivector Force Law Experiment

**Location**: `v58/pregeometric/unified_multivector_force/`  
**Date Started**: 2026-05-18  
**Goal**: Discover and rigorously verify a single multivector equation (or small family of equations) whose consequences, in appropriate limits and projections, recover:

**Key Reference Documents** (created for long-term continuity):
- `TERMINATION_CRITERIA_AND_CURRENT_STATUS.md` — Precise, measurable definition of success and the current concrete gap.
- `RESUME_HANDOFF.md` — Practical guide for anyone resuming this experiment later.
- Newtonian gravity with ambient-density dependence (as developed in the v58 notes)
- Coulomb's law
- Maxwell's equations

while naturally incorporating causal structure / regularization so that interactions remain local at speed \(c\) and do not produce unphysical infinite tails.

This experiment runs two tightly coupled tracks in parallel:
- **Python / Geometric Algebra track** — rapid discovery, numerical testing, and shape exploration.
- **Lean 4 track** — formal encoding and proof of implications / equivalences to known physics.

The two tracks are expected to give continuous feedback to each other.

---

## 1. Overall Philosophy and Approach

We are working inside the multivector / projective geometric algebra ontology developed across the v58 series (especially `MULTIVECTOR_DENSITY_GRAVITY.md`, `MEDIUM_DYNAMICS_TAILS_AND_WAKES.md`, `PARTICLES_AS_DENSITY_ACHIEVERS.md`, and `MULTIVECTOR_FORCE_LAW.md`).

The guiding intuition is that both gravity and electromagnetism are different projections or response modes of the **same underlying multivector field** \(M\) (or the pre-geometric substrate that gives rise to it). Gravity is primarily carried by scalar-like density gradients; electromagnetism is carried by phase, rotational, and bivector degrees of freedom (chiral currents / twists).

A successful unified equation should be **single** at the fundamental level, with the separation into "gravity" and "EM" emerging only after taking limits or projecting grades.

Because we are ultimately aiming at a pre-geometric or causal-set-like foundation, the equation must not presuppose a fixed background manifold with built-in infinite-range propagation. Some form of causal regularization (retarded kernels, algebraic partial orders, or quantum-like operators) is required.

---

## 2. Two Parallel Tracks and Their Roles

### Python / Geometric Algebra Track (Discovery & Intuition)
**Purpose**: Explore the *shape* of possible unified equations quickly.

**Primary activities**:
- Implement multivector fields (especially in \(\mathbb{R}(3,1,1)\) or similar PGA).
- Write candidate single equations (see `BACKGROUND_AND_SPECULATIVE_EQUATIONS.md`).
- Symbolically project grades to recover Newtonian + Maxwell limits.
- Numerically evolve small systems (lattices, graphs, or causal-set-like structures) with retarded interactions to observe emergent \(1/r^2\) gravity and Maxwell-like behavior.
- Test different forms of ambient-density modulation \(f(\rho_{\rm ambient})\) and chiral current definitions \(\mathcal{J}_\chi\).
- Use visualization (ganja.js, custom 3D) to develop intuition about how density gradients and bivector twists interact.

**Output artifacts**:
- Python scripts / notebooks in `python/` sub-folder.
- Regular updates to `PYTHON_FINDINGS.md` (structured findings, promising equation forms, numerical results, dead ends).

### Lean 4 Track (Rigor & Constraint)
**Purpose**: Formally encode promising equations and prove what they do (or do not) imply.

**Primary activities**:
- Encode the multivector algebra (or a sufficient fragment) in Lean.
- State the candidate unified equation as an assumption.
- Prove theorems of the form: "Under these conditions, the equation implies the Maxwell equations" and "implies the Newtonian limit with density dependence".
- Explore the space of necessary extra assumptions (specific kernels, specific forms of \(f(\rho)\), specific regularization operators) by proving what must hold for the known laws to emerge.
- Maintain a "Known Physics" module that states Newtonian gravity + Maxwell as separate statements, then prove implication or equivalence where possible.

**Output artifacts**:
- Lean code in `lean/` sub-folder (following the project's existing Lean conventions).
- Regular updates to `LEAN_FINDINGS.md`.

---

## 3. Coordination Protocol (How the Tracks Talk to Each Other)

The two tracks must not work in isolation. They are expected to give each other actionable feedback on a regular cadence.

**Mechanism**:
- Both agents (and the coordinator) write to shared files in the experiment root:
  - `PYTHON_FINDINGS.md`
  - `LEAN_FINDINGS.md`
  - `COORDINATION_LOG.md` (high-level synthesis and next priorities)
- When one track makes a significant discovery or blockage, it should explicitly call out implications for the other track (e.g., "This form of kernel appears to produce the correct Coulomb law numerically; Lean should try to prove whether it also yields the full inhomogeneous Maxwell equations under these algebraic assumptions").

**Cadence** (suggested):
- After every major Python experiment or Lean proof, the agent writes a short structured update.
- The coordinator (or the other agent on the next cycle) reads the other's latest findings and incorporates them into its next prompt or work plan.

**Example feedback loops**:
- Python finds numerically that a certain nonlinear term is needed to get the correct \(1/r^2\) fall-off for gravity while preserving Maxwell linearity → Lean adds that term as an assumption and tries to prove the limits.
- Lean proves that a particular form of ambient-density function \(f(\rho)\) is required for consistency with the equivalence principle → Python tests whether that same function still allows interesting large-scale modifications (e.g., for the dark-matter-like effects discussed earlier).

---

## 4. Initial Speculative Equations

Detailed speculative candidate equations, background motivation, and decomposition strategies are collected in the companion document:

**`BACKGROUND_AND_SPECULATIVE_EQUATIONS.md`**

Agents should treat that document as the current working set of ideas. New candidates discovered during the experiment should be added there (with date and rationale) and cross-referenced in the findings files.

---

## 5. Tools Stack

### Python / Geometric Algebra
- `python3` + `pip install clifford sympy numpy matplotlib`
- `ganja.js` (browser-based PGA visualizer) — extremely useful for intuition about bivectors and rotors.
- Optional: custom graph / causal-set simulators using `networkx` + custom multivector class.
- Jupyter notebooks or plain `.py` scripts (both acceptable).

### Lean 4
- Follow the project's existing Lean setup (see root `lean/` directory).
- Use `lake` for building.
- Mathlib4 is available and should be used for algebra, analysis, and geometry where helpful.
- Keep files in `lean/` sub-folder with clear module structure.

### Coordination & Documentation
- All long-form writing goes in Markdown files in the experiment root.
- Agents are expected to write clearly structured updates so that a human (or the other agent) can quickly understand the current state without re-reading everything.

---

## 6. Success Criteria (Phased)

**Phase 1 (Discovery)**: Identify at least one candidate single multivector equation (or small parametric family) that numerically produces both inverse-square gravity and Maxwell-like propagation on small test systems, with plausible causal regularization.

**Phase 2 (Partial Proof)**: In Lean, prove that the candidate implies the inhomogeneous Maxwell equations and the Newtonian limit (with at least a simple form of ambient-density modulation) under clearly stated extra assumptions.

**Phase 3 (Causal Structure)**: Show (numerically + ideally with some Lean support) that the equation can be made to respect a local causal structure (no superluminal influences, finite propagation speed \(c\)) without destroying the classical limits.

**Phase 4 (Integration)**: Demonstrate consistency with the other v58 concepts (particles as density + chiral achievers, tails/wakes, ambient-density effects on both gravity and EM strength).

---

## 7. Immediate Next Actions for the Agents

Both agents should begin by reading the following two documents in full:

1. `BACKGROUND_AND_SPECULATIVE_EQUATIONS.md` (current speculative landscape)
2. This `EXPERIMENT_OUTLINE.md`

Then:
- **Python agent**: Begin by implementing a basic multivector field + at least one of the speculative equations from the background document. Run a minimal symbolic projection or numerical test and report findings.
- **Lean agent**: Begin by setting up the module structure and encoding the multivector algebra (or a useful fragment). Then pick one speculative equation and start stating the first implication theorems.

Agents should write their first structured updates to their respective findings files within the first few cycles.

---

*This experiment is deliberately open-ended and exploratory. The value lies as much in the process and the constraints discovered as in any final "unified equation." All findings, dead ends, and partial results are valuable.*