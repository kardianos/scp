# HFKT v3 — Task Authoring Protocol

Every investigation task MUST be written to `v3/tasks/open/<name>.md` before implementation begins. When a task is completed (positive or null result), move it to `v3/tasks/done/` with the results appended.

---

## Required Sections

Every task file must contain ALL of the following sections:

### 1. Goal
A single paragraph stating the physics question being answered. What would we know after this task that we don't know now?

### 2. Foundation References
Explicit references to the relevant sections and equations in `FOUNDATION.md`. Quote the exact assumptions being tested or extended. If the task depends on results from other completed tasks, reference those by filename in `tasks/done/`.

### 3. Technical Requirements
Precise specification of:
- What code must be written (new files, copied infrastructure, modified algorithms)
- What mathematical derivations are needed before coding
- What input data is required (profiles, parameters, prior results)
- What numerical methods will be used (and why those methods)

### 4. Success Criteria
Concrete, numerical conditions that define success. Examples:
- "Q_EM = 1.000 ± 0.001 for proton isospin orientation on the B=1 hedgehog"
- "The linearized bivector equation reduces to □F = J with J^0 matching the known baryon density to 4 significant figures"

### 5. Null Result Criteria
Equally concrete conditions that define a null or negative result. A null result is still a result — it constrains the theory. Examples:
- "If the coupling matrix M²(r) has no negative eigenvalues, the bivector sector does not bind to the soliton"
- "If Q_EM depends on an arbitrary gauge choice, charge is not topologically determined"

State what the null result would MEAN for the theory if it occurs.

### 6. Implementation Plan
Ordered steps. For each step, state what it produces and how to verify it before proceeding to the next step.

### 7. Reusable Infrastructure
List any code being copied or adapted from v2, with file paths and the specific functions/structs being reused. Note any known bugs in the v2 code that must be fixed in the copy (check `v2` memory files for the bug list).

---

## Principles

- **No implementation without a task file.** The task file is the contract. If the goal changes during implementation, update the task file first.
- **Null results are first-class.** A well-documented null result is worth more than an undocumented positive one. Always record what was tried, what failed, and what it rules out.
- **Reference FOUNDATION.md by section number.** The foundation document is the source of truth for what the theory claims. Every task should trace back to a specific claim, open problem, or caveat.
- **State assumptions explicitly.** If the task assumes σ-model (λ→∞), or a specific profile, or hedgehog symmetry, say so. These are the conditions under which the result holds.
- **One physics question per task.** If a task naturally splits into independent sub-questions, split it into separate task files. Tasks may depend on each other (state dependencies explicitly).
- **No v2 source modification.** v2 code may be copied or referenced but never modified in place.
- **Results go to RESULTS.md.** FOUNDATION.md is never modified. When a task resolves an open problem or caveat from FOUNDATION.md, record the result in `RESULTS.md` with the relevant section reference. RESULTS.md shows actionable items only: what was resolved, what remains open.
- **Task files get results appended.** When a task is completed (positive or null), move it from `tasks/open/` to `tasks/done/` and append a `## Results` section with the numerical outcomes, conclusions, and any follow-up items identified.
