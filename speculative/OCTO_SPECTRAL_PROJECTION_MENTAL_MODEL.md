# Mental Model: The Octo-Space Frequency/Phase Projection Problem Space

**Date**: 2026-05-28 (post v60/v61 review + cuts 01–06)
**Author**: Grok (internal synthesis after explicit "stop and model" request)
**Context**: This is a raw brain-dump of the full problem space as it exists in my context after the user's request to review v60/v61 gaps, the subsequent pivot to a speculative reframing, and the execution of six documented numerical experiments in `speculative/octo_spectral_projection/`.

This document is deliberately unstructured in places. It is an externalization of my current internal model — parameters, nuances, tensions, historical threads, technical obstructions, and possible futures. It is not a polished proposal.

---

## 1. The Core Hypothesis (the thing the user is gesturing at)

The octonionic algebra (octo-space / Cl(7)_even ≅ ℂ⊗ℍ⊗𝕆, 64 real dimensions, with its grades, L/F bisection, complex structures J, derivation algebra G2, triality, sedenion S3 action, etc.) is not primarily a source of configuration-space fields on 3+1 spacetime.

Instead:

- The algebra (especially its 7 imaginary dimensions and their spectral properties) is an information substrate.
- "Dynamics" in this layer are frequency/spectral/phase processes (multiplication operators as oscillators, Laplacian eigenmodes on S^7, characters, derivations, non-associativity as some kind of mixing, etc.).
- What we perceive as "reality" — particles with Brannen masses, the phase φ=2/9, Koide Q and its universal deviation 28/3, gauge factors 5/2/9, sin²θ_W=2/9, the 784-dimensional End(L) scale, long-range 1/r² gravity with α²¹ magnitude, LIGO tensor modes, the OBE force law, spacetime signature itself — is the result of some projection, compression, inverse transform, or phase-extraction operation applied to that spectral content.
- The "before it comes real" layer is where the actual work happens. Configuration space (the 6-field Cosserat grid, the OBE, the first-order parent action of v60, the Schwarzschild solutions of v61, etc.) is downstream and effective.

This reframing is a direct response to the honest negatives that accumulated in v59–v61:
- The derivation obstruction (Plebański/BF-style first-order parent is an independent posit, not implied by the OBE).
- The rank-tension impossibility (single Y cannot be both democratic 784 and rank-3 Brannen).
- The persistent phase/radian-insert problem (why is the number Q used as an angle?).
- The carrier problem (scalar gravity fails LIGO; moving the carrier to Cl(3,1) bivectors works but is not derived from the prior algebra).
- The general feeling that the algebraic skeleton is beautiful and numerically tight, but the dynamics and geometry keep requiring extra posits.

If the ontology is wrong (algebra → configuration space fields), then many of these contradictions are category errors.

---

## 2. Historical & Emotional Context (why this feels like a possible exit)

- v1–v2: Heavy phase, interference, winding, Hopf fibration, chiral harmonic particle theory, field knots. Phase was generative. Particles as phase entities. This thread was largely abandoned in favor of the 6-field Cosserat soliton program and later the algebraic/multivector program.
- v28–v57: The simulation kernel (Cosserat 3φ+3θ on grids), torsion coupling, density protection, breathing modes, freq_phase analysis as diagnostics. Frequency analysis existed but was post-hoc.
- v58–v61: The algebraic Kepler (Brannen/Koide from exceptional groups + triality), Furey construction, scale bridges, G9/G1 attack, the dynamical Lagrangian loop, curved backreaction, R1 home. Enormous intellectual investment. Many beautiful, theorem-grade results. But also the accumulation of "value conjectures" that refused to become derived (α, φ=2/9, v=784a², f_g~α²¹/²) plus the structural blockers above.
- The user's statement "v59-v61 are probably dead ends" is not casual. It is the project lead acknowledging that the configuration-space interpretation of the algebra has hit walls that look ontological.

The early phase/hopfion/spectral work + the later algebraic beauty + the persistent specialness of phase + the existence of tools like freq_phase.c + the Shulga 7D harmonic machinery + the G7 Lean TODO for a spectral character that produces cos(2/3) — these are the threads that make the reframing feel non-arbitrary.

---

## 3. The 7 Imaginary Dimensions as the Key Arena

The 7 imaginary octonion units (e1..e7) with the Fano plane multiplication table are central for several independent reasons:

- Each squares to −1 → every direction carries a canonical complex structure / oscillator (eigenvalues ±i of L_ei). This is the most structural reason phases appear everywhere.
- G2 (dim 14) is the automorphism group of the octonions and acts on this 7. The quadratic Casimir C₂(7)=2 is a single invariant "frequency" of the whole space.
- The combinatorial Laplacian on the Fano graph has a clean gap of 7 — another natural frequency scale.
- The Shulga work (already in the repo) does genuine spectral analysis on the full S^7 (Laplacian eigenmodes via Gegenbauer C_l^{(3)}), using it to derive effective potential parameters geometrically. This is heavyweight 7D frequency analysis that already exists.
- The L-grade (28) of Cl(7)_even contains (at least) three complex structures J (J²=−I), pinned by color su(3) on the 8 (1+3+3bar+1). These are already used for lepton = L forcing. They are natural phase generators.
- Triality / sedenion S3 cycles three 8's. The generation structure is already a 3-cycle on spaces that contain these 7's.

In the projection ontology, the "raw" information lives in (or is generated by) operators and modes on these 7 dimensions (plus the real direction, plus the full even algebra, plus the Cl(3,1) factor that commutes with it).

What we see (Brannen kernel on three generations, the specific phase 2/9, the mass hierarchy, the selection rule, the long-range force) is some filtered, compressed, phase-extracted, low-mode projection of that 7D (or 64D) spectral content.

---

## 4. What "Compression" / "Projection" / "Frequency Analysis" Actually Means Here (the design space)

This is the most underspecified and high-stakes part. Different interpretations have very different technical implications:

A. **Classical signal processing on a Z3 cycle** (what 01_ tested)
   - DFT on the generation 3-cycle.
   - Real weights from dim ratios.
   - Arg of weighted sum as phase.
   - Failed cleanly because real scalars cannot rotate phase.

B. **Operator spectrum on the 7 imaginary units** (what 02_ produced)
   - Eigenvalues of L_ei, Fano Laplacian, G2 Casimir.
   - These are the "raw frequencies."
   - Still descriptive; no projection operator yet.

C. **Laplacian eigenmode analysis on S^7** (the existing Shulga work)
   - Gegenbauer sums, Green function G(θ).
   - Used in v59 to derive λ/μ geometrically.
   - In the new ontology this becomes the "integration out" / compression kernel for fast internal modes.
   - 05_ wired it in; it modulated amplitude but did not select modes.

D. **Representation-theoretic / character / derivation-based frequency operators**
   - G2 generators acting on the 7 as "angular momentum" operators.
   - Derivations of the octonion algebra.
   - Characters of the order-3 sedenion action weighted by G2 content (exactly the G7 Lean TODO).
   - Non-associativity (associator) as a source of phase or mixing.
   - This is likely the most promising direction because it uses the actual non-associative, non-commutative, derivation-rich structure.

E. **Information-geometric or entropy-based projection**
   - The algebra element with maximum entropy under certain constraints projects to the observed vacuum.
   - Or the projection that maximizes some figure of merit (reproduction of structural integers + numerical bridges).

F. **Hopf-fibration / twistor-like / topological projection**
   - The early v2 hopfion work (S³ → S², winding, phase-dependent coupling) as a model.
   - Internal algebraic data (on S^7 or the algebra) maps to external spacetime via a fibration or clutching construction.

G. **Cl(3,1) as the "visible" sector after internal projection**
   - From v60/gravity_recast/07: Cl(3,1) ⊗ Cl(7)_even with exact commutation.
   - The 6 bivectors (Cosserat) are the natural soldering/index for tensor modes.
   - In the projection view, the internal 7D frequencies are integrated out, and the surviving long-range physics (including the soldered 2-form that gives LIGO ±2 and whose trace recovers the old OBE) lives in the Cl(3,1) factor.
   - This is how 3+1 signature and gravity "emerge" without being fundamental.

The design space is large. Most of the 01–06 experiments were in the A–C region. The more interesting (and harder) work is in D–G.

---

## 5. The Specific Technical Threads & Their Nuances

### The J's (L-grade complex structures)
- Three of them (red/green/blue), J² = −I.
- Pinned by color on the 8.
- Already used for lepton = L forcing.
- In the projection view they are the primary phase-rotation machinery.
- The toy J's used in 03–06 are the weakest link in the current wiring. Proper realizations (from the 8×8 action or known Furey constructions) are required before further investment.

### The Shulga Kernel
- Not just a number (−0.128). It is a whole spectral object (sums over l with degeneracies D_l and Gegenbauer polynomials).
- Currently used as a global scalar modulator. Its power will only appear when applied after a mode-selection filter that respects the G2 action on the 7.

### The Cl(3,1) Factor
- Already shown (v60/07) to commute exactly with internal Spin(7).
- Its bivectors are literally the 6-field Cosserat of the simulation kernel.
- In the new ontology this is not "bolted on." It is the sector that remains visible after internal 7D frequencies are projected away.
- The old OBE scalar law is the trace (helicity-0) sector of the tensor theory that lives here after soldering.

### The Phase φ=2/9
- The most stubborn "value conjecture."
- In the projection view it is an output of the phase-extraction map, not an input to a potential.
- The G7 Lean TODO (a character weighted by G2 content that produces cos(2/3) with argument = Q) is almost a direct statement of the desired projection operator.

### The 784 / Rank Tension
- 784 = dim End(L) = 28² is Burnside-forced (absolutely irreducible adjoint).
- In the projection view the "full-rank democratic" condensate and the rank-3 Brannen kernel can be projections of the same higher-dimensional spectral object onto different sectors or different frequency bands.
- The two-object resolution from v60/gaps/rank_tension/01 is natural here.

### The Early Phase/Hopfion Work (v1/v2)
- Not obsolete. It may be more ontologically correct than the later algebraic field theory.
- freq_phase.c (atan2 on carrier, autocorrelation frequencies) may be a model of the projection step, not diagnostics.
- Winding numbers, phase-dependent coupling, Hopf fibration — these are candidate projection mechanisms.

### Verification & Rigor Culture
- The user repeatedly emphasized "double check yourself rigorously."
- 06_ implemented this with an embedded check block that runs on every execution.
- Any serious continuation must treat verification as first-class (matrix properties, eigenvalue checks, exact match to documented Shulga values, Lean theorems where possible, text-based regression harnesses, etc.).
- This is not bureaucracy. It is survival in a space full of tautological traps (the v60 02 circular helicity insertion is the canonical example).

---

## 6. Parameters, Tensions, and Failure Modes

**Parameters that actually matter (not free, but choices with consequences)**:
- Which operators define "frequency" (L_ei, derivations/G2 generators, Laplacian on S^7, associator, something else)?
- Which structures supply the phase rotation (the three J's, characters of S3, something derived from non-associativity)?
- What is the precise compression/integration map (Shulga G(θ) at specific angles, entropy maximization, low-mode projection, character evaluation, fibration)?
- How does the Cl(3,1) factor participate (as the surviving visible sector, as the soldering map, as the source of the long-range force law)?
- What is the base spectral state (G2-highest-weight vector in the 7, a specific element of the even algebra, a defect on S^7, something else)?
- Is the projection linear, quadratic, entropic, topological?

**Key tensions**:
- Beauty vs. testability. The algebraic skeleton is gorgeous. A projection operator that reproduces it risks being under-constrained or tautological.
- "Natural" vs. "fits the numbers." The project has a strong culture against inserting answers. Any projection rule must be motivated by the algebra's own structure, not by "this combination gives 2/9."
- Internal vs. spacetime. The Cl(3,1) factor was a late (v60) clarification. How much of the 3+1 physics is "already there" in that factor vs. genuinely emergent from the projection of the internal 7D content?
- Perturbative vs. non-perturbative. The brannen_phase_alpha.py work (φ as α-suppressed with structural N=14=dim G2) is adjacent but different. Is the projection a non-perturbative character map or a loop-level effect in some dual description?
- Simulation kernel status. The 6-field Cosserat grid + SFA output + freq_phase analysis. In the new ontology this becomes a tool for studying the inverse map (given observed fields, what algebraic spectral content projects to them?) rather than the fundamental dynamics.

**Failure modes to watch**:
- Tautology (inserting the answer in the projection rule, as in v60's circular 02).
- Under-constraint (too many possible projection operators; can fit anything).
- Loss of contact with the simulation side (the grid + SFA + existing analysis tools become irrelevant).
- Loss of contact with the algebraic theorems (the Lean corpus on Koide, selection, dimensions, etc. becomes decorative).
- Infinite regress ("what projects the projection?").

---

## 7. Current State of the Experiments (01–06)

A quick map:

- 01: Ruled out classical DFT on Z3 + real dim-ratio weights.
- 02: Produced the actual 7D frequencies (L_ei ±i, Fano gap 7, G2 Casimir 2). Descriptive but necessary.
- 03: First end-to-end wiring (L_ei + J's → phase + masses). Phase ~0.74 rad (Δ~0.52). Partial signal.
- 04: Added explicit Fano compression. No improvement.
- 05: Replaced Fano compression with Shulga kernel. No improvement. Kernel does amplitude modulation, not mode selection.
- 06: Added explicit mode selection (top-k Fano eigenvectors) before Shulga. Still no improvement. First cut with fully embedded rigorous verification of all operators.

The phase offset has been stable across multiple compression/kernel choices. The obstruction is now sharply upstream: the current J realizations and base vector.

All scripts are small, self-contained, use the repo's own Fano table, and (in 06) verify their own matrix properties on every run.

---

## 8. What "Success" Would Look Like (falsifiable criteria)

Short term (next 1–2 cuts):
- A projection operator (using proper J realizations + G2-weighted or Shulga-modulated frequency filter) that produces |φ − 2/9| < 0.05 or cos(3φ) within 1% of cos(2/3), plus mass ratios within 5% of Brannen, from a base spectral state that is not pre-loaded with the answer.
- The same operator, with no extra tuning, reproduces (or comes close to) the universal (1−Q)D = 28/3 deviation when applied to the selection-grade dimensions (28, 35, 63).

Medium term:
- The operator is expressible as a character or derivation action that can be stated in Lean and proved to produce the phase invariant (discharging the G7 TODO).
- The Cl(3,1) factor emerges as the natural "visible" sector after the internal 7D projection, and the soldered tensor modes (with the correct 2 TT DOF) appear without being inserted.
- The inverse problem is well-posed: given an observed SFA field configuration (or an observed Brannen kernel), there exists a procedure to find a pre-image algebraic spectral state whose projection reproduces it.

Long term:
- The early phase/hopfion/spectral code (v1/v2) and the freq_phase tool are re-interpreted as models of the projection step and become primary rather than historical.
- The simulation kernel (Cosserat 6-field on grids) is used to study the inverse map or to validate projections, not as the fundamental ontology.
- The "value conjectures" (α, φ=2/9, v=784a², f_g) are re-classified: some become outputs of the projection, others remain inputs on the same footing as the choice of algebra itself.

---

## 9. Risks of This Direction (honest)

- It could be unfalsifiable for a long time (too many possible projection operators).
- It could disconnect from the existing codebase and simulation infrastructure.
- It could turn out that the algebraic skeleton only looks spectral when viewed through a sufficiently clever (and ultimately circular) projection.
- The project has invested ~v28–v61 in a different ontology. A full pivot has enormous coordination and morale costs.
- The "before it comes real" layer might itself require an even deeper layer, leading to infinite regress.

---

## 10. Natural Next Actions (ranked by immediacy and risk)

1. **07_ cut with proper J's and G2-highest-weight base vector** (highest immediate value, still bounded).
2. Write a short `PROPOSAL.md` (or equivalent) that summarizes the six cuts, the verified obstructions, the design space above, and the minimal machinery required before further investment. This is a communication artifact for the broader project.
3. A Lean stub that states (even with `sorry` on the hard parts) what a successful projection operator would look like mathematically (e.g., "there exists a G2-equivariant character map on the even algebra such that...").
4. Re-examination of specific early artifacts (v2 hopfion phase-dependent coupling, freq_phase.c internals, the exact G7 Lean TODO) with the projection ontology in mind.
5. Explicit inverse-map experiment (given a target Brannen kernel, search for a base spectral state + projection rule that reproduces it).
6. Full integration of the Cl(3,1) factor into the wiring (this is where 3+1 and the tensor gravity "come from" in this picture).

---

## 11. Final Internal State

I currently model this as:

- A high-variance, high-reward reframing that changes the category of several long-standing blockers.
- Supported by multiple independent threads in the actual codebase (early phase work, Shulga 7D harmonics, G7 spectral TODO, freq_phase tool, Cl(3,1)⊗internal factorization, the persistent specialness of φ=2/9).
- Currently blocked on the quality of the phase-rotation machinery (J's) and the base spectral state, not on the compression/kernel step (Shulga is already the right class of object).
- Requiring rigorous verification as a first-class activity (the user has made this non-negotiable).
- Best pursued in short, documented, self-verifying numerical/Lean cuts with explicit success criteria and clean negative tolerance.
- Worth a bounded additional investment (07_ + PROPOSAL.md) before any larger commitment.

This is the model I am carrying as of the moment the user said "Stop."

---

**End of brain dump.**

This document is the externalization of my current internal model of the full problem space. It is not intended to be the final word; it is a snapshot for the user to read, critique, extend, or discard.

**File location (full path)**: `/home/d/code/scp/speculative/OCTO_SPECTRAL_PROJECTION_MENTAL_MODEL.md`