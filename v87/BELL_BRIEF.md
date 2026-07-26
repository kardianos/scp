# v87 Rung B0 — the Bell problem, stated as a solvable symbolic task
**Shared brief. Read this first. Both working agents receive it.**

---

## 0. What is being asked, in one paragraph

The SCP programme's ontology is a **deterministic sub-quantum field theory**:
there is a definite fabric state at all times, and quantum probability is
**epistemic** — it arises because uncontrolled degrees of freedom of the fabric
are not held fixed between runs. The task is to establish, **symbolically and
then numerically, and completely**, whether a model obeying the programme's two
ontological requirements can reproduce CHSH correlations **above the local bound
of 2 and bounded at or below the Tsirelson value 2√2** — and if so, exactly what
structure and what event source that requires, and exactly what it costs.

This is **not** a request to refute Bell's theorem. Bell's theorem is a theorem
and it is correct. It is a request to locate, with complete precision, which of
its assumptions the fabric ontology relaxes, to build an explicit model in that
relaxation, and to quantify the relaxation required. An answer that claims to
exceed 2 while keeping all three assumptions is wrong by construction and will
be discarded.

---

## 1. The two ontological requirements (hard constraints; a model violating either is out of scope)

**O1 — Definite substrate, epistemic probability.** There is a single definite
fabric configuration evolving by definite dynamics. There is one history. No
amplitudes, no branching, no fundamental randomness. Observed probabilities come
from ignorance of, or lack of control over, fabric degrees of freedom — in the
same way that a car driving a circular track at unknown phase yields a hit
probability for a laser fired at an uncontrolled time, while the underlying
description remains "the car is driving in a circle."

**O2 — Field monism.** Particles are not entities separate from the fabric; they
are manifestations of it — localised, persistent field structures. There is no
second substance, no separate stage, and no observer standing outside the
fabric. Detectors, settings-choosers and experimenters are themselves fabric.

Note that O2 is not decorative: it bears directly on the argument below.

---

## 2. Bell / CHSH, stated exactly

Two spacelike-separated measurements. Settings a, b. Outcomes A, B ∈ {−1, +1}.
Assume:

- **(R) Realism / hidden variables.** Outcomes are fixed by a variable λ drawn
  from ρ(λ): A = A(a, λ), B = B(b, λ).
- **(L) Locality / parameter independence.** A does not depend on b; B does not
  depend on a.
- **(MI) Measurement independence / free settings.** ρ(λ | a, b) = ρ(λ).

Define E(a,b) = ∫ dλ ρ(λ) A(a,λ) B(b,λ) and

    S = E(a,b) + E(a,b′) + E(a′,b) − E(a′,b′).

Then **(R) ∧ (L) ∧ (MI) ⟹ |S| ≤ 2** (CHSH). Quantum mechanics gives
|S|max = 2√2 (Tsirelson). Loophole-free experiments (Hensen 2015, Giustina 2015,
Shalm 2015) confirm |S| > 2. Purely no-signalling correlations would permit
|S| = 4 (Popescu–Rohrlich boxes); nature does not do this.

**Therefore exactly one of (R), (L), (MI) must be relaxed, and the model must
explain why the result stops at 2√2 rather than running to 4.** Both halves are
required. A model that reaches 2√2 but has no principle bounding it below 4 has
solved half the problem and must say so.

---

## 3. The three doors, and why O1/O2 pre-select

**Door A — relax (MI): measurement dependence.** O2 says the setting-chooser is
made of fabric and shares a causal past with the source. Statistical
independence of λ from (a, b) then looks like a smuggled dualism rather than a
neutral assumption. This is the door the stated ontology motivates on its own
terms. Retains determinism and locality; costs "free choice."
 *Known anchor:* M. J. W. Hall, PRL 105, 250404 (2010) and PRA 84, 022102 (2011)
 — an explicit local deterministic model reproducing singlet correlations with
 only partial measurement independence, with the required degree QUANTIFIED
 (a fraction of ~1/15, i.e. ≈0.0663 bits of mutual information, for the full
 singlet). Barrett & Gisin, PRL 106, 100406 (2011), quantify how much free
 randomness is needed. **These are the results to reproduce, generalise, and
 then attempt to DERIVE from fabric dynamics rather than postulate.**

**Door B — relax (L): a nonlocal element.** de Broglie–Bohm, and specifically de
Broglie's *double solution* (1927), in which a physical field carries a
localised singular/solitonic region that IS the particle — precisely O2. Costs a
preferred foliation. **Note: the SCP kernel already exhibits a preferred frame
at the 1–5% group-velocity level (v86 THEORY A2), currently recorded as a
defect.** Under Door B that is the required structure, not an artifact.
 *Known anchors:* de Broglie's double solution; its modern revival (Colin, Durt,
 Willox, 2017); Couder–Fort walking-droplet analogues.

**Door C — relax the temporal ordering: retrocausality / all-at-once.**
Two-state-vector formalism; Price; Wharton & Argaman, Rev. Mod. Phys. 92, 021002
(2020). Retains locality in a spacetime-block sense. Costs the assumption that λ
is fixed prior to the settings.

**The task is not to pick a favourite. It is to work each door far enough to say
which one O1 ∧ O2 actually admits, and at what cost, with the algebra done.**

---

## 4. The Tsirelson half of the problem

Why 2√2 and not 4? Any candidate answer must be examined:

- **Information causality** — Pawłowski, Paterek, Kaszlikowski, Scarani,
  Winter, Żukowski, Nature 461, 1101 (2009).
- **Macroscopic locality** — Navascués & Wunderlich, Proc. R. Soc. A 466, 881
  (2010).
- **Local orthogonality / exclusivity** — Fritz et al., Nat. Commun. 4, 2263
  (2013).

A substrate model that *derives* 2√2 from its own structure — rather than
importing it — would be the strongest possible outcome of this rung. If that is
not achieved, say so explicitly and state what is missing.

---

## 5. Required deliverables

1. **Symbolic.** For the chosen door(s), an explicit model: the fabric variable
   λ and its state space; ρ(λ) and, for Door A, ρ(λ | a, b); the response
   functions A(a, λ), B(b, λ); the event source (what "produces" a run and what
   is or is not controlled between runs). Then E(a,b) evaluated **in closed
   form**, and S(a, a′, b, b′) maximised **symbolically**. Use Maxima (installed:
   5.46.0) or SymPy (1.12) and **commit the actual worksheets**, not just
   results.
2. **The bound.** Prove or bound max|S| for the model. State whether it reaches
   2√2, falls short, or overshoots toward 4. Overshoot is a *failure mode* to be
   reported, not hidden — a model that gives 4 is unphysical.
3. **The cost, quantified.** For Door A: the mutual information I(λ : a,b) in
   bits, versus the Hall/Barrett–Gisin bounds. For Door B: the required
   foliation and its signature. For Door C: what boundary condition is imposed.
4. **Numerical.** Monte-Carlo the model with an explicit pseudorandom event
   source, recover S numerically, and confirm it matches the symbolic result to
   stated precision. Include the code.
5. **The derivation question — the crux.** Can the required structure be
   *derived* from fabric dynamics obeying O1 ∧ O2, or must it be postulated? A
   postulated correlation between λ and settings is a parametrisation of the
   answer, not an explanation. Say clearly which you achieved.
6. **A falsifiable consequence.** Any real relaxation leaves a fingerprint —
   a residual setting-correlation, a preferred-frame signature, a bound on
   detector-randomisation timing. State at least one, quantitatively.
7. **Honest negative results are full deliverables.** If a door closes under
   O1 ∧ O2, prove the closure. That is as valuable as an opening.

---

## 6. Standards

- Quote what you attack; number your findings.
- Distinguish **derived**, **measured**, **postulated**, **conjectured** on every
  claim. The parent programme tags these [D]/[M]/[C]; do the same.
- No result without its error bar, residual, or explicit exactness statement.
- Cite real literature with real identifiers. **If you are not certain a
  reference exists, say so** — do not manufacture citations. A fabricated
  reference invalidates the work it supports.
- If you find the task impossible as posed, say that and prove it.

## 7. Context

The parent programme is at `/home/d/code/scp`. Relevant reading, in order:
`v86/THEORY_v86.md` (the ontology and the measured state — see A1, A2, A9),
`v86/PART0_RESULTS.md` (what has actually been measured, and to what precision),
`v86/NEXT_PROGRAM.md` (the open walls: fermions, ħ, charge quantisation, Bell,
and the observation that topology is the common ingredient in the first three).

You do **not** need to run the simulation kernel. This rung is symbolic and
numerical analysis. Do not modify anything under `sfa/` or the v86 results.
