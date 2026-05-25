# The Lepton Phase Law  φ = Q/3

*Theory document — the current cohesive understanding of why the charged-lepton
Brannen phase equals one third of the Koide ratio, what is established, and what
a physical explanation would have to look like.*

## 1. Statement

The charged-lepton masses are reproduced by the **Brannen circulant mass matrix**

>   M = a (I + ξ S + ξ̄ S²),   ξ = t·e^{iφ},   S = the 3-cycle (Z₃ generation shift),

whose eigenvalues are the √m_k.  The matrix has two dimensionless shape
parameters: the **amplitude** t = |ξ| and the **phase** φ = arg(ξ).  Empirically
both are fixed by the single Koide ratio Q:

>   **amplitude:**  t² = (3Q − 1)/2    (this *is* Koide: Q = (1+2t²)/3, so t²=1/2 ⇔ Q=2/3)
>   **phase:**      φ = Q/3            (the law studied here; lepton value φ = 2/9)

Equivalently, in terms of the complex order parameter ξ:

>   **9·arg(ξ) = 1 + 2|ξ|²**     ( = 2 at the Koide point |ξ|² = 1/2 ).

So the law locks the *phase* of ξ to its *modulus*.  The factor 9 = (generations)²,
or read as 3·3 — one 3 from the definition Q = (1+2t²)/3, one from the law φ = Q/3.

## 2. Empirical status — positive

`φ = Q/3` is **not a loose coincidence; it is as tight as Koide itself.**  From the
PDG-2024 masses (see `study.py`):

| relation        | measured       | target | deviation | significance |
|-----------------|----------------|--------|-----------|--------------|
| Koide Q = 2/3   | 0.66666051     | 2/3    | 6.2×10⁻⁶  | 0.91σ        |
| phase  φ = 2/9  | 0.22222963     | 2/9    | 7.4×10⁻⁶  | 0.89σ        |

Both agree at ~10⁻⁵ and ~0.9σ, *both limited by the same m_τ uncertainty*.  So the
two charged-lepton shape parameters are both functions of one structural number
Q = dimG₂/dimSpin7 = 2/3.  In that precise sense `φ = Q/3` is a **second Koide
relation**.  (Machine-checked in `lean/LeptonPhaseEmpirical.lean`:
`phase_as_tight_as_koide`, `generation_count_pinned`.)

The divisor is **uniquely the number of generations**: among φ = Q/n only n = 3
matches (n = 1,2,4 miss by > 10⁻²).

## 3. Empirical status — the limits (what the law is NOT)

These clarify the law by excluding readings of it (machine-checked in
`lean/PhaseAmbiguity.lean`, `lean/PhaseExclusions.lean`, `lean/BrannenPhase.lean`):

* **Lepton-specific.**  Quarks fit the same circulant form but with unrelated
  phases (down 0.110 ≠ Q_d/3 = 0.244; up 0.074 ≠ Q_u/3 = 0.283).  So the law is a
  property of the *lepton* (color-singlet) sector, not universal.
* **Convention-fixed, not a free observable.**  The masses determine the phase only
  up to the generation S₃ (shift φ→φ+2π/3 and conjugation φ→−φ); the only invariant
  is cos 3φ = cos(2/3) ≈ 0.7859.  "φ = 2/9" is the principal-branch representative
  (a naive fit returns 2/9 + 2π/3 ≈ 2.317).
* **Not a geometric angle.**  2/9 is rational *in radians*, hence not a rational
  multiple of π — so it is not a holonomy / rotation / root angle.
* **Not the Weinberg angle.**  The natural "phase = gauge mixing" identity
  cos 3φ = cos²θ_W = 7/9 is *false*: cos(2/3) > 7/9 (proven).  The gauge 2/9
  (= sin²θ_W, tree/Pati-Salam) shares only the *value* Q/3, not the cosine.
* **Not a naive vacuum.**  The lowest Z₃-invariant phase potentials (cos 3φ, cos 6φ)
  have no extremum at 2/9.

## 4. The central question

Q is independent of φ in the kernel (`BrannenPhase.Q_phase_independent`): nothing in
the Brannen matrix *forces* 3φ = Q.  So `φ = Q/3` is **genuine extra physics** — a
relation between the family CP phase and the mass-spread that holds to 10⁻⁵ but has
no derived mechanism.  It sits exactly where the Koide relation Q = 2/3 sat for
decades: empirically precise, structurally suggestive, mechanism unknown — except
that here the *easy* mechanisms are already excluded (§3), so any real explanation
must come from the family/ideal structure, not geometry, gauge identity, or a simple
potential.  Candidate physical mechanisms are developed in `physical_mechanisms.md`.

## 5. Is it alone?

At the 10⁻⁵ level, the *only* observable precisely at 2/9 is the lepton phase.  The
value Q/3 = 2/9 does recur as the **tree-level / Pati-Salam sin²θ_W** (a structural
claim, not the running M_Z value 0.231), tying mass and gauge sectors at the
*value* level (`lean/TwoNinthsUnification.lean`).  The CKM ~0.22 numbers (Cabibbo,
Wolfenstein λ ≈ 0.225) are percent-level near-misses, not siblings.  So:

* the *precise* (10⁻⁵) law `φ = Q/3` is, so far, **alone** (lepton-only);
* the *value* 2/9 = Q/N_gen is **not** alone — it is a node in the broader v59
  pattern "SM dimensionless numbers = (Lie-dimension ratios)/(generation powers)".

The open question this study exists to resolve: is `φ = Q/3` the lepton instance of a
genuine law `phase = (sector Koide ratio)/(generations)` that the quarks satisfy in a
*branch- or form-corrected* way — or is the lepton case singular?  (Current evidence:
singular, but the quark phases have not been exhausted across branches/orderings.)
