# Scrutiny of "P2 — G₂ inert": stress-testing `t² = (D − dimG₂)/D`

*2026-05-24.  Adversarial audit of the assumption behind
`06_koide_from_G2_maximal_mixing.md` / `g2_koide_derivation.py`.  Computations in
`07_g2_inert_scrutiny.py` (run it: `python3 07_g2_inert_scrutiny.py`); it reuses the
octonion + G₂ construction from `g2_koide_derivation.py` and adds the Λ^k branchings,
the spinor branching, and the commutant computation.  All rep-theory facts are verified
numerically against textbook G₂ branchings.*

## What was tested

The claim under test (P2): **G₂ = Aut(𝕆) is an "inert core" carrying no
mass-splitting information, so the Koide amplitude is the maximally-mixed weight
*outside* G₂: `t² = (D − dimG₂)/D = (D − 14)/D`** — giving lepton 1/2, d-quark 3/5,
u-quark 7/9, hence Koide Q = 2/3, 11/15, 23/27.

I reproduced `g2_koide_derivation.py` (G₂ built as the 14-dim derivation algebra of a
verified 𝕆; so(8) = 14⊕7⊕7) and then ran four new computations.

---

## (1) Is the "−14" the removal of a genuine G₂-subspace in every sector?

**Verified G₂-content of each grade `Λ^k ℝ⁷`** (invariant counts and commutant
dimensions, both stable rank computations; they match the textbook branchings exactly):

| grade | dim | #G₂-invariants | #irreps | branching |
|---|---|---|---|---|
| Λ¹ | 7 | 0 | 1 | **7** |
| Λ² | 21 | 0 | 2 | **7 ⊕ 14** |
| Λ⁴ | 35 | 1 | 3 | **1 ⊕ 7 ⊕ 27** |
| Λ⁶ | 7 | 0 | 1 | **7** |

Now the per-sector "−14":

- **Lepton (D = 28).**  In `g2_koide_derivation.py` the lepton sector is `so(8) = Λ²ℝ⁸`,
  and `so(8) = 14 ⊕ 7 ⊕ 7` with the **14 the G₂ adjoint subalgebra**.  So `28 − 14`
  *does* remove a genuine G₂-invariant subspace.  (Caveat below: this `28` is `so(8)`,
  *not* the `Λ²ℝ⁷ ⊕ Λ⁶ℝ⁷ = 21 + 7` that `Predictions.lean` Tier 7 names as the lepton
  ambient.  I checked these two "28"s are at least **G₂-isomorphic** — both decompose as
  `14 ⊕ 7 ⊕ 7` — so the lepton's `1/2` is robust to which one you mean.  Still, the two
  documents name *different vector spaces* by the same dimension, which should be
  reconciled.)

- **d-quark (D = 35 = Λ⁴ℝ⁷).**  `Λ⁴ = 1 ⊕ 7 ⊕ 27`.  **No subset of the constituents
  `{1,7,27}` sums to 14.**  The G₂ adjoint (14) is **not even a constituent of Λ⁴**.  So
  there is **no 14-dimensional G₂-invariant subspace to remove** — `35 − 14 = 21` is bare
  arithmetic, *not* "the sector minus its automorphism core."  The honest reading of the
  d-quark split would be `1 ⊕ 27` (the 28 G₂-trivial-plus-big part) vs. the `7`, or
  `1` vs. `7 ⊕ 27`, etc. — none of which is "remove G₂."

- **u-quark (D = 63 = Λ²⊕Λ⁴⊕Λ⁶).**  Branching `(7⊕14)⊕(1⊕7⊕27)⊕7 = 1 ⊕ 14 ⊕ 7⊕7⊕7⊕7 ⊕ 27`.
  Here a **14 (the adjoint, living inside the Λ² piece) IS present**, so `63 − 14` *can*
  be read as removing it.  But it is exactly the copy that the d-quark sector lacked
  entirely — so the "subtract the G₂ adjoint" story is available for the u-quark but
  was unavailable for the d-quark, despite both using the identical formula `(D−14)/D`.

**Verdict (1): the uniform "−14 = remove dim G₂" reading holds ONLY for the lepton
(and, weakly, the u-quark).  For the d-quark it is provably NOT a G₂-subspace removal —
the formula `(35−14)/35` is arithmetic that fits, with no G₂-geometric content.**  This
is the single most damaging finding: the SAME integer 14 is given a G₂-subspace meaning
in one sector where it genuinely is one, and applied as bare subtraction in another where
it demonstrably is not.

---

## (2) Does G₂ genuinely fix the lepton (the "inert" claim)?

Built the 8-dim spinor of Spin(7) from explicit Cl(7) gammas (Clifford relations
verified), expressed the 14 G₂ generators in the spinor rep, and counted invariants:

- **8-spinor under G₂ has exactly 1 invariant** (so `8 = 1 ⊕ 7`, as expected), and G₂
  annihilates that singlet vector to `max‖X·v‖ = 3.4×10⁻¹⁶`.  **The G₂ singlet is
  genuinely inert.**  This corroborates `ColorSU3.lean` (color Z₃ fixes the lepton) and
  extends it to the *full* G₂: the singlet is fixed by all 14 generators, not just Z₃.

- **But there are TWO leptons (N=0, 3), and only ONE G₂-singlet.**  The second lepton sits
  in the **7** as the SU(3)-singlet (the 7 → 1⊕3⊕3̄ under SU(3) ⊂ G₂).  I verified the
  **7 is G₂-irreducible** (0 invariants, commutant dim 1).  So G₂ acts *non-trivially* on
  the subspace that hosts the second lepton; only the SU(3) subgroup fixes it.  "G₂ fixes
  the lepton" is therefore true for *one* of the two lepton states and false for the other
  at the full-G₂ level.

**Verdict (2): G₂ truly fixes the lepton G₂-singlet (inert, confirmed), but this is a
statement about ONE singlet direction, not about the lepton generation structure or the
second lepton.  "G₂ carries no info that could split the leptons" does NOT follow from
fixing the singlet — see (3).**

---

## (3) The crux: does "G₂-invariant ⟹ no mass-splitting" hold?

This is the load-bearing logical step, and it **fails as stated**.  The argument in P2 is
"automorphisms can't tell generations apart, so the splitting lives in the complement."
Representation theory says otherwise:

- A mass operator that *respects* the inert core means one that **commutes with G₂**.  By
  **Schur's lemma**, such an operator is a *scalar on each G₂-irrep* — but it may take a
  **different scalar on each distinct irrep**.

- I computed the commutant dimension of G₂ on each sector grade:
  `Λ² → 2`, `Λ⁴ → 3`, `Λ⁶ → 1`.  These count the **free parameters of a G₂-commuting
  operator** (= number of distinct irreps, since these are multiplicity-free).

- **Explicit demonstration:** I constructed a G₂-commuting operator on `Λ⁴` and found its
  eigenvalue multiset has **3 distinct values** — i.e. it assigns *different* masses to
  the `1`, `7`, `27` blocks.  It is manifestly **not** a multiple of the identity.

So a "G₂-inert" (G₂-respecting) mass operator **can** split a sector into its irrep
blocks; it only fails to split states *within a single irrep*.  The three generations are
**not** a single G₂-irrep (and in the Fock picture the generation structure is the
external `S₃`, orthogonal to G₂), so **generation/mass splitting is not forbidden by
G₂-inertness.**

**Verdict (3): "G₂-invariant ⟹ no mass-splitting" is FALSE as a blanket claim.  It
forbids only INTRA-irrep splitting.  The crucial step — that the Brannen 3×3 off-diagonal
amplitude `t²` equals the maximally-mixed non-G₂ *weight* of the sector — is an ADDED
identification, not a consequence of G₂ being an automorphism group.  This is exactly the
"generation map" gap that `06_…md` itself flags as open; my computation shows the gap is
real and not closable by representation theory alone (the inert core leaves a nonzero
commutant that can carry splitting).**

---

## (4) Is "14" privileged as dim G₂?

Three exact readings of 14 collide:
`14 = dim G₂ = dim Aut(𝕆)` ; `14 = 2·7 = 2·dim(Im𝕆) = dim(7⊕7)` ;
`14 = 21 − 7 = dim so(7) − dim(vector 7)`.

Sector-by-sector, which reading carries G₂-subspace meaning:

- **lepton** (`so(8) = 14⊕7⊕7`): all three coincide; the non-G₂ complement *is* `7⊕7 = 2·7`,
  and it *is* the orthogonal complement of the G₂ adjoint.  Genuinely G₂.
- **d-quark** (`Λ⁴ = 1⊕7⊕27`): **neither** `14 = dimG₂` nor `2·7` is a sub-constituent;
  `35 − 14 = 21` matches no block sum.  The only invariant content of "14" here is the
  bare integer.
- **u-quark** (`1⊕14⊕7⊕7⊕7⊕7⊕27`): the adjoint-14 is present (favoring `dimG₂`), but
  `2·7` is *also* realizable (two of the four 7's), so the reading is ambiguous.

**Verdict (4): the G₂ interpretation is privileged ONLY for the lepton.  Everywhere else,
"14" is an integer that happens to equal both `dim G₂` and `2·dim Im𝕆` (and `21−7`).  The
claim that G₂ specifically (rather than the integer 14, or `2·dim Im𝕆`) is the origin is
therefore not supported by the quark sectors; it is supported by the lepton sector, where
the geometry genuinely picks out the adjoint.**

---

## Overall verdict on P2

**P2 is, in its strong global form, a post-hoc fit; in its lepton-only form it is a
well-motivated principle but not a theorem.**  Concretely:

- The **lepton** case is the strong point.  `dim G₂ = 14` is *constructed*, not asserted;
  `so(8) = 14 ⊕ 7 ⊕ 7` is genuine; G₂ fixes the lepton singlet (verified to 1e-16); and
  the non-G₂ fraction is exactly `1/2 ` *because* `D = 28 = 2·dim G₂`.  For the lepton, P2
  is a **well-motivated principle (b)** — though still resting on the unproven
  identification "Brannen `t²` = non-G₂ weight" (the open generation map), so not yet a
  theorem.

- The **quark** cases are the weak point and reveal the formula as **arithmetic (c)** for
  d-quarks: there is **no 14-dim G₂-invariant subspace in Λ⁴** to "remove," so `(35−14)/35`
  has the form `(D − dimG₂)/D` only as numerology.  The uniform application of the *same*
  "subtract dim G₂" across three sectors is not backed by a uniform G₂ decomposition —
  it works geometrically for the lepton, ambiguously for the u-quark, and not at all for
  the d-quark.

- The **central logical claim** "G₂-invariant ⟹ no mass-splitting" is **false** as stated
  (Schur permits inter-irrep splitting; demonstrated explicitly on Λ⁴ with a 3-distinct-
  eigenvalue G₂-commuting operator).  So even granting maximal mixing (P1), P2 does not
  *force* `t²` to be the non-G₂ weight; that is a separate, unproven identification.

### What would upgrade P2 from (c)/(b) to a theorem

1. **A uniform G₂-geometric reading of "−14" valid in ALL three sectors.**  As it stands,
   the d-quark falsifies the "remove the G₂ core" reading.  Either (a) reformulate the
   quark denominators so the subtracted 14 is a real G₂-invariant subspace in each, or (b)
   admit that the quark formula is a *fit* and only the lepton is derived.

2. **The generation map** (`Aut(𝕊) = G₂ × S₃`): an explicit proof that the Brannen 3×3
   off-diagonal² `t²` *equals* the maximally-mixed non-G₂ weight.  My (3) shows this cannot
   come from G₂-inertness alone — the G₂-commutant is nonzero and can itself carry the
   splitting — so the proof must come from the `S₃`/maximal-mixing side, and must rule out
   the inter-irrep splitting that G₂-inertness permits.

3. **Reconcile the two "28"s** (`so(8) = Λ²ℝ⁸` in `06_…md` vs `Λ²ℝ⁷ ⊕ Λ⁶ℝ⁷` in
   `Predictions.lean` Tier 7).  They are G₂-isomorphic (both `14⊕7⊕7`), so the lepton
   number survives, but the documents should name the same space, and the u-quark `63`
   (which mixes `Λ²ℝ⁷`, `Λ⁴ℝ⁷`, `Λ⁶ℝ⁷` of ℝ⁷, not ℝ⁸) cannot then borrow the so(8)
   reading of the 14.

**Bottom line: P2 is a well-motivated principle for the LEPTON (where every piece checks
out except the still-open generation map), and a post-hoc arithmetic fit for the QUARKS
(where the subtracted 14 is provably not a G₂ subspace).  The blanket statement "G₂ carries
no mass-splitting" is not a theorem and is in fact false for inter-irrep splitting.**
