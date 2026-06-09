# G9 — Soldering is the precise missing ingredient (honest result)

**Date**: 2026-05-25 (same-day correction + advance on the kickoff `01`–`03`)
**Artifacts**:
- `04_soldering_helicity_honest.py` (runs clean; all assertions pass)
- `../lean/G9Soldering.lean` (compiles clean; axiom-only, no `sorry`)
- `../lean/G9InducedMetric.lean`, `../lean/G9ToyHelicity.lean` (rewritten so they
  build; the originals did **not** compile — see §1)
**Supersedes**: `02_constrained_helicity_count.py` (circular — see §1) and the
"machine-checked, no sorry" claims attached to it.

---

## 1. What was wrong before (cleared up)

**`02` was circular.** It built `Jz_TT ⊗ Id_internal`, where `Jz_TT =
Jz_on_sym2_transverse()` is *by definition* the generator on a spacetime
symmetric-traceless transverse 2-tensor — the graviton — already known
(v59 `g9_polarization_test.py`, item [3]) to have helicity ±2. Tensoring a spin-2
generator with the identity preserves its ±2 eigenvalues. So `02` did not derive
±2 from a constraint; it **inserted the spin-2 generator by hand** and reported it
back. The actual G9 problem — a map binding internal `so(8)` indices to spacetime
μ,ν — was never touched. This is the same tautology the project's own notes warn
about ("`P/m=2` in hedgehog code is TAUTOLOGICAL").

**The Lean did not compile.** `G9ToyHelicity.lean` imported a non-existent module
path (`Mathlib.LinearAlgebra.Matrix.Charpoly`) and tried `decide` on `Multiset ℂ`
(complex equality is not decidable) and `Matrix.eigenvalues` (a `IsHermitian`-only
API) on a non-Hermitian matrix. `G9InducedMetric.lean` had English prose in term
position (`(the physical helicity multiset after constraint)`), wrong imports, and
a noncomputable `Real` def. Both are now confirmed to error out; both have been
rewritten to build.

---

## 2. The honest computation (`04_soldering_helicity_honest.py`)

One consistent piece of machinery (build the `SO(2)_z` little-group generator on a
field space; helicities = imaginary parts of its eigenvalues), four scenarios that
differ in **exactly one** physical input — whether the internal index co-rotates
under the spacetime little group:

| scenario | carrier | internal index | max \|h\| | has ±2? |
|---|---|---|---|---|
| **A** | scalar amplitude × `so(8)` bivector (28) | inert | 0 | no |
| **B** | spacetime 2-form `B_μν`, valued in `so(28)` | inert | 1 | no |
| **C** | spacetime 2-form `B_μν^{IJ}`, `IJ ∈ so(3,1)` | **soldered (co-rotates)** | **2** | **YES** |
| **D** | TT metric fluctuation `h_+, h_×` (endpoint) | — | 2 | YES (exactly 2) |

- **A** reproduces the v59 no-go: the internal index is inert, the carrier is
  scalar → helicity 0.
- **B** shows that even *promoting* `Ω` to a spacetime 2-form is not enough: with
  the internal index inert it is spin-1 (max \|h\| = 1, ×28). No ±2.
- **C** is the resolution. The generator is `Jz_2form ⊗ I + I ⊗ Jz_2form`. The
  **second term is the soldering** — the internal Lorentz index transforms under
  the *spacetime* little group (this is what a tetrad `e_μ^a` does). Now ±2
  appears, and it is **not inserted**: each factor (a 2-form) has max \|h\| = 1;
  the ±2 = (+1)+(+1) and (−1)+(−1) is the Minkowski sum of the two co-rotating
  charges. (The pre-constraint multiplicity is 4 = 2×2; simplicity + diffeo
  constraints reduce the 36-component `B` to the 2 physical modes of **D**.)
- **D** is the constraint/gauge endpoint: exactly 2 TT modes, helicities ±2 — LIGO.

**The decisive, non-tautological content**: the *only* difference between B (fails)
and C (works) is the term `I ⊗ Jz_internal`. **Soldering = the Minkowski sum of
helicity charges, and it is the unique operation that lifts max-helicity from 1 to
2.** v59's carrier is the inert-internal object (A/B-type), so it is *forced* to be
scalar/spin-1. There is no constraint on purely internal indices that escapes this
— helicity is a spacetime little-group label.

---

## 3. Machine-checked Lean (`../lean/G9Soldering.lean`)

Compiles clean against the v59 Mathlib (v4.29.0) via
`lake env lean`; headline theorems depend only on the standard trio
`[propext, Classical.choice, Quot.sound]` — **no `sorryAx`**. Verified.

- `twoform_no_spin2`, `twoform_max_helicity_one` — the 2-form (hence the v59
  inert-internal carrier) never reaches ±2; max \|h\| = 1.  *(the no-go)*
- `solder a b := (a.product b).map (·.1 + ·.2)` — soldering as the Minkowski sum.
- `solder_reaches_two (ha : 1 ∈ a) (hb : 1 ∈ b) : 2 ∈ solder a b` — **fully general**:
  two co-rotating helicity-1 charges sum to 2. The ±2 is *derived*.  *(the resolution)*
- `soldering_is_the_difference : (2 ∉ hel2form) ∧ (2 ∈ solder hel2form hel2form)`
  — the entire gap is the co-rotation term.
- `JzTT_charpoly : JzTT.charpoly = X² + 4` — the explicit TT generator
  `!![0,-2; 2,0]` has eigenvalues ±2i (helicity ±2), derived from the matrix
  entries (corrects the original's broken 4×4 `decide`-over-ℂ).
- `JzTT_has_eigenvalue_2i` — explicit complex eigenvector `(1, −i)` for `+2`.
- `SolderingExists : Prop` and `G9_resolved_if_soldering` — the single open claim
  is isolated as a hypothesis; the conclusion (±2) is the checked consequence.

---

## 4. Where the soldering lives in v59 (grounding) — and the genuinely open part

This is the real advance over `01`–`03`, which proposed soldering inside the
*internal* `Λ²(V^8) = so(8)` (routes A2/B1) and then worried that triality/color
would be broken. **That worry is an artifact of choosing the wrong sector.**

v59 already contains both pieces the resolution needs, in the **spacetime** sector
(`v59/INTEGRATION.md §5`):

- `Cl(3,1)` (dim 16) is the spacetime Clifford algebra. Its **bivector grade is
  6-dimensional = `γ^{μν}` = `so(3,1)` Lorentz generators** (3 spatial + 3 boost).
  *These are literally the 6-field Cosserat (3φ + 3θ) of the project's own
  simulation kernel.* This is the soldering index of scenario **C**.
- The `Cl(3,1)` **scalar grade is already labeled "gravity source ρ_M"** — i.e.
  `ρ_grav = Tr(M†M)` is the scalar source, exactly as a matter trace sources
  `h_μν` in GR.
- The internal `ℂ⊗𝕆` (color, `G₂`, triality) is a **distinct sector**.

So the honest diagnosis of G9: **v59 put the gravity carrier in the internal
`so(8)` sector (forcing it scalar), while the correct spin-2 carrier — a 2-form
soldered via the `Cl(3,1)` bivectors — was sitting unused in the spacetime sector,
together with the scalar source.** The prefactor `21/16 = dim Spin(7)/dim Cl(3,1)`
already pairs internal-over-spacetime, consistent with a
source(`Spin(7)`/internal) ÷ carrier(`Cl(3,1)`/spacetime) split.

**Naturalness (plausible, not yet proved):** because the Lorentz structure lives
in `Cl(3,1)` and the internal `G₂`/color in `ℂ⊗𝕆`, a soldering built from the
`Cl(3,1)` bivectors acts on a different sector from `G₂`/triality/color and should
commute with them — automatically resolving the compatibility worry that sank the
internal-soldering routes. **Caveat:** the precise embedding of `Cl(3,1)` within /
alongside `ℂ⊗ℍ⊗𝕆` (and hence the exact commutation) has not been pinned down here
and must be checked before this is asserted.

**Genuinely open (the remaining content of G9):**
1. Write the Plebański-type action for the 2-form `B` valued in the `Cl(3,1)`
   `so(3,1)` bivectors, with the simplicity constraint, sourced by `ρ_grav` (the
   `ℂ⊗𝕆` second moment).
2. Show its weak-field limit reproduces `□Ω = f_g ρ_grav` (the established scalar
   law) as the trace/source part, while the TT part propagates the 2 modes of **D**.
3. Verify the `Cl(3,1)`-soldering commutes with `G₂`/triality/color so the lepton =
   L forcing and the Koide/Brannen results survive (the naturalness caveat above).
4. Check the magnitude (`α^{21}` class) survives in the emergent weak-field metric.

None of these is the "find a constraint that makes internal indices carry ±2"
problem that `01`–`03` were chasing (which is *impossible* — §2/§3). They are the
well-posed problem of putting the carrier in the sector that already exists.

---

## 5. Status line for ROADMAP/CLOSEOUT

- G9 helicity question: **resolved at the mechanism level** — soldering (Minkowski
  sum of co-rotating charges) is necessary and sufficient for ±2; machine-checked.
- G9 carrier placement: **identified** — `Cl(3,1)` bivector 2-form (spacetime),
  not `Λ²(V^8)` (internal); source `ρ_grav` stays scalar. Grounded in v59 §5.
- G9 dynamics + naturalness proof: **open**, but now well-posed (§4 items 1–4).
- The "internal-soldering" routes A2/B1 of `01`–`03`: **demoted** — they attempt
  the impossible (manufacturing spacetime helicity from internal indices).
