# Combining v58 dynamics with v59 Cl(7): the F-grade channel lifts the flatness (a candidate dynamics)

*2026-05-24.  Per the steer "combine with the v58 equation to get dynamics."  Key finding: lifting
the v58 field equation from Cl(3,0) to v59's **Cl(7)** supplies a channel v58 structurally lacked —
the **grade-4 (F = Λ⁴) sector** — and `⟨M²⟩₄` is **non-flat** on the lepton grade `L`, exactly where
v58's `{0,2}` projection was flat.  So the combined dynamics **can select the vacuum alignment**, and
the selecting term physically **couples the lepton (L) mass to the quark (F) sector**.  This is the
first concrete vehicle for the "missing physical relation."*

> **⚠ CORRECTION (2026-05-24, `20_constrained_ofe_result.md`): the F-channel is INERT on the lepton.**
> The non-flat `⟨M²⟩₄` here uses a *generic 28-dim L-config*.  The actual lepton is the **color
> singlet**, confined to the 4-dim `ℍ`-slice `{1,e₀₁,e₀₂,e₁₂}` (closed `≅ℍ`, grades `{0,2}`), where
> the grade-4 (F) channel is **identically zero** (verified).  So lifting v58 to Cl(7) adds **no**
> working channel for the singlet lepton; the OFE on `ℍ` is the XiVacuum (`{0,2}`), giving `|ξ|²=½`
> only as the Mexican-hat parameter (input) and the phase as a free Goldstone.  Read this doc as "the
> F-channel exists in Cl(7) but is inert on the color-singlet lepton."  The combined-dynamics line is
> therefore **negative for fixing the lepton couplings** (`20_…md`).

## Why v58 alone was flat — and why Cl(7) fixes it

The v58 equation `⟨DΩ + λΩ² + μ⟨Ω,Ω⟩⟩_{0,2} = …` projects to grades `{0,2}` **because that is the
entire even subalgebra of Cl(3,0)** (`Cl⁺(3,0)=ℍ`, grades 0,2 only).  On the lepton grade
`L=skew`, `M²` is symmetric ⇒ `⟨M²⟩₂ ≡ 0` (computed) ⇒ the alignment is a flat modulus (`05_…md`).

But the lepton actually lives in **Cl(7)**, whose even subalgebra is `{0,2,4,6}`.  Lifting the v58
equation to Cl(7) **naturally includes grades 4 and 6** — and:

> `M ∈ L` (skew) ⇒ `M²` symmetric `= grade-0 ⊕ grade-4`.  `⟨M²⟩₂ ≡ 0` (flat), but
> **`⟨M²⟩₄ ≠ 0` and configuration-dependent** — verified: democratic `‖⟨M²⟩₄‖²=5.71`, random `≈5`,
> pure complex structure `=0`.

So the alignment information v58 couldn't see is sitting in **grade 4 = F = Λ⁴**, the channel that
*does not exist* in Cl(3,0).  **The flatness was an artifact of the small algebra.**  In Cl(7) the
`{0,2,4,6}` (full even) dynamics has a non-flat, alignment-dependent grade-4 term.

## The physical relation it supplies: lepton(L) ↔ quark(F) coupling

The grade-4 part of the *lepton* (`L`) mass operator's square lands in **F = Λ⁴ = the quark
grade**.  So a v58-type self-interaction in Cl(7), `⟨Ω²⟩₄`, **couples the lepton mass alignment to
the quark (F) sector** — a genuine **dynamical relation between the lepton and quark sectors**, of
exactly the *physical* (not merely algebraic) kind the steer suspected was missing.  Schematically
the new energy term is `E_F ∝ ‖⟨M²⟩₄‖²` (the F-channel self-energy), absent in v58, non-flat on `L`.

## The combined object (needs a name)

This is **v58's dynamical form** (retarded multivector field equation, quadratic self-interaction
`λΩ²+μ⟨Ω,Ω⟩`, Higgs-like vacuum manifold `MM̃=v²`) **on v59's kinematics** (the Furey `Cl(7)_even`
fermion states, `L`/`F` grades), projected to the **full even subalgebra `{0,2,4,6}`** (not just
`{0,2}`).  Proposed names (pick one):
- **OFE — the Octonionic Field Equation** (the dynamical law on `ℂ⊗𝕆 ≅ Cl(7)_even`);
- **the Cl(7) multivector dynamics**;
- **the L–F dynamics** (named for its signature: the grade-4 lepton↔quark coupling).

## Honest status — necessary, not yet sufficient

- **Established:** the combination is **non-flat** — the grade-4 (F) channel, present in Cl(7) and
  absent in Cl(3,0), is alignment-dependent on `L`.  This removes the *structural* obstruction
  (`05_…md`'s flatness) and gives a concrete coupling-selecting term.  *This is the prerequisite the
  whole "missing dynamics" analysis said was needed.*
- **Not yet shown:** that minimising/extremising the F-channel energy on the `L`-vacuum manifold
  **selects the Koide point** (`t²=½`, `φ=2/9`).  Two sub-questions: (i) what is the actual energy
  functional (sign of the F-term, the `λ`, the `DΩ` kinetic part); (ii) does its extremum, mapped to
  the Brannen amplitude/phase (the generation↔internal map, `14_…md`), land on Koide?  Note the
  F-content is `0` for a pure complex structure and maximal for the democratic config — so the
  extremum is a real, computable selection, not flat.

## Net

**Combining v58 with v59's Cl(7) gives, for the first time, a non-flat coupling-selecting channel —
the grade-4 (F) sector that Cl(3,0) lacked — and it physically couples the lepton mass to the quark
sector.**  This is the concrete candidate for the missing physical (dynamical) relation.  The next
step is to write the combined energy functional on `Cl(7)_even` (grades `{0,2,4,6}`) and test whether
its vacuum extremum selects the Koide couplings — the first time the dynamics could *fix* `t²`, `φ`
rather than leave them free.
