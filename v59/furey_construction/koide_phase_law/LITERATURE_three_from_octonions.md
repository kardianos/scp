# Literature: is the "/3" (three generations) emergent from octonion space?

*Search May 2026.  Question: is "the factor 3 / three generations emerges from octonionic
structure" discussed?  Bottom line: **yes, extensively — and the precise finding explains our
holonomy-negative**: the generation `Z₃` is NOT in the octonions proper (`Aut(𝕆)=G₂`), it
emerges in the **sedenion** extension (`Aut(𝕊)=G₂×S₃`).*

## What the literature says

1. **Furey, "Generations: Three Prints, in Colour" (arXiv:1405.4601).**  From the complex
   octonions `ℂ⊗𝕆 ≅ ℂℓ(6)` (64-ℂ-dim), the `SU(3)` color generators partition the algebra
   into "six triplets, six singlets, and their antiparticles," a structure that *mirrors* three
   generations.  **Recognition, not derivation** — the "3" is the color triplication, and it was
   later debated whether these are genuine (inter-mixing) generations.

2. **Gresnigt et al., "Algebraic realisation of three fermion generations with S₃ family … from
   ℂℓ(8)" (arXiv:2407.01580, EPJC 2024) — the most on-point.**  Exactly three generations are
   grounded in the **automorphism group of the sedenions**:
   > `Aut(𝕊) = Aut(𝕆) × S₃`,   with `S₃` "the algebraic source for the existence of exactly
   > three generations."
   The order-3 generator `ψ` (mixing octonion parts with a `√3` factor) **permutes the three
   generations** — this is precisely our Brannen generation `Z₃ ⊂ S₃`.  Crucially, this `S₃` is
   **absent from octonion automorphisms** (`Aut(𝕆)=G₂`); it appears only on extending one
   Cayley–Dickson step beyond the division algebras, to the sedenions `𝕊`.

3. **Dixon** (`T = ℂ⊗ℍ⊗𝕆`, Leech-lattice motivated): three generations encoded in the
   tensor algebra; Furey & Hughes develop the Dixon-algebra route to the SM gauge group.

## Why this matters for us — it explains the holonomy-negative

Our result (`holonomy_result.md`): the generation `Z₃` is **external to the single octonion
ideal** — the color `Z₃` inside `ℂ⊗𝕆` fixes the lepton singlet, so `φ=Q/3` is not an
intra-octonion holonomy.  The literature gives the structural reason:

> the generation `Z₃` is the `S₃` factor of `Aut(𝕊)=G₂×S₃` — it lives in the **sedenion
> extension**, not in `Aut(𝕆)=G₂`.

So "the `/3` is emergent from octo-space" is **true with a precise caveat**: it is not internal
to the octonions; it emerges when octonions are extended to sedenions, exactly the
"expand the algebra" step.  **Concrete next move:** rerun the generation-cycle holonomy test in
the sedenion structure, where the generation `Z₃` is a genuine automorphism `ψ` (with its `√3`)
— there the closed-cycle phase is an honest holonomy and `arg(ξ³)=Q` becomes a well-posed,
non-tautological test.  (Note the `√3` in `ψ`: a `√3`/angle structure at the generation level is
suggestive given that `cos(2/3)` — the genuine invariant — is the target.)

## The "3 real dimensions" point (your note)

Pulling **three real spatial dimensions** out is a *different, well-established* "3":

- `Im(ℍ) = 3` — the three quaternionic imaginaries `(i,j,k)` are the 3 spatial directions
  (real part = time); the division-algebra → spacetime ladder is `ℝ,ℂ,ℍ,𝕆 ↦ 3,4,6,10` dims
  (the SUSY/Kugo–Townsend correspondence): octonions ↦ 10D, quaternions ↦ 6D.
- **This is the same `ℍ` our Koide structure lives in.**  The Brannen ℍ-slice is
  `Cl⁺(3,0) ≅ ℍ` (`FORWARD_PROPOSAL_v58_dynamics.md`), and the v58 bivector (chirality/space)
  sector is `Im(ℍ)` — so the lepton mass quaternion and the 3 emergent spatial dimensions are
  plausibly the **same `ℍ`**.  v58's "emergent effective 3D" program is essentially realising
  this `Im(ℍ)=3`.

So there are (at least) **three different "3"s** in the algebra, and the literature mostly keeps
them separate:

| "3" | source | role |
|---|---|---|
| 3 spatial dims | `Im(ℍ) = 3` (quaternion imaginaries) | v58 emergent space; the Koide ℍ-slice |
| 3 colors | `SU(3) ⊂ G₂ = Aut(𝕆)` | the color triplet (in one octonion ideal) |
| 3 generations | `S₃` of `Aut(𝕊)=G₂×S₃` (sedenions) | the Brannen `Z₃`; the `/3` in `φ=Q/3` |

A **minority view worth flagging** (P. Fré et al., *Pure Spinors → Fermion Physics*,
hep-th/0107158): quaternions "generate the SU(2) of isospin and are at the origin of 3
families" — i.e. a hint that the spatial/isospin `ℍ`-3 and the family-3 might be linked.
Mainstream treatments derive the family-3 from the sedenion `S₃` and the spatial-3 from `Im ℍ`
independently.

## Honest verdict

- **Yes**, "3 from octo-space" is a real, actively-developed idea.  For *generations* the sharp
  result is `Aut(𝕊)=G₂×S₃`: the `/3` is the **sedenion** `S₃`, not internal to the octonions —
  which is exactly why our intra-octonion holonomy test came back trivial.
- Your instinct ("any 3 from octo-space is good news for 3 real dims") is supported by the
  `Im(ℍ)=3` spacetime ladder, and the Koide structure already sits in that very `ℍ`.  Whether the
  **family-3 and the spatial-3 are the same 3** is open and tantalising; the literature leans
  "no" (different sources: `S₃`-sedenion vs `Im ℍ`), with one minority claim of a link.
- **Forward direction:** lift the φ=Q/3 search to the **sedenion** layer, where the generation
  `Z₃` is a genuine automorphism `ψ` (carrying a `√3`), and test the holonomy there.

## Sources
- Furey, *Generations: Three Prints, in Colour* — https://arxiv.org/abs/1405.4601
- Gresnigt et al., *Three fermion generations with S₃ family from ℂℓ(8)* — https://arxiv.org/html/2407.01580
- Furey, *Three generations, two unbroken gauge symmetries, and one 8-dimensional algebra* —
  https://link.springer.com/article/10.1007/JHEP10(2014)046 (and academia.edu mirror; exact
  arXiv id not captured in this search — verify before citing)
- *Why are there three generations of fermions?* — https://arxiv.org/pdf/1412.7658
- Fré et al., *Geometry of Pure Spinors … to Fermion Physics* — https://arxiv.org/pdf/hep-th/0107158
- Baez–Huerta, *Division Algebras and Supersymmetry* — https://arxiv.org/pdf/0909.0551
