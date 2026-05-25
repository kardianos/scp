# Goal (1) build: Koide from G₂-content of the maximally-mixed vacuum — RESULT

*2026-05-24.  The key build (`g2_koide_derivation.py`): construct `G₂=Aut(𝕆)` from scratch,
decompose the lepton sector `L=so(8)` under it, and test whether the maximally-mixed vacuum's
non-G₂ weight is `t²=(D−dimG₂)/D`.  **Outcome: the G₂ side is now rigorous (constructed, not
assumed); the Koide value follows from two clean principles; one bridge (the generation map)
remains, now sharply posed.***

> **⚠ CORRECTION (2026-05-24, `08_generation_map_result.md`): the bridge does NOT close.**  The
> follow-up `G₂×S₃` generation-map build shows the generation symmetry provably **leaves the
> Brannen amplitude `t²` free** (S₃-covariance forces `M` circulant but fixes no `ξ`), and no
> maximal-mixing over the generation structure gives `1/2` (it gives `0` or `1`).  So the `(D−14)/D`
> below is a **dimensional coincidence in the so(8) Clifford grade, NOT a mechanism for the
> generation amplitude** — P2's bridge has no realization.  The "derivation modulo P1+P2" framing
> of this note is **retracted**; what survives is the *arithmetic* + the *built* `dimG₂=14`.  Read
> this doc as "the rigorous G₂ side," and `08_…md` for the (negative) verdict on the mechanism.

## What was constructed and verified (rigorous, from first principles)

1. **Octonions** built from Fano triples; verified **alternative & normed** (a genuine 𝕆).
2. **G₂ = derivation algebra of 𝕆** — `{D∈so(7) : D(uv)=(Du)v+u(Dv)}` — solved as a null space.
   **dim = 14, computed, not assumed.**  So the `14` in `t²=1−14/D` *is* `dim Aut(𝕆)`, grounded.
3. **Decomposition of the lepton sector** `L = Λ²⊕Λ⁶ = so(8)` (D=28, proven `L=skew=so(8)`):
   under G₂, `so(7)=21 = G₂(14) ⊕ 7`, and `so(8)=28 = G₂(14) ⊕ 7 ⊕ 7`, so
   **`L = G₂(14) ⊕ nonG₂(14)` — non-G₂ fraction exactly `14/28 = 1/2`.**
4. **Uniformity (the maximal-mixing signature):** *every* `so(7)` bivector has the **same**
   G₂-fraction `14/21 = 2/3` — no preferred direction.  So G₂ sits as a "diagonal" 14-dim
   subspace, and a maximally-mixed (uniform) vacuum has weight `(D−dimG₂)/D` outside G₂ — not by
   choosing special directions, but as the uniform average.

## The derivation (modulo two principles, both clean)

> **maximal mixing (P1) + G₂ carries no mass-splitting (P2)  ⟹  t² = (D−dimG₂)/D  ⟹  Koide.**

| sector | `D` | `t²=(D−14)/D` | `Q=(1+2t²)/3` |
|---|---|---|---|
| lepton (`L`, D=28=2·dimG₂) | 28 | **1/2** | **2/3** |
| d-quark (`F=Λ⁴`, D=35) | 35 | 3/5 | 11/15 |
| u-quark (`L⊕F`, D=63) | 63 | 7/9 | 23/27 |

- **(P1) maximal mixing** — the vacuum is uniform over the sector.  *Vindicated* by
  `05_v58_vacuum_alignment.md`: the v58 energy is flat on the vacuum manifold, so maximal
  symmetry is the (consistent, primitive) selection principle, not an energy accident.
- **(P2) G₂ is the inert core** — the octonion automorphisms preserve the algebra structure, so
  they carry **no generation-distinguishing (mass-splitting) information**; the mass splitting is
  the complement.  Physically well-motivated (an automorphism cannot tell generations apart), and
  the `14` is now verified as exactly `dim Aut(𝕆)`.

Neither principle is a fit: `14=dimG₂` is computed, `D` are grade dimensions, P1/P2 are
structural.  The lepton lands on the symmetric `1/2` *because* `D_lepton = 28 = 2·dimG₂` — the
sector is exactly double the automorphism core.

## The one remaining bridge (honest — and now sharply posed)

The build rigorously fixes the **G₂ side** (the `14`, the `(D−14)/D` decomposition).  It does
**not** by itself prove that the Brannen amplitude `t²` (the off-diagonal² of the *3×3 generation*
mass matrix) **equals** the L-vacuum's non-G₂ weight.  That identification is (P2) made explicit,
and it requires the **generation map** — how the 3 generations sit relative to the D-dim sector.

This is now sharply posed by the sedenion automorphism **`Aut(𝕊) = G₂ × S₃`** (`sedenion_S3.py`,
`koide_phase_law/`): the octonion `G₂` (inert core) and the **`S₃` = generation permutation** are
*orthogonal factors*.  The Brannen kernel `M=a(I+ξS+ξ̄S²)` is the `S₃`-circulant; its amplitude
`t²` is the `S₃` content, and the claim `t² = (D−dimG₂)/D` is that the `S₃`-amplitude equals the
non-`G₂` weight of the maximally-mixed sedenion vacuum.  **Constructing the `G₂×S₃` factorization
explicitly (the inter-generational map) is the remaining step** — the same one flagged open in
`FORWARD_PROPOSAL_v58_dynamics.md`, now with both factors identified.

## Honest status

- **Rigorous:** `G₂=Aut(𝕆)` has dim 14 (built); `L=so(8)=G₂⊕nonG₂` with non-G₂ fraction `1/2`
  (built); uniform per-direction G₂-fraction `2/3` (built); `t²=(D−14)/D ⇒` the three Koide `Q`s
  (arithmetic).  The `14` is grounded as `dim Aut(𝕆)`, removing the "fitted integer" objection.
- **Principled (not fitted):** P1 (maximal mixing, vindicated by v58 flatness), P2 (G₂ inert).
- **Open:** the explicit `G₂×S₃` generation map proving `t²_{Brannen} = ` non-G₂ weight.  Sharply
  posed; not yet built.

**Net:** this is **not** a from-nothing derivation of Koide, but it is a *real* one **modulo two
clean structural principles**, with the central integer `14=dimG₂` now constructed rather than
asserted and the lepton's `1/2` explained by `D_lepton=2·dimG₂`.  The remaining gap is a single,
sharply-posed map.  Good enough to formalize the solid chain in Lean and wire it into the system.
