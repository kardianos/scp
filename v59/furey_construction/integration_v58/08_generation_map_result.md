# Goal (1), the G₂×S₃ generation map: RESULT — the bridge does NOT close (negative)

*2026-05-24.  Built the `Aut(𝕊)=G₂×S₃` generation map to test whether the maximally-mixed
vacuum's so(8) non-G₂ weight `(D−dimG₂)/D = 1/2` is realized as the Brannen amplitude `t²`.
**Outcome: it is not.** The generation symmetry provably leaves `t²` free, and no maximal-mixing
over the generation structure yields `1/2`.  So Koide `2/3` is **not** derived — `t²=1/2` is a
residual coupling magnitude, exactly parallel to the phase `φ=2/9`.  This **walks back the
"derivation modulo P1+P2" claim of `06_…md`.***

## What was built and verified

- **Sedenions `𝕊`** (16-dim Cayley–Dickson double of 𝕆): built, verified (identity, `eᵢ²=−1`,
  zero divisors present ⇒ genuine non-division sedenions).
- **`Aut(𝕊) = G₂ × S₃`** with the explicit `S₃ = ⟨ε,ψ⟩` (Gresnigt et al.; `sedenion_s3.py`):
  `ε²=ψ³=I`, `εψ=ψ²ε`, both genuine automorphisms.  **`ψ` = the Z₃ generation cycle** (a 120°
  rotation of the two Im-octonion blocks) — this *is* the Brannen shift `S`.  The `/3` is
  structural (confirmed, as in `sedenion_S3_analysis.md`).

## The decisive test: does maximal mixing fix `t² = 1/2`?  No.

The mass operator is the ψ-covariant circulant `M = a(I + ξψ + ξ̄ψ²)`; `t²=|ξ|²`, `Q=(1+2t²)/3`.
The candidate "maximal-mixing" notions over the **generation** structure give:

| maximal-mixing notion | `t²` | `Q` |
|---|---|---|
| uniform over the Z₃ group `{I,ψ,ψ²}` (democratic circulant) | **1** | 1 |
| all weight on `I` (degenerate spectrum) | **0** | 1/3 |
| maximally-mixed state `ρ=I/3` on the 3 generations | **0** | 1/3 |

**None give `1/2`.**  The Koide value sits *between* these extremes, at no symmetric point of the
generation structure.  And — the root cause — **S₃-covariance forces `M` circulant but leaves `ξ`
(hence `t²`) completely free** ("any `ξ` allowed," `sedenion_S3_analysis.md`).  The amplitude is a
*coupling*, not fixed by the generation symmetry — **exactly** as the phase `φ` is not.

## Why the `(D−14)/D` dimensional fact does not transfer

The `(D−dimG₂)/D = 1/2` of `06_…md` is a fact about the **so(8) Clifford grade** `L` (D=28): under
G₂, `so(8)=28 = 14 ⊕ 7 ⊕ 7`, non-G₂ fraction `14/28=1/2`.  But:

- the **sedenion** fractions are different: `Im(𝕊)=15 = 7⊕1⊕7` under G₂, non-singlet fraction
  `14/15`, **not** `1/2`;
- the **Brannen `t²`** is the off-diagonal of the *generation* circulant, a different object again;
- there is **no symmetry or maximal-mixing map** carrying the so(8) `1/2` to the generation
  amplitude — the three structures (so(8) grade, sedenion, generation circulant) do not align, and
  the one bridge between sector and generations (the `S₃`) leaves the amplitude free.

So the elegant "`t² =` non-G₂ weight" identification (P2 made concrete) **has no realization**:
`(D−14)/D` is a dimensional coincidence in the Clifford grade, not the Brannen amplitude.

## Honest reconciliation with `06_…md`

`06_…md` claimed "Koide derived modulo (P1) maximal mixing + (P2) G₂ inert."  This build tests P2's
*bridge* — that the maximally-mixed non-G₂ weight **is** the Brannen `t²` — and finds **it does not
hold**.  What survives from `06`:

- **Solid:** `dimG₂=14` built from scratch (octonion derivations); `L=so(8)=14⊕14` under G₂; the
  *arithmetic* `(D−14)/D ⇒` Koide ratios.  These are real and remain.
- **Retracted:** the *mechanism* claim.  `(D−14)/D` is **not** realized as the generation amplitude
  by maximal mixing — the generation symmetry leaves `t²` free.  So `06`'s "derivation modulo two
  principles" overstated it: P2 has no bridge to the observable amplitude.

## Net: Koide 2/3 is a residual coupling, parallel to the phase

The result unifies the two pieces of the lepton mystery — but in the **negative** direction:

> **Both the Koide amplitude `t²=1/2` and the phase `φ=2/9` are coupling magnitudes that the
> `Aut(𝕊)=G₂×S₃` structure provably leaves free.**  The generation symmetry fixes the *form*
> (circulant) and the *count* (`/3`); it fixes **neither** the amplitude nor the phase.  `t²=1/2`
> matches the structural `(D−dimG₂)/D`, and `φ=2/9` matches `Q/3`, but **neither match has a
> mechanism** — they are precise (`10⁻⁵`) residual couplings, orthogonal to the symmetry.

This returns Koide to "Kepler, not Newton," now with a *specific, proven* reason: the generation
couplings are free under the full automorphism group.  The maximal-symmetry principle (vindicated
as the vacuum *selector* by the v58 flatness, `05_…md`) does **not** reach down to fix the
generation couplings.  Deriving `t²=1/2` would require a dynamical input *beyond* the `G₂×S₃`
symmetry — the same wall the phase `φ=2/9` hit.

**Recommendation:** record Koide `2/3` honestly as a residual coupling (with `φ=2/9`), demote the
`06` mechanism claim to "dimensional coincidence `(D−14)/D`, no generation-level mechanism," and
keep only the rigorous parts (`dimG₂=14` built; the so(8) decomposition; the arithmetic).
