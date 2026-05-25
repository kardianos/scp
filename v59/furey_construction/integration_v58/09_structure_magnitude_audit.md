# Structure/magnitude audit + the missing multivector representation

*2026-05-24.  A conservation-of-structure audit prompted by "we are missing a constraint": list
the physical core quantities and the octo-space structures, map them, read off the residue.  The
residue reveals a clean **discrete-vs-magnitude dichotomy**, and — per the steer — a **methodological
failure**: we have been using octo-space's *symmetry/grade* content (magnitude-blind) and **missing
its multivector representation** (the geometric product / associator / forms), which is where
magnitude lives.  Companion to the bridge map in `08_…md`.*

## The bridge (recap)

The bridge must connect the **algebra** (Cl(7)/𝕆/𝕊, grades, `G₂`, `S₃`, dimensional invariants) to
the **physics** in the Brannen kernel `M = a(I + ξS + ξ̄S²)`: scale `a`, amplitude `t²=|ξ|²`
(Koide `⟺ t²=½`), phase `φ=arg ξ` (`⟺ 2/9`).  The lepton sector splits along one seam: the `/3`
and the circulant **form** are fixed by `S₃`; the **magnitudes** `t², φ` are **free** under `S₃`.

## Side A — physical core quantities (lepton + EW + Higgs)

`D` = discrete (count/rep/quantized); `M` = continuous magnitude.

| # | quantity | type | value |
|---|---|---|---|
| A1 | generation count | D | 3 |
| A2 | electric charges (ν,e,d,u) | D | 0,−1,−⅓,+⅔ |
| A3 | gauge group SU(2)×U(1) | D | — |
| A4 | mass-bearing grade / complex structure | D | `L=Λ²⊕Λ⁶` |
| A5 | `sin²θ_W` | M (rational) | 2/9 |
| A6 | Koide amplitude `t²` | **M** | ½ |
| A7 | Brannen phase `φ` | **M** | 2/9 |
| A8 | `α` | **M** | 1/137 |
| A9 | `v_Higgs/a_lepton²` | **M** | 784 |
| A10 | Higgs quartic `λ` (`m_H/v`) | **M** | untouched |
| A11 | `a_lepton` | **M, dimensionful** | the unit |

## Side B — octo-space structures

`S` = symmetry (group/automorphism, magnitude-blind); `T` = tensor (specific numbers, magnitude-carrying); `#` = count.

| # | structure | type |
|---|---|---|
| B1 | grades Λ⁰,Λ²,Λ⁴,Λ⁶ = 1,21,35,7 | # |
| B2 | `L=so(8)=28`, `F=Λ⁴=35` | # |
| B3 | `Spin(7)=21`, `h∨=5` | S/# |
| B4 | `G₂=Aut(𝕆)=14` | **S** |
| B5 | `SU(3)_color ⊂ G₂` | **S** |
| B6 | `S₃=⟨ψ(Z₃),ε(Z₂)⟩`, `Aut(𝕊)=G₂×S₃` | **S** |
| B7 | Witt/Fock `8=Λ•ℂ³`, N-grading, αᵢ | # + T |
| B8 | complex structure `J∈Λ²` | T (fixed ±) |
| B9 | **octonion product / associator `[x,y,z]`** | **T** |
| B10 | **G₂-invariant 3-form `φ₃`, coassociative 4-form `∗φ`** | **T** |
| B11 | **the norm form** | **T** |
| B12 | `dim Cl(3,1)=16` | # |

## The map (with mechanism status)

| physical | ← octo | status |
|---|---|---|
| A1 gen=3 | B6 ψ | **mechanism ✓** |
| A2 charges | B5 `Q_em=⅓·(color#)` | **mechanism ✓ (discrete!)** |
| A3 gauge grp | B3 Spin(7) decomp | **mechanism ✓** |
| A4 mass∈L | B8 `J∈Λ²` forcing | **mechanism ✓ (proven)** |
| A5 `sin²θ_W=2/9` | B3,B4 Pati-Salam (5,2) | semi-mechanism |
| A6 `t²=½` | B4 `dimG₂/dimSO7` | **numerical only — FAILED** |
| A7 `φ=2/9` | B6 (`/3`) + ? | `/3` ✓; **magnitude FAILED** |
| A8 `α` | B3,B12 | **numerical — FAILED** |
| A9 `v/a²=784` | B2 `dim(L)²` | **numerical — partial/failed** |
| A10,A11 | — | unmapped |

## The residue

- **Unmapped physical:** A6–A11 — and **every one is a magnitude (M)**.  *No discrete quantity is
  left unmapped.*
- **Unused octo-structure:** B9 (associator/non-associativity), B10 (the invariant forms — used
  only to *define* G₂, never as dynamical objects), B11 (norm form), and the `ε` of `S₃` — and
  **every one is a tensor (T), not a symmetry.**

## Diagnosis 1 — the discrete/magnitude dichotomy

> **Everything octo-space has fixed (with mechanism) is discrete; everything it has failed to fix
> is a continuous magnitude.**  Charges, generation count, gauge group, mass grade, the structural
> *ratios* — all **discrete/rep-theoretic** — are mapped.  All **magnitudes** (`t²`, `φ`, `α`) are
> unmapped or numerical-only.

This is not accidental: we used only the **symmetry** content of octo-space (`G₂, S₃, Spin(7),
SU(3)`), and **symmetries are magnitude-blind by construction** — they fix *which* and *how many*,
never *how much*.  Every failed attempt to fix `t²`/`φ` (holonomy, maximal-mixing, G₂-inert,
energy-min) reached for symmetry to do a job symmetry cannot do.  v58 failed the same way: its
energy was *flat* on the vacuum (magnitude-blind on alignment).

## Diagnosis 2 (the deeper one) — a failure of method: the MISSING MULTIVECTOR REPRESENTATION

The audit's two residues **type-match**: the leftover *physical* quantities are all magnitudes;
the leftover *octo* structures are all magnitude-carrying tensors.  This is the steer's point:

> **We have been looking at octo-space the wrong way.**  We represented it by its **symmetry
> groups** (`G₂, S₃`) and its **grade-dimensions** (1,21,35,7) — i.e. by *what is invariant* — and
> that view is **magnitude-blind**.  We have **not** used the **multivector representation within
> octo-space**: the geometric/Clifford **product** itself, with its grade-mixing and the
> magnitudes relating the grades, and the **associator** (the non-associativity = the genuinely
> octonionic content).  The symmetry group `G₂` is only the *automorphisms* of the product — a
> shadow.  Working with the shadow guaranteed we could never fix a magnitude.

The magnitude-carrying content of octo-space **is** the multivector structure (B9–B11): the
geometric product of two vectors splits into scalar + bivector + … with definite *relative
magnitudes*; the associator `[x,y,z]` is a specific multivector (not a symmetry); the invariant
forms are tensors with numerical values.  These are exactly the inputs that could fix `t²`, `φ`.

**Concretely, the methodological correction:** treat the order parameter `ξ` / the mass operator
`M` not as an abstract scalar coupling on an abstract 3-generation index, but as a genuine
**multivector** in the geometric algebra, and seek a **geometric-product constraint** — an
idempotent (`M²=M`-type), a norm/reversion condition (`MM̃=…`), or an associator/3-form condition —
that fixes its grade magnitudes.  That is a *magnitude* constraint of the right type, and it is the
class of structure we have never engaged.  (It also re-legitimises v58 as the dynamical home: v58
*is* a multivector field theory — the flatness we found came from projecting to symmetry-invariant
scalars `⟨M²⟩₀,₂`, which discards exactly the multivector magnitude content.)

## Forward direction (the first positive lead the audit produces)

1. **Build the lepton order parameter as a full multivector** in Cl(7) (all grades, the geometric
   product), not a scalar coupling on a generation index.
2. **Impose a geometric-product constraint** native to the multivector view (idempotent / norm /
   reversion / associator), and compute whether it **fixes `t²=½` and/or `φ=2/9`** — a magnitude
   from a magnitude-carrying structure.
3. If yes: the missing constraint is found, and it is the multivector product (B9–B11), not any
   symmetry.  If no: we will have tested the *one* remaining class of octo-structure and can say
   the magnitudes are genuinely external (the deflationary verdict), now on firm ground.

**Net:** the audit converts "we are missing a constraint" into a *typed* statement — **the missing
constraint is magnitude-carrying, and octo-space's only unused magnitude-carrying content is its
multivector representation (geometric product / associator / forms)**.  Every prior attempt used
the magnitude-blind symmetry view; the multivector representation is the untried, type-correct
direction.
