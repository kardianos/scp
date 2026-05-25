# Octomath: the "missing 1/2" is the L-grade complex-structure root-product (Lean-locked)

*2026-05-24.  Pushing connection (3) with "octomath" (the grade/octonionic square-sign law, not
ordinary algebra), per the steer.  **Result: the `1/2` is the root-product of the *complex-structure*
half-element, fixed entirely by the grade's square-sign — and the lepton mass is *proven* to be
that grade (L).** Ground relations machine-checked: `lean/OctoHalf.lean` (axiom-clean, in
`AxiomCheck`).  Computation: `g2_koide_derivation.py` sibling (inline grind).*

## The basic relations (ground, then Lean)

For *any* element `u` in *any* ring — the **half-element law** (`half_element_law`):

> `(1+u)² = 2(1+u) + (u² − 1)`.

Everything about the half-element `P = (1+u)/2` is controlled by **`u²`** — and `u²` is fixed by
the **grade** (`BladeSquareSign`: even `k`-blade squares to `(−1)^{k(k+1)/2}`):

| `u` | grade | `u²` | `P=(1+u)/2` | monic poly | **root-product** |
|---|---|---|---|---|---|
| complex structure | **L = Λ²⊕Λ⁶** | `−1` | `P² = P − ½` | `X² − X + ½` | **½** |
| real structure | **F = Λ⁴** | `+1` | `P² = P` (idempotent) | `X² − X` | **0** |

Verified in the octonions (`e₁²=−1`), in the Cl(7) L-blade (`e₀₁²=−I`) and F-blade (`e₀₁₂₃²=+I`),
and **machine-checked over any field** (`complex_half_field`, `real_half_field`, `root_product_*`).

## The point (why this answers "are you sure it's ½ / what's the missing ½?")

The `½` is **not** an arbitrary ordinary-algebra value, and there *was* a missing piece — it is the
**root-product (constant term) of the complex-structure idempotent-half**, and its value is forced
by `u² = −1`:

- **An F-grade (real, `u²=+1`) mass would give root-product `0`** — a genuine idempotent.
- **The L-grade (complex, `u²=−1`) mass gives root-product `½`** — `P² = P − ½`.

And we **proved** (`LeptonRealityForcing`, `BladeSquareSign`) that the lepton mass / complex
structure lives in **L** (`u²=−1`).  So:

> **The `½` is the signature of the lepton mass being L-grade (a *complex* structure).**  It is the
> octonionic/Clifford square-sign (`u²=−1`) feeding through the universal half-element law into the
> root-product `½`.  This is genuinely "octomath": ordinary real algebra (`u²=+1`) gives `0`; the
> complex-structure grade gives `½`.

This also explains the recurring **`⅓`** from the earlier forcing tests (`11_…md`): those minimized
the *real*-idempotent defect (`ζ²=I`, the `u²=+1` branch), which lands at `⅓`/`0`; the physical
`u²=−1` (L-grade) branch is the one carrying `½`.

## Status (honest)

- **Locked (Lean, axiom-clean):** the half-element law and the complex/real root-product split
  (`½` vs `0`), signature-independent.  Tied to the *proven* `mass ∈ L`.  ⇒ **the `½` is no longer
  arbitrary — it is the L-grade complex-structure root-product.**
- **Still open (the same gap):** this is the *internal* `½` (the complex-structure half-element).
  Transporting it to the *generation* Brannen amplitude `t²` is the generation↔internal map
  (`08_…md`).  What changed: we now know *what the `½` is* (the complex-structure root-product) and
  *why it is `½` not `0`* (L-grade, `u²=−1`, proven) — so the map's job is sharply defined: carry
  the L-grade root-product `½` to `t²`.

## Files
- `lean/OctoHalf.lean` — `half_element_law`, `complex_half`/`real_half` (ring, signature-indep),
  `complex_half_field`/`real_half_field` (the explicit `½`), `root_product_complex`/`real`.
  Built, axiom-clean, imported in `AxiomCheck`.
- Grind: inline in the session (octonion `e₁`, Cl(7) blades).
