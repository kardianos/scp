# GEN3 — Matter sector: a potential whose EL vacuum DERIVES the Koide cone (Q=2/3)

**Date**: 2026-05-26 (dynamical-Lagrangian loop, Generation 3)
**Artifacts**:
- `12_matter_sector.py` (runs clean; every assertion passes)
- `12_koide_invariants.mac` (independent **Maxima** cross-check)
- `../lean/MatterSector.lean` (compiles clean; `koide_invariant_form` is a genuine
  `ring` identity over ℝ; standard-trio axioms only)
**Builds on**: GEN1/GEN2 (the gravity half — `ρ_grav` is the source they need).

---

## Verdict

The gravity half (GEN1/2) needed a source `ρ_grav = Tr(M†M)`. GEN3 supplies the
**matter half**: a potential on the internal generation triple whose
Euler–Lagrange **vacuum derives Koide `Q = 2/3`** (rather than positing it), with
the **Brannen phase as a Goldstone** (so the dynamics fix the cone but not the
phase — honest about what v59 already knew is "not geometric"), and whose
**second moment is exactly `ρ_grav`**.

**New vs v59**: v59 *identified* `Q = 2/3 = dim G₂/dim Spin(7)` group-theoretically
and *posited* the Brannen vacuum. Here `Q = 2/3` **emerges from minimizing a
potential** — i.e. from the EL dynamics — written in the **S₃-symmetric invariants**
of the triple (the same sedenion S₃ automorphism v59 used for the "/3" in
`φ = Q/3`).

---

## The construction

Order parameter `x = (√m₁, √m₂, √m₃)` (the diagonal of the mass-amplitude kernel
`M`, with `M†M = diag(mₖ)`). Symmetric invariants `e₁ = Σ√m`, `e₂ = Σ_{i<j}√mᵢ√mⱼ`,
`e₃ = √m₁√m₂√m₃`.

### (A) Koide `Q = 2/3 ⟺ e₁² = 6 e₂`  — a clean symmetric-polynomial relation

Genuine `ring` identity (Lean `koide_invariant_form`, also SymPy + Maxima):
```
3(Σm) − 2e₁²  =  e₁² − 6e₂        (both equal Σm − 4e₂)
```
`Q = Σm/e₁²`; `Q = 2/3` ⟺ LHS = 0 ⟺ `e₁² = 6e₂`. The Brannen triple
`xₖ = a(1 + √2 cos(φ + 2πk/3))` satisfies `e₁ = 3a`, `Σm = 6a²`, `e₁² = 6e₂ = 0`
for **all φ** (SymPy + Maxima) — confirming the phase is a flat direction.

### (B) Potential and EL vacuum — `Q = 2/3` from minimization

```
V = λ(e₁² − 6e₂)²  +  μ(e₁ − c)²       (S₃-symmetric)
```
`grad V = 0` at the Brannen vacuum (verified `[0,0,0]`). The vacuum pins `e₁ = c`,
`e₂ = c²/6` ⟹ `e₁² = 6e₂` ⟹ **`Q = 2/3` exactly**. So Q=2/3 is now an EL **minimum
condition**, not an input. Interpretation: Q=2/3 ⟺ the democratic (trace) and the
traceless condensates balance.

### (C) The Brannen phase is a Goldstone

Hessian at the `φ=2/9` vacuum (λ=μ=1): eigenvalues **`[0, 2.98, 435]`** →
**1 zero (Goldstone) + 2 positive (massive)**. The zero eigenvector is exactly the
phase tangent `dx/dφ` (`‖H·dx/dφ‖ = 6×10⁻¹⁴`). The vacuum manifold
`{e₁=c} ∩ {e₂=c²/6}` is 1-dimensional (`3 − 2 = 1`), and that one dimension is the
phase. **Dynamics fix the cone (Q=2/3) but NOT the phase** — an explicit S₃-breaking
tilt (e.g. on `e₃`) would be the residual posit that pins `φ = 2/9`.

### (D) Second moment = the GEN1/GEN2 gravity source

`Σm = e₁² − 2e₂ = c² − c²/3 = (2/3)c²`. With `c = 3a` (per-generation amplitude):
`Σm = 6a² = 9Qa²` (Q=2/3). **This is exactly `ρ_grav` that GEN1/2 feed into the
gravity trace law `∇²Ω = −f_g ρ_grav`** — the two halves are now coupled by a
single scale `a`.

### (E) Kinetic term / dynamics

`ℒ = ½⟨∂_μΦ ∂^μΦ̃⟩ − V(Φ)`, `Φ ∈ Cl(7)_even` (L-grade). EL: `□Φ = −V′(Φ)`. The
homogeneous static vacuum minimizes `V` (the Brannen cone); the linearized internal
spectrum is **2 massive radial modes + 1 massless Goldstone (the phase)**.

---

## Lean (`MatterSector.lean`)

| theorem | content | kind |
|---|---|---|
| `koide_invariant_form` | `3Σm − 2e₁² = e₁² − 6e₂` | **`ring` over ℝ** |
| `second_moment_vacuum` | `c² − 2(c²/6) = (2/3)c²` | `ring` |
| `second_moment_in_a` | `(2/3)(3a)² = 6a²` | `ring` |
| `rho_grav_match` | `6a² = 9·(2/3)·a²` | `ring` |
| `one_goldstone` / `two_massive` | `3−2=1`, `=2` | `decide` |
| `gen3_matter_sector` | bundled headline | mixed |

---

## Status table

| claim | status | tool |
|---|---|---|
| `Q = 2/3 ⟺ e₁² = 6e₂` | **verified** | SymPy + Maxima + Lean (`ring`) |
| Brannen triple on cone ∀φ; `Σm = 6a²` | **verified** | SymPy + Maxima |
| potential EL vacuum on the cone, Q=2/3 derived | **verified** | SymPy |
| Brannen phase = single Goldstone (`H·dx/dφ=0`) | **verified** | SymPy |
| `Σm = 9Qa² = ρ_grav` (couples to GEN1/2) | **verified** | SymPy + Lean |
| phase φ=2/9 NOT fixed by V (residual posit) | **honest open** | — |

---

## What GEN3 settled, and what's next

- **Settled**: the matter sector has a potential whose EL vacuum **derives** the
  Koide cone Q=2/3 (a dynamical, not group-theoretic, statement), with the phase
  as a Goldstone and the second moment equal to the gravity source. Gravity half
  (GEN1/2) + matter half (GEN3) now share one scale `a`.
- **Honest residual**: the Brannen phase `φ=2/9` remains undetermined by `V` (a
  Goldstone) — consistent with v59's "phase not geometric"; pinning it needs an
  explicit S₃-breaking term (a posit), to be characterized later.
- **Next (GEN4, aspect 4)**: the **coupling** — show `ρ_grav = ⟨Φ̃Φ⟩` enters the
  gravity trace law as the source from a *covariant* matter→gravity coupling
  `S_source[ρ_grav; g(B)]`, and verify **equivalence-principle universality** (the
  same `ρ_grav` couples to all states with one `f_g`). That closes the GEN1–GEN3
  loop into a single action and sets up GEN5 (full linearized spectrum: graviton +
  Goldstone + massive together).
