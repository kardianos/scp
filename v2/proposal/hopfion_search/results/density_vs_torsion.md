# Density Conservation vs Field Torsion: The Missing Mechanism

## The Disconnect

The CHPT narrative spec (Ch 02, Ch 07) builds gravity on **density conservation**:

> A knot concentrates density locally, depleting the surrounding field.
> Because total density is conserved (Axiom 3), the excess inside the knot
> must be exactly compensated by a deficit in the surrounding field.
> — spec/02_energy_and_density.md

This gives 1/r² gravity via geometric dilution of a conserved deficit:

    δρ(r) ~ M / (4πr²)

The mathematical spec (Lagrangian, field equations, numerical solutions) does NOT
produce this. Two regimes were tested:

### σ-model (λ → ∞): No depletion at all

|q(x)| = ρ₀ everywhere. The field norm is constant by constraint. The soliton
is a *twist* in the field direction, not a *concentration* of field density.
Total ∫ρ d³x is the same with or without the soliton.

### Finite-λ: Localized, exponentially recovering

ρ(r) varies near the core but recovers exponentially:

    ρ(r) ~ ρ₀ - Δρ · e^{-√(2λ)(r-R_core)}

At λ = 8000: ρ(0) = 0.803ρ₀, recovery length ~ 0.01 code = 0.006 fm.
No 1/r tail. No long-range depletion.

### The fundamental problem

The Lagrangian has no global density conservation law. Energy is extensive
(E = ∫H d³x), not redistributive. A soliton ADDS energy to the vacuum — it
does not borrow density from the surroundings. The vacuum far from any soliton
is always ρ = ρ₀ regardless of how many solitons exist.

This was confirmed across all six gravity paths investigated (see gravity_paths.md).
None produced 1/r gravity from the existing Lagrangian.

---

## The Key Difference: What Is Conserved?

The narrative spec assumes:

    ∫ ρ(x) d³x = const  (total density conserved over all space)

The math spec has:

    ∫ B⁰(x) d³x = B ∈ Z  (topological charge conserved, exactly integer)
    ∫ T₀₀(x) d³x = E     (energy conserved, but not a fixed constant)

**Density is NOT conserved as a global constraint.** There is no equation forcing
∫ρ d³x = V·ρ₀. The sigma model enforces |q(x)| = ρ₀ pointwise (a local
constraint), not globally. Finite-λ has no constraint at all — ρ(x) is
dynamically determined, and the vacuum has more total ρ than any soliton state.

The narrative mechanism requires a conserved scalar quantity that:
1. Is concentrated by the soliton
2. Must be depleted elsewhere (conservation)
3. Decays as 1/r² (geometric dilution in 3D)
4. Produces an attractive force (pressure gradient toward depletion)

No scalar quantity in the Cl⁺(3,0,1) theory satisfies all four requirements.

---

## The Torsion Hypothesis

### Core idea

What if the non-local effect is not scalar density but **directional structure**?

The field Ψ ∈ Cl⁺(3,0,1) has 8 components: a scalar amplitude (s) and 7
directional components (f₁,f₂,f₃ bivector + j₁,j₂,j₃,p degenerate). The
σ-model constraint |q| = ρ₀ fixes the scalar amplitude but leaves the
directional structure free. A soliton is a *twist* — it rotates the field
direction, not the field magnitude.

**Hypothesis**: The soliton creates a non-local *torsion* in the field —
a 3D directional twist pattern that extends to infinity, even though the
scalar density is uniform everywhere.

### The wire analogy

An electron moving through a wire creates a magnetic field:
- The charge density of the wire is unchanged (neutral wire)
- The current creates a non-local *torsion* (magnetic field lines circling the wire)
- This torsion decays as 1/r (3D geometric dilution of the source current)
- Test charges respond to the torsion, not to any density change

In CHPT terms:
- The soliton (knot) is the current source
- The field density ρ = ρ₀ everywhere (neutral medium)
- The non-local effect is a *twist pattern* in the field direction
- This twist decays as 1/r or 1/r² depending on the multipole structure
- Other solitons respond to the twist, not to density

### What "torsion" means mathematically

The Cl⁺(3,0,1) field has a natural torsion interpretation. At each point,
the field defines an orientation frame via the rotor q = s + f₁e₂₃ + f₂e₃₁ + f₃e₁₂.
Far from any soliton, q → ρ₀ (identity rotation, no twist). Near a soliton,
q twists — the local frame rotates relative to the asymptotic frame.

The *torsion* of this frame field is the connection ∇q · q⁻¹, which measures
how fast the local orientation changes from point to point. For the B=1
hedgehog Skyrmion:

    Frame rotation angle: f(r) (radial profile, f: π → 0)
    Torsion magnitude: ~ f'(r) (derivative of profile)
    Asymptotic: f'(r) ~ -2e^{-r}/r → torsion decays exponentially

**Problem**: The hedgehog torsion decays exponentially (Yukawa), not as 1/r.
This is because the Skyrme field is *massive* — the vacuum has a mass gap
set by the soliton's natural frequency.

### The massless channel

For 1/r decay, need a MASSLESS torsion mode. In Cl⁺(3,0,1):

- Bulk sector (s,f): massive — gap from soliton curvature
- Degenerate sector (j,p): at e₀² = 0, these are NON-DYNAMICAL (algebraic
  constraints, no propagation). But if given small kinetic terms (e₀² > 0),
  the pseudoscalar p would be massive (μ > 0) while the pseudovector j has
  mass set by the coupling structure.

The constraint propagator (e₀² = 0 limit) is 1/μ² — instantaneous, not 1/r.
This is the Coulomb-like propagator that gives 1/r potential (Path 3: B⁰p coupling).
The torsion interpretation adds geometric content: the pseudoscalar p IS the
frame torsion component (scalar part of the Clifford torsion tensor).

### Connection to existing results

| Mechanism | Scalar/Tensor | Range | In current math? |
|-----------|--------------|-------|-------------------|
| Density depletion | Scalar | Exponential | No (σ-model: ρ=ρ₀) |
| B⁰p coupling (Path 3) | Scalar | 1/r (massless p) | Yes, but g_top free |
| BLV effective metric (Path 4) | Tensor | Core-scale | Yes, from L₂+L₄ |
| **Frame torsion** | **Vector/Tensor** | **1/r if massless** | **Not yet tested** |

The torsion hypothesis adds a new geometric interpretation to Path 3: the
pseudoscalar p is not just an algebraic field — it is the component of frame
torsion that couples to the topological charge B⁰. The coupling g_top measures
how strongly the soliton's twist pattern sources the torsion field.

---

## What Would Need To Be True

For the torsion mechanism to produce gravity:

1. **The torsion must be physical**: The frame field q/|q| must have independent
   degrees of freedom beyond the soliton profile. In the σ-model, the frame IS
   the field (fully determined by the configuration). No extra torsion.

2. **The torsion must be massless**: For 1/r decay. The current theory has no
   massless modes in the bulk sector. The degenerate sector at e₀² = 0 has
   instantaneous propagation (effectively massless Coulomb).

3. **The coupling must be universal**: Every soliton must source the torsion
   proportional to its mass. Path 3 has B⁰ ∝ mass for single species, but
   the coupling g_top is a free parameter.

4. **The torsion must produce tensor gravity**: 1/r scalar potential gives half
   the light bending of GR. Need spin-2 for full GR.

### The density conservation question

If we IMPOSE density conservation as an axiom (not derived from the Lagrangian),
then:

    ∫ ρ(x) d³x = V·ρ₀  (V = total volume)

A soliton with local ρ > ρ₀ FORCES ρ < ρ₀ somewhere else. In the σ-model
(|q| = ρ₀ pointwise), this is automatically satisfied — no depletion, no
gravity. In a theory with variable ρ, this constraint would be a new GLOBAL
conservation law not present in the current Lagrangian. It would require either:

- A constraint (Lagrange multiplier) enforcing constant total density
- A modified Lagrangian where the field equation implies ∫ρ d³x = const
- Or: reinterpret "density" as something other than |q|

The third option connects to torsion: if "density" means the FULL field content
(including directional structure, not just amplitude), then a soliton that
concentrates *twist* locally must deplete *twist capacity* elsewhere. This is
closer to the narrative mechanism but operates on frame orientation, not scalar
density.

---

## Numerical Results (torsion.c)

### Test A: Decay Rates (σ-model)

All quantities computed on the B=1 hedgehog at σ-model (ρ=ρ₀=1):

| Quantity | Definition | Decay law (fit r=5-15) | Expected |
|----------|-----------|----------------------|----------|
| Frame deviation D(r) | 2\|sin(f/4)\| | **2.63/r^2.11** | ~2/r² |
| Frame torsion \|ω(r)\| | √(f'²+2sin²f/r²) | **10.4/r^2.99** | 1/r³ |
| Baryon density B⁰(r) | sin²f\|f'\|/(2π²r²) | **10.8/r^9.16** | 1/r⁹ |
| **Path 3 potential p(r)** | Q_enc/(4πr) | **0.0795/r^1.00** | **1/r** |

**Key**: Only the Path 3 constraint potential gives 1/r. The bulk frame
torsion decays as 1/r³ — at r=10, the ratio p/|ω| ≈ 0.8 but Path 3
overtakes at r > 33 (crossover: 0.08/r = 10.4/r³ → r ≈ 33).

### Test B: Density Conservation

The soliton has LESS density at its core, not more:

| λ | ρ(0) | Deficit ΔQ | ΔQ/Q_vac | Yukawa range |
|---|------|-----------|----------|--------------|
| 100 | 0.970 | 0.517 | 4.6×10⁻⁶ | 0.07 code = 0.04 fm |
| 10000 | 0.9997 | 0.00516 | 4.6×10⁻⁸ | 0.007 code = 0.004 fm |

Deficit scales as ~1/λ → vanishes in σ-model limit (ρ=ρ₀ everywhere).
The narrative spec says "knot concentrates density" — the math says the
opposite: the knot *depletes* density at its core. The field amplitude is
suppressed where the twist is strongest.

### Test C: Density-Constrained Soliton

Adding chemical potential μ to enforce ∫ρ·4πr²dr = V·ρ₀:

At λ=100:
- μ* = 9.15×10⁻⁴
- Far-field shift: ρ_∞ - 1 = **4.58×10⁻⁶** (unmeasurable)
- ∫ρ matches vacuum to 10⁻¹⁶ (perfect conservation)
- Approach to ρ_∞ is **Yukawa** (exponential, NOT 1/r²)

The density mode mass m = √(2λ(3ρ_∞²-1)) ≈ 2√λ makes all density
perturbations Yukawa. Conservation redistributes ~0.5 code³ of density
over volume ~10⁵ code³ → per-point shift of ~5×10⁻⁶. No long-range
force gradient.

### Test D: Boosted Soliton

| d (transverse) | D(d) frame dev | \|v_boost\|/v | p₃(d) Path 3 |
|---|----|----|----|
| 2 | 4.59×10⁻¹ | 2.23×10⁻¹ | 3.98×10⁻² |
| 10 | 2.08×10⁻² | 2.08×10⁻³ | 7.96×10⁻³ |
| 20 | 3.80×10⁻³ | 1.90×10⁻⁴ | 3.98×10⁻³ |

The boost velocity field decays as ~1/d³, NOT 1/d (unlike EM magnetic
field from wire current). The wire analogy holds ONLY through the Path 3
constraint sector: □p = g_top·B⁰ gives both 1/r scalar and 1/r vector
(gravitomagnetic) potential.

---

## Experimental Signature

The torsion hypothesis makes a specific prediction different from density depletion:

**Density depletion**: A propagating soliton leaves a density wake (ρ < ρ₀ behind
it, ρ > ρ₀ in front). Gravitational effects correlated with density variations.

**Frame torsion**: A propagating soliton drags its twist pattern with it. No
density wake, but a rotating/twisting field pattern that moves with the soliton.
Gravitational effects correlated with frame gradients, not density gradients.

In the 3D simulations (scatter.c), this difference is testable: measure
both ρ(x) and the frame orientation θ(x) during soliton propagation. Does the
long-range perturbation look like a density wave or a twist wave?

**Preliminary evidence**: In the B+B̄ scattering simulation (v=0.5c), the
topological charge Q = 0 throughout — the fields PASS THROUGH each other.
This means the long-range field perturbation is in the frame, not the density.
A density-based mechanism would predict different scattering dynamics than a
frame-based one.

---

## Summary

| | Narrative Spec | Math Spec | Torsion Hypothesis |
|---|---|---|---|
| Conserved quantity | Scalar density ρ | Topological charge B | Frame orientation |
| Soliton effect | Concentrates ρ | Winds the frame | Twists the frame |
| Far-field signal | ρ < ρ₀ (depletion) | ρ = ρ₀ (no signal) | Asymptotic frame rotation |
| Range | 1/r² | Exponential | 1/r if massless mediator |
| In current Lagrangian? | No | Yes | Partially (Path 3) |
| New physics needed? | Global ρ conservation | None | Massless torsion mode |

The key open question: does the Cl⁺(3,0,1) algebra have a massless torsion
mode that couples universally to topological charge? This would unify the
narrative mechanism (non-local field distortion from solitons) with the
mathematical structure (constraint propagator, B⁰p coupling) and potentially
explain why g_top has the specific value it does.

---

## Verdict

### Hypothesis 1 (density conservation → 1/r² gravity): NEGATIVE

Density is NOT conserved in the Lagrangian. The soliton has LESS density
at its core (ρ₀ = 0.970 at λ=100), not more. Even with a forced conservation
constraint (chemical potential), the density mode is massive (m = 2√λ), so
perturbations are Yukawa. The constraint shifts the far field by ~10⁻⁶ ρ₀.
No 1/r² depletion zone exists.

### Hypothesis 2 (frame torsion → non-local 1/r effect): PARTIALLY CONFIRMED

The soliton IS a frame twist (not a density concentration) — this is correct.
The twist extends non-locally with power-law decay (1/r² for frame deviation,
1/r³ for torsion connection). But these decay too fast for gravity. Only the
Path 3 constraint channel (p sourced by B⁰) gives 1/r.

### Hypothesis 3 (wire analogy → magnetic torsion): NEGATIVE for bulk, POSITIVE for Path 3

The boosted soliton's frame velocity decays as 1/d³ (not 1/d). The wire
analogy DOES work but only through the degenerate sector constraint field:
□p = g_top·B⁰ gives both static 1/r (gravitoelectric) and moving 1/r
(gravitomagnetic) potentials. The pseudoscalar p IS the "torsion vector"
the theory needs — it is the scalar component of the Cl⁺(3,0,1) frame
torsion. But g_top remains free.

### The path forward

The density conservation mechanism CANNOT work with a massive density mode.
Two options:

1. **Flat-bottom potential**: Replace V = (λ/4)(ρ²-ρ₀²)² with a potential
   having V''(ρ₀) = 0. This makes the density mode massless → 1/r depletion
   → density-based gravity. But: flat potential may not support solitons
   (Derrick's theorem requires a mass gap for stability).

2. **Accept Path 3**: The constraint propagator (e₀²=0 → p sourced by B⁰)
   already gives 1/r gravity with the correct sign, gravitomagnetic effects,
   and universal coupling. The only problem is that g_top is a free parameter.
   The torsion interpretation adds geometric meaning: gravity = frame torsion
   of the degenerate Clifford sector, sourced by topological charge.

Code: `src/torsion.c` (self-contained, ~500 lines)
Run: `./bin/torsion [-lambda L]`
