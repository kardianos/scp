# Lengyel's PGA series — takeaways for v56

Read across:
- *Projective Geometric Algebra Done Right* (2020)
- *Symmetries in Projective Geometric Algebra* (2020)
- *Space-Antispace Transform Correspondence* (2022)
- *Poor Foundations in Geometric Algebra* (2024)
- *The Relativistic Inertia Tensor in PGA* (2024)
- *Relativistic Quaternions* (2024)

## Bottom line

Lengyel's series is a complete and honest **kinematic** framework. None
of the articles propose a dynamical equation, a wave equation, a field
theory, or a Lagrangian. The author's project is "the algebra, done
correctly" — fixing inconsistent definitions of inner products,
contractions, duals, and how rigid motions are represented. He is not
proposing a new physics theory.

This means **adopting his framework gives us better mathematical
machinery, not a new physics**. The dynamics are still ours to write.

## What we get from ℝ(3,1,1)

### Algebra structure

5 basis vectors with metric signature `(−1, +1, +1, +1, 0)`:

```
e₀  timelike     e₀² = −1
e₁  spacelike    e₁² = +1
e₂  spacelike    e₂² = +1
e₃  spacelike    e₃² = +1
e₄  projective   e₄² =  0
```

Total dimension = 2⁵ = 32 multivector components. The interesting
subspaces are smaller:

- **Vectors** (5 components): point in spacetime as
  `r = ct·e₀ + x·e₁ + y·e₂ + z·e₃ + e₄`
- **Bivectors** (10 components): lines in spacetime, momentum-energy,
  electromagnetic field tensors
- **Trivectors** (10 components): planes in spacetime
- **Even subalgebra** (16 components): rotations + boosts + screws
- **"Relativistic quaternion"** (8 components): even-grade elements
  that are linear combinations of `eᵢⱼ₀` and the spatial pseudoscalar
- **"Relativistic dual quaternion"** (12 components): general rigid
  spacetime motions

### Sandwich product is the universal rule

Any transformation of any object:
```
x' = Q ⟇ x ⟇ Q̃
```
where `⟇` is the geometric antiproduct and `Q̃` is the antireverse of
`Q`. The same formula transforms points, lines, planes, screws,
bivectors. No casts, no special cases. This is the elegance Lengyel is
selling.

### Lorentz invariance is not imposed — it's a property of the algebra

A pure spatial rotation, a Lorentz boost, a spacetime translation, and
combinations thereof are all sandwich products by elements of the even
subalgebra. The mass-shell relation `E²/c² − p² = m²c²` falls out of the
**bulk norm** of the energy-momentum bivector — it's not an axiom, it's
a definition.

### Spin and bulk-vs-weight split

Every multivector splits into:
- **Bulk** — the parts that don't involve `e₄` (the spatial geometry)
- **Weight** — the parts that do (the projective scaling)

The bulk norm has Minkowski signature `(+,−,−,−)` (or its dual,
depending on grade), so subluminal motion has real bulk norm,
luminal/lightlike has zero, superluminal has imaginary. This embeds the
causal structure into the algebra.

### Spinor structure

The relativistic quaternion (8 components, even subalgebra) is a
double cover of the proper orthochronous Lorentz group: under a 2π
rotation it picks up a sign, returning to itself only after 4π. **This
is the algebraic prerequisite for half-integer spin** — fermions, Pauli
exclusion, the Dirac equation.

The current 6-field Cosserat kernel does not have spinor structure.
The θ field is a pseudovector under SO(3); a 2π spatial rotation maps
it back to itself (no sign flip).

## What we don't get

Lengyel offers no equations of motion. Reading the series, we have to
build the dynamics ourselves. Three natural candidates:

### Candidate 1: Klein–Gordon on the multivector

The simplest Lorentz-invariant equation: each component of `M` satisfies

```
□ M = ∂²M/∂t² − c²∇²M = −m²M − ∂U/∂M
```

with a Lorentz-invariant potential `U(|M|², |M|⁴, …)`. This is closest
to what we currently do — replace the 6 Cosserat scalars with the 8
multivector components and write a similar potential.

Pros: minimal change to kernel; immediate mapping; no new operators.
Cons: pure scalar dynamics doesn't naturally produce fermions; the
algebra structure is barely used (only for invariance of the
potential).

### Candidate 2: Dirac equation on the relativistic quaternion

The natural equation when the field IS a spinor:

```
γ^μ ∂_μ ψ + m ψ = ∂L/∂ψ̄
```

where `ψ` is the 8-component relativistic quaternion field and `γ^μ`
are gamma matrices encoded as algebra elements. First-order in time
(needs a different integrator than current Verlet, or 2nd-order via
the squared Dirac equation `(□ + m²) ψ = …`).

Pros: this is the equation for fermions, gives spin-1/2 automatically.
Cons: requires major kernel redesign; first-order time; the gamma
matrix structure has to be encoded.

### Candidate 3: Nonlinear sigma model on the unit-norm subspace

Constrain `|ψ|² = ρ₀` and let ψ map spacetime → 7-sphere (the unit
relativistic quaternions). The map gets a topological winding number
which is conserved. This is the **Skyrme model** in its proper PGA
formulation:

```
Lagrangian = (½) ∂_μ ψ ∂^μ ψ + (Skyrme term: |[ψ⁻¹∂_μ ψ, ψ⁻¹∂_ν ψ]|²)
```

Pros: solitons are guaranteed by topology (π₃(S³) = Z), so
"particles" exist as a matter of mathematics.
Cons: the relevant homotopy group for relativistic quaternions
is not as well-studied as SU(2); we'd need to compute it; the Skyrme
term is fourth-order and computationally heavier.

## Recommended path for v56

A two-stage progression, both leveraging the existing foam kernel:

### Stage A: Klein–Gordon multivector field (Candidate 1)

- 8 scalar fields per cell instead of 6 (rename, don't rewrite)
- Same Laplacian operator, same Verlet integrator
- Lorentz-invariant potential built from algebra-derived scalars:
  bulk norm `|M|²`, geometric `M ⟑ M̃`, etc.
- Self-interaction tunable: linear term, quartic, sextic
- Look for solitons; expect Q-balls (non-topological solitons of
  Lorentz-invariant complex scalars) as the simplest case

### Stage B: Dirac multivector field (Candidate 2)

If Stage A produces interesting dynamics, upgrade to first-order
Dirac. This requires:

- New integrator (split-operator or Crank–Nicolson)
- Encoded gamma matrices (as 8×8 reals using the algebra)
- The non-relativistic limit reproduces Schrödinger
- Solitons here are *real* fermions if any exist

## Mapping to v55 foam kernel

Direct concordance — the cell layout is unchanged:

| v55 layout                | v56 layout                          |
|---------------------------|-------------------------------------|
| 6 scalar fields per cell  | 8 scalar fields per cell (mv components) |
| Laplacian on cells (FV)   | Laplacian on cells (same)           |
| Curl coupling η × curl(θ) | Algebra product M ⟑ M̃ self-coupling |
| V(P) = (μ/2) P²/(1+κP²)   | V(|M|²) = −μ |M|² + λ |M|⁴          |
| Verlet                    | Verlet (Stage A) or split-op (Stage B) |
| FMSH+FCEL+FCEP frames     | unchanged (just n_columns=8)        |

The SFA cell-native format already supports any `n_columns`. volview
needs minor changes (component-to-RGB mapping, derived scalar
visualisation).

## Concrete questions to answer before coding v56

1. **Which 8 components?** Lengyel's relativistic quaternion picks out
   `eᵢⱼ₀` and `e₃₂₁`. Are those the right components for our particle
   search, or should we use the full 16-component even subalgebra?
   (Probably the latter — gives us boosts AND rotations on equal
   footing.)

2. **How does the bulk norm enter the potential?**
   `V = −μ |M|²_bulk + λ |M|⁴_bulk` is the obvious choice; it's
   Lorentz-invariant by construction.

3. **What replaces the η × curl coupling?** In current SCP, this is
   the φ–θ exchange that gives EM-like behaviour. In algebra form,
   this might become a coupling between bulk and weight parts of the
   same multivector — a single field instead of two.

4. **Topological charge candidate?** Even if we go with Stage A
   (Klein–Gordon), can we identify a topological invariant of the
   multivector field? Lengyel's bulk-weight split + the antiproduct
   structure suggest a Hopf-like invariant might exist for
   appropriately constrained `M`.
