# V25: Three-Field Elastic Gravity in 3D

## Vision

The field IS spacetime. Three scalar fields φ₁, φ₂, φ₃ are the displacement
of an elastic medium in the three spatial directions. The metric is not
imposed — it is COMPUTED from the fields at each point. Causality (ds = c·dt)
is not a separate equation — it is the propagation speed of the fields
themselves. Spin-2 gravity is not added — it is the transverse-traceless
strain of the displacement field, which automatically has spin-2 character.

Everything from V24 is retained: the oscillon (matter), the pairwise coupling
(forces), the mass gap (stability). What changes is the INTERPRETATION and
the DIMENSIONALITY: 3D fields with field index = spatial direction.

## The Principle: ds = c·dt FROM the Fields

In V24, we bolted on □Φ = -αρ as a separate PDE for gravity. This worked
(causal, correct static limit) but was artificial — Φ was an independent
field with its own dynamics.

In V25, there is NO separate Φ. The metric g_{ij} is derived from the
fields themselves:

    g_{ij}(x,t) = δ_{ij} + h_{ij}(x,t)

    h_{ij} = ∂_i φ_j + ∂_j φ_i    (linearized strain)

The fields φ_a propagate at speed c (from the Lagrangian). Therefore h_{ij}
propagates at speed c. Therefore changes in the metric propagate at c.
Therefore ds = c·dt. No separate equation needed.

## The Lagrangian

### Matter Sector (from v21)

    L_matter = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)(∂_i φ_a) - ½m²φ_a²]
               - V(P)

    V(P) = (μ/2)P²/(1+κP²),  P = φ₁φ₂φ₃

This gives the oscillon: a breathing bound state of three fields.

### Pairwise Coupling (from V24-ME)

    L_pw = -λ(φ₁φ₂ + φ₂φ₃ + φ₃φ₁)

This splits the mass spectrum: m²_S = m²+2λ (symmetric, oscillon),
m²_A = m²-λ (antisymmetric, mediator).

### Elastic Coupling (NEW in V25)

In standard elasticity, the Lagrangian depends on the strain INVARIANTS:

    I₁ = tr(ε) = ∂_i φ_i = div(φ)          (compression)
    I₂ = tr(ε²) = ε_{ij}ε_{ij}              (shear + compression)

where ε_{ij} = ½(∂_i φ_j + ∂_j φ_i) is the symmetric strain tensor.

The elastic Lagrangian:

    L_elastic = -½λ_L(div φ)² - μ_S·ε_{ij}ε_{ij}

where λ_L is Lamé's first parameter and μ_S is the shear modulus.

Note: the existing gradient term -½(∂_i φ_a)(∂_i φ_a) already contains
part of the elastic energy. The ADDITIONAL elastic terms are:

    L_cross = -½η(∂_i φ_j)(∂_j φ_i)     (cross-gradient, from V24-B)
    L_lame  = -½λ_L(∂_i φ_i)²           (compression stiffness)

These couple field index to spatial index — the elastic interpretation.

### Self-Consistent Metric (the core of V25)

The fields evolve on the metric they create. The effective EOM:

    g^{μν} ∇_μ ∇_ν φ_a - m²_eff φ_a - forces = 0

where g^{μν} depends on the strain (and thus on φ_a), creating a
NONLINEAR, SELF-REFERENTIAL system.

In the weak-field limit (|h_{ij}| << 1):

    ∂²φ_a/∂t² = (δ_{ij} + h_{ij})∂_i∂_j φ_a
                 + (∂_i h_{ij})∂_j φ_a
                 - m²(1 + ½h_{kk})φ_a
                 - forces

The metric correction h_{ij} = ∂_i φ_j + ∂_j φ_i enters the Laplacian,
the gradient coupling, AND the effective mass. All three effects are
computed from the SAME fields — no separate Φ.

## Why Spin-2 Emerges

### The Strain Decomposition

The strain tensor ε_{ij} decomposes into:

    ε_{ij} = ⅓θδ_{ij} + σ_{ij}

where θ = ε_{kk} = div(φ) is the trace (compression, spin-0) and
σ_{ij} = ε_{ij} - ⅓θδ_{ij} is the traceless part (shear, spin-2).

In 3D, σ_{ij} has 5 independent components — matching the 5 polarizations
of a massive spin-2 field. For a wave propagating in the z-direction,
the two transverse-traceless (TT) polarizations are:

    h₊ = σ_{xx} - σ_{yy}     (plus polarization)
    h× = 2σ_{xy}              (cross polarization)

These are EXACTLY the gravitational wave polarizations.

### The Key Question

Do the TT strain modes propagate as INDEPENDENT degrees of freedom?

In standard elasticity: YES. Shear waves propagate at speed c_T = √(μ_S/ρ),
independently from compression waves at c_L = √((λ_L+2μ_S)/ρ).

In our model: the mass term m²φ² gives ALL modes the same gap. But the
pairwise coupling (V24-ME discovery) SPLITS the spectrum. And the elastic
couplings (η, λ_L) modify the speeds differently for compression and shear.

The combination of pairwise coupling + elastic terms could produce:
- Compression (θ): massive, speed c_L, gap m_S
- Shear (σ_{ij}): lighter or gapless, speed c_T, gap m_A

If m_A → 0 (from large λ, as in V24-P1): the shear mode becomes
MASSLESS and SPIN-2 — the graviton.

## Implementation Plan

### Phase 1: 3D Elastic Oscillon

1. Extend v21/triad3d.c to include:
   a. Cross-gradient term: η(∂_i φ_j)(∂_j φ_i)
   b. Lamé compression term: λ_L(∂_i φ_i)²
   c. Pairwise coupling: λ_pw(φ₁φ₂+φ₂φ₃+φ₃φ₁)
2. Grid: N=128, L=20 (same as V24-S4 which confirmed 3D survival)
3. Parameters: μ=-20, κ=20, m=1.0, λ_pw=0.5, η=0.1, λ_L=0.1
4. Verify: does the 3D oscillon survive with elastic couplings?
5. Measure: breathing frequency, energy, profile

### Phase 2: Strain Field Diagnostics

6. At each diagnostic step: compute the full strain tensor ε_{ij}
7. Decompose into θ (trace) and σ_{ij} (traceless)
8. Compute the TT projection σ^TT_{ij} in the far field
9. Measure: is there spin-2 content in the oscillon's radiation?
10. Compare with V24-S6 (which found spin-0 WITHOUT elastic terms)

### Phase 3: Self-Consistent Metric

11. Compute h_{ij} = ∂_i φ_j + ∂_j φ_i at each grid point
12. Modify the Laplacian: ∂²φ_a/∂x_i² → (δ_{ij}+h_{ij})∂_i∂_j φ_a
13. Modify the mass: m² → m²(1 + ½tr(h))
14. This makes the evolution SELF-REFERENTIAL: fields create the metric,
    metric modifies the fields
15. Use weak-field (|h| << 1) and ramp up gradually (from V24-F lessons)

### Phase 4: Gravitational Wave Detection

16. Perturb the 3D oscillon (boost or deform)
17. Track the TT strain σ^TT_{ij} at distance r=20 from the oscillon
18. Does it propagate as a wave at speed c_T?
19. Does it have quadrupolar (l=2) angular pattern?
20. This is the DEFINITIVE test for emergent spin-2 gravity

### Phase 5: Two-Oscillon Gravitational Interaction

21. Two 3D oscillons at separation D=30
22. With self-consistent metric from Phase 3
23. Does the TT strain mediate attraction between them?
24. Is the force proportional to M₁M₂/r² (Newtonian)?
25. Is the force TENSOR (anisotropic) or SCALAR (isotropic)?

### Phase 6: Oscillon Lattice in 3D with Elastic Gravity

26. Small 3D lattice (2×2×2 = 8 oscillons)
27. With pairwise coupling + elastic terms + self-consistent metric
28. Phonon spectrum: compression branch + shear branch
29. Is the shear branch gapless? (This would be the graviton)
30. Speed of shear phonon = speed of gravity

## Computational Requirements

3D at N=128: 2.1M grid points × 3 fields × 3 arrays = 50 MB
Plus strain tensor (6 components): +100 MB
Plus metric correction: +50 MB
Total: ~200 MB per run, ~30 min at N=128 with OpenMP

For the lattice (Phase 6): N=256 needed for 8 oscillons at spacing 16
→ 16.8M points × ~200 MB = 1.3 GB. Expensive but feasible.

## What V25 Must Demonstrate

To claim emergent spin-2 gravity, ALL of these must hold:

1. □ The TT strain σ^TT propagates as an independent mode (not coupled
   to compression)
2. □ σ^TT is gapless (or nearly so) while compression is gapped
3. □ σ^TT has quadrupolar (l=2) radiation pattern
4. □ σ^TT mediates an ATTRACTIVE force between oscillons
5. □ The force is proportional to energy (universal coupling)
6. □ Propagation speed = c (from ds = c·dt, built into the fields)
7. □ All of this from a SINGLE Lagrangian with no separate Φ field

If even ONE of these fails, it's not spin-2 gravity. But the V24
results give reason for optimism: the pairwise coupling provides the
mass splitting (point 2), the elastic interpretation provides the
tensor structure (points 1, 3), and ds = c·dt is automatic (point 6).

## Parameters

μ=-20, κ=20, m=1.0
λ_pw = 0.5 (pairwise, from V24-E)
η = 0.1 (cross-gradient, to be scanned)
λ_L = 0.1 (Lamé compression, to be scanned)

Grid: N=128, L=20 (single oscillon)
      N=128, L=30 (two oscillons)
      N=256, L=80 (lattice)

Compile: `gcc -O3 -fopenmp -Wall -o v25 src/v25.c -lm`

## Reference Code

- v21/src/triad3d.c (3D solver, primary base)
- v24/fundamental/testE_lattice/src/lattice.c (pairwise lattice)
- v24/fundamental/testF_selfref/src/selfref.c (metric backreaction)
- v24/proca_3d/src/proca_3d.c (3D with pairwise coupling)
- v24/fundamental/combo_1257/src/combo1257.c (self-gravitating lattice)

## Success Criteria

**Minimum success**: 3D oscillon survives with elastic couplings and
self-consistent metric. TT strain is measurable in the far field.

**Moderate success**: TT strain propagates independently from compression.
Shear branch exists in the lattice phonon spectrum.

**Full success**: Shear branch is gapless (or gap controllable via λ_pw).
TT strain mediates attractive force between oscillons. Quadrupolar
radiation pattern confirmed. All from one Lagrangian, no separate Φ.

**Breakthrough**: Equivalence principle verified (all oscillons accelerate
identically in same strain field, regardless of mass). This would be
emergent General Relativity.
