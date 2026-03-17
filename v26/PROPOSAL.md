# V26: Braided Solitons — Topological Binding from Knot Theory

## Reference

Algebraic topology of knots, links, and braids:
https://aeb.win.tue.nl/at/algtop-5.html

Key concepts applied: braid group B₃, knot invariants (Jones polynomial,
writhe, genus), Seifert surfaces, Borromean linking, prime knot factorization.

## Vision

Replace the OSCILLON (breathing, scalar, monopolar, no topology) with a
BRAIDED SOLITON (static, topological, aspherical, protected by knot
invariants). The three displacement fields φ₁, φ₂, φ₃ form three strands
of a braid in physical 3D space. The braid topology prevents unwinding,
eliminates breathing, breaks spherical symmetry, and is NOT constrained
by the P/m=2 identity that killed the Skyrmion.

## Why Braids Instead of Skyrmions

| Property | Skyrmion | Oscillon | Braided Soliton |
|----------|---------|----------|-----------------|
| Topology | π₃(S³)=Z | None | Braid group B₃ |
| Ground state | Static | Breathing | Static |
| Radiation | None (topological) | Monopole (breathing) | None (topological) |
| Spherical symmetry | Yes (hedgehog) | Yes (breathing) | NO (braid breaks it) |
| P/m constraint | P/m=2 (kills gravity) | N/A | NO (B₃ ≠ SU(2)) |
| Spin | From isospin rotation | None | From writhe w |
| Mass | From topological degree B | From binding energy | From genus g |

## Mathematical Framework

### The Braid Group B₃

Three strands with generators s₁ (strand 1 crosses over strand 2) and
s₂ (strand 2 crosses over strand 3). Relations:

    s₁ s₂ s₁ = s₂ s₁ s₂      (Yang-Baxter relation)
    s_i s_j = s_j s_i  for |i-j| > 1   (far commutativity, trivial for n=3)

The braid word σ = s₁^{a₁} s₂^{a₂} ... encodes the crossing sequence.
The CLOSURE of the braid (connecting top endpoints to bottom) gives a
LINK in 3D — this is Alexander's theorem.

### Braid Closure → Link → Knot Invariants

Closure of braid σ₁σ₂ (one full twist) → trefoil knot 3₁
Closure of σ₁σ₂σ₁σ₂ → figure-eight knot 4₁
Closure of σ₁² → Hopf link (two-component)

Each link has:
- **Jones polynomial** V(t): computed from Kauffman bracket ⟨D⟩
- **Writhe** w: sum of crossing signs → angular momentum / spin
- **Genus** g: minimal genus of Seifert surface → mass/energy
- **Linking number**: pairwise linking of components → charge

### The Kauffman Bracket

For a crossing in the braid diagram, two resolutions:
- A-channel (left turn): coefficient A
- A⁻¹-channel (right turn): coefficient A⁻¹
- Loop removal: factor (-A² - A⁻²)

    ⟨D⟩ = A⟨D'⟩ + A⁻¹⟨D''⟩

The Jones polynomial: V(L) = (-A)^{-3w} ⟨D⟩, with t = A⁻⁴.

### Application to Three Displacement Fields

The three fields φ₁, φ₂, φ₃ in 3D define three "strands" — curves in
physical space along which each field has its maximum gradient. When these
strands BRAID (twist around each other), the crossing structure defines
a braid word, and the closure gives a link.

**Physical picture**: inside a tube of radius R around a core axis, the
three field directions rotate as you move along the axis. One full
rotation = one braid twist. The fields at any cross-section form a
pattern like three spokes of a wheel, rotated by an angle that depends
on position along the axis.

### The Borromean Configuration

Three rings, pairwise unlinked, triply linked. Removing any one frees
the other two. This mirrors the triple-product coupling:

    P = φ₁ φ₂ φ₃ ≠ 0 only when ALL THREE are present

A Borromean configuration of three field tubes would be:
- Individually: each φ_a is a smooth, non-knotted loop (unknot)
- Pairwise: any two φ_a, φ_b are unlinked
- Triply: all three together are LINKED (non-trivial linking number)

This is topologically protected by the triple linking, not by individual
field topology. It's the RELATIONSHIP between the three fields that creates
the topology — matching the triple-product coupling physics.

## Initialization: Braided Field Configuration

### Method 1: Twisted Tube

Place a tube of radius R along the z-axis. Inside the tube, the three
fields form a twisted pattern:

    φ₁(r,θ,z) = A(r) · cos(nz/L + 0)
    φ₂(r,θ,z) = A(r) · cos(nz/L + 2π/3)
    φ₃(r,θ,z) = A(r) · cos(nz/L + 4π/3)

where A(r) is a radial profile (Gaussian or sech), n is the number of
twists, and L is the tube length. With periodic BC in z: the closure
gives a closed braid.

For n=1: the three fields do one full 120° rotation along the tube.
This is a σ₁σ₂ braid → trefoil closure.

### Method 2: Toroidal Braid

Place the braid on a TORUS. The three field lines wind around the torus
with (p,q) winding numbers:

    φ_a on a path: (θ, ψ) = (t, q_a·t/p)  for t ∈ [0, 2πp]

where θ is the toroidal angle and ψ is the poloidal angle. Three
strands with different q_a values form a braid on the torus surface.

For (p,q) = (2,3): trefoil knot.
For (p,q) = (3,2): trefoil (same knot, different presentation).

### Method 3: Borromean Rings

Three interlocking rings in 3D, none passing through any other, but
triply linked. Each ring is a tube of field amplitude for one φ_a:

    Ring 1 (φ₁): circle in xy-plane, center at origin, radius R₀
    Ring 2 (φ₂): circle in xz-plane, center at origin, radius R₀
    Ring 3 (φ₃): circle in yz-plane, center at origin, radius R₀

Inside each ring tube (radius r_tube): the corresponding φ_a is large.
Outside: φ_a ≈ 0. The triple product P is nonzero only where all three
tubes overlap (near the origin, where all three rings intersect).

## The Lagrangian

Same as V25, but the INITIAL CONDITION is different (braided, not Gaussian):

    L = Σ_a [½(∂_t φ_a)² - ½(∂_i φ_a)(∂_i φ_a)]
        - ½m²Σ φ_a²
        - (μ/2)P²/(1+κP²)
        - λ_pw(φ₁φ₂+φ₂φ₃+φ₃φ₁)
        - ½η(∂_i φ_j)(∂_j φ_i)
        - ½λ_L(∂_i φ_i)²

The SAME Lagrangian as V25. What changes is the TOPOLOGY of the initial
field configuration. The dynamics then determines whether this topology
is preserved (stable braided soliton) or destroyed (unwinds to oscillon
or disperses).

## What Must Be Demonstrated (All in 3D)

### Phase 1: Braided Configuration Survival

1. Initialize a trefoil braid (Method 1 or 2) on N=128, L=20 grid
2. Initialize a Borromean rings configuration (Method 3)
3. Evolve each for t=1000
4. Does the topological structure persist?
5. Measure: crossing number, linking number at t=0 and t=1000
6. Compare energy with the oscillon at the same parameters

### Phase 2: Non-Breathing Verification

7. Track the total energy density ρ(t) at the soliton center
8. For the oscillon: ρ oscillates at ω (breathing)
9. For the braided soliton: ρ should be CONSTANT (no breathing)
10. Compute the DFT of ρ(0,t): the breathing peak should be ABSENT

### Phase 3: Aspherical Strain — Spin-2 Content

11. Compute the strain tensor ε_{ij} of the braided soliton
12. Decompose into l=0,1,2 multipoles on a sphere at R=10
13. The braid breaks spherical symmetry → l=2 should be LARGE
14. Compare with V25 Phase 2 (oscillon: 0.5% l=2)

### Phase 4: Writhe and Spin

15. Compute the writhe w of the field configuration numerically:
    w = (1/4π) ∫∫ (dr₁ × dr₂) · (r₁-r₂)/|r₁-r₂|³
    where the integral is over pairs of points on the braid
16. The writhe = angular momentum content = spin
17. Is w quantized (integer or half-integer)?

### Phase 5: Self-Consistent Metric on Braided Soliton

18. Compute h_{ij} = α_g(∂_i φ_j + ∂_j φ_i) from the braided fields
19. The strain field of a braid has NON-TRIVIAL tensor structure
    (unlike the oscillon's spherical strain)
20. Evolve with backreaction: fields on the metric they create
21. Does the braided soliton survive self-consistent gravity?

### Phase 6: Two-Soliton Interaction

22. Two braided solitons at separation D=30
23. They don't breathe → no monopole radiation → cleaner force measurement
24. Does the strain field mediate attraction?
25. Is the angular pattern of the force quadrupolar (l=2)?

## Computational Approach

All simulations MUST be 3D. The braided topology exists only in 3D —
it requires all three spatial dimensions for the strands to cross
over/under each other.

Grid: N=128, L=20 minimum (N=192 preferred for the Borromean rings
which need more spatial resolution).

The initialization is the hardest part — constructing a braided field
configuration on a discrete grid that has the correct topology. The
field values must be set carefully at each grid point to produce the
desired braid structure.

### Topological Verification

At each diagnostic step, verify the topology is preserved:
- Compute the linking number between field lines (iso-surfaces)
- Track the maximum of P (triple product) — should stay localized
- Count the number of "crossings" by tracking field-line intersections

### Numerical Knot Invariants

The Jones polynomial can be computed numerically from the field
configuration by:
1. Extract the braid word from the field line crossings
2. Apply the Kauffman bracket recursion
3. Multiply by (-A)^{-3w}

This gives a TOPOLOGICAL quantum number for the soliton.

## Why This Might Work for Gravity

The Skyrmion has P/m=2 because the SU(2) target space metric is
proportional to the kinetic metric. The braid group B₃ has NO METRIC —
it's a discrete algebraic structure. The "metric" experienced by waves
propagating through a braided medium depends on the CROSSING STRUCTURE,
not on a target-space metric.

The TT strain of a braided soliton is:
- Non-zero at l=2 (aspherical braid breaks symmetry)
- Static (no breathing → no monopole contamination)
- Topologically protected (can't unwind → persistent l=2 structure)

If the TT strain mediates attraction between braided solitons, and the
force has quadrupolar angular structure, that IS emergent spin-2 gravity
from a topological field configuration.

## Parameters

μ=-20, κ=20, m=1.0, λ_pw=0.5, η=0.1, λ_L=0.1, α_g=0.001
Braid: n=1 (trefoil), R_tube=3.0, A=0.8
Grid: N=128, L=20 (Phase 1-5), N=192, L=40 (Phase 6)

## Reference Code

- v25/src/v25.c (3D elastic solver with all couplings)
- v21/src/triad3d.c (base 3D solver)
- v2/hopfion_composition/ (3D field initialization with topological structure)

## Output

- `src/v26.c` — braided soliton 3D simulator
- `data/` — time series, strain multipoles, writhe, topology checks
- `RESULTS.md`

Compile: `gcc -O3 -fopenmp -Wall -o v26 src/v26.c -lm`

## Success Criteria

**Minimum**: braided configuration survives t=1000 without unwinding.
**Moderate**: non-breathing (no oscillation in ρ(0,t)), l=2 > 10%.
**Full**: self-consistent metric, strain-mediated attraction between two
braided solitons with quadrupolar angular pattern.
**Breakthrough**: Jones polynomial computed, writhe = spin, genus = mass.
