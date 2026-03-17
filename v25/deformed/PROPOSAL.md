# V25 Phase 6: Deformed Oscillon — Tidal Spin-2 Radiation

## Thesis

An alternative to the binary: apply a TIDAL deformation to a single
oscillon and measure whether it radiates spin-2. A tidal force creates
a quadrupolar (l=2) deformation of the oscillon — exactly the shape
needed to emit spin-2 GW.

In GR, tidal forces are described by the Weyl tensor C_{ijkl}. A mass
in a tidal field deforms into an ellipsoid and radiates quadrupolar GW
as it oscillates.

In our model: apply an external quadrupolar potential and measure the
TT strain response.

## Setup

### The Tidal Potential

Add an external quadrupolar potential:

    V_tidal = -½ε_T · (x² - y²) · Σ_a φ_a²

This squeezes the oscillon along y and stretches it along x (or vice
versa depending on the sign of ε_T). The l=2, m=0 harmonic.

Alternatively, use a time-dependent tidal field:

    V_tidal = -½ε_T · cos(Ω_T t) · (x² - y²) · Σ_a φ_a²

This drives quadrupolar oscillations at frequency Ω_T. The radiated
TT strain should be at Ω_T with l=2 angular pattern.

### Why This Works

The spherical oscillon has Q_{ij} = 0 (no quadrupole moment). The tidal
field creates Q_{ij} ∝ ε_T · diag(+1, -1, 0). This oscillating Q
radiates TT strain at Ω_T.

The AMOUNT of radiation depends on the oscillon's TIDAL DEFORMABILITY
(Love number k₂ in GR language). A rigid oscillon (high k₂) deforms
little → weak radiation. A soft oscillon (low k₂) deforms more → strong.

## Method

### Phase 6a: Static tidal deformation

1. Equilibrate 3D oscillon (t=500, from V25 Phase 1)
2. Turn on static tidal field: V_tidal = -½ε_T(x²-y²)Σφ_a²
3. Scan ε_T = {0.001, 0.01, 0.1}
4. Measure the oscillon shape: compute the quadrupole moment
   Q_{ij} = ∫ (3x_ix_j - r²δ_{ij}) ρ(x) d³x
5. The tidal Love number: k₂ = Q_{22} / (ε_T · R⁵) where R is the
   oscillon radius

### Phase 6b: Oscillating tidal field → GW

6. Apply oscillating tidal field at frequency Ω_T = 0.1 (well below
   the breathing frequency ω ≈ 1.4 to avoid resonance)
7. Evolve t=500
8. Measure TT strain on sphere at R=15
9. Decompose into multipoles: expect l=2 dominance
10. Measure: h₊ and h× at equator (θ=π/2) and poles (θ=0)
    - Equator: maximum (tidal stretching is in xy-plane)
    - Poles: zero (no quadrupole moment along z)

### Phase 6c: Resonant driving

11. Scan Ω_T to find RESONANT response: the oscillon's quadrupole
    normal mode frequency
12. At resonance: the deformation (and GW emission) should be maximum
13. The resonant frequency is related to the shear mode frequency
    from V24-Combo 2+5 (ω ≈ 0.51 for the antisymmetric optical branch)

### Phase 6d: Self-consistent tidal field (from another oscillon)

14. Instead of external V_tidal: place a SECOND oscillon at distance D
    along the x-axis (NOT orbiting, just sitting there)
15. The second oscillon's gravitational field creates a tidal field at
    the first oscillon's location
16. Measure: does the first oscillon develop a quadrupole moment
    pointing TOWARD the second?
17. This is the TIDAL INTERACTION — the precursor to binary dynamics

## Connection to Spin-2

The tidal deformation directly measures the oscillon's coupling to
l=2 perturbations. If the response is:
- Strong and at l=2: the oscillon couples to spin-2 → emits spin-2 GW
- Weak or at l=0: the oscillon is rigid to tidal forces → weak GW

The Love number k₂ quantifies this coupling. In GR, neutron stars have
k₂ ≈ 0.05-0.15. Black holes have k₂ = 0 (maximally rigid). Our
oscillon's k₂ determines its "gravitational softness."

## Implementation Notes

The tidal potential V_tidal = -½ε_T(x²-y²)Σφ² adds a force:
    acc_a += ε_T · (x² - y²) · φ_a    (for static)
    acc_a += ε_T · cos(Ω_T t) · (x² - y²) · φ_a    (for oscillating)

For the quadrupole moment measurement:
    Q_{ij} = Σ_gridpoints (3x_ix_j - r²δ_{ij}) · ρ · dV

where ρ = energy density and dV = dx³. The key components:
    Q_{20} ∝ ∫ (2z² - x² - y²) ρ dV    (prolate/oblate)
    Q_{22} ∝ ∫ (x² - y²) ρ dV          (xy-plane quadrupole)

For the TT projection at distance R:
    h₊ = (σ_{θθ} - σ_{φφ}) / 2
    h× = σ_{θφ}

where σ is the traceless strain in spherical coordinates.

## Reference Code

- v25/src/v25.c (Phase 1-4, strain diagnostics, self-consistent metric)
- v24/fundamental/testB_inertia (deformation measurement)

## Output

- `src/deformed.c` — tidal deformation + GW emission code
- `data/tidal_static_eps{E}.tsv` — quadrupole moment vs ε_T
- `data/tidal_driven_omega{O}.tsv` — GW emission vs driving frequency
- `data/tidal_gw_multipoles.tsv` — angular decomposition
- `data/tidal_love.tsv` — Love number measurement
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ_pw=0.5, η=0.1, λ_L=0.1, α_g=0.001
ε_T scan: {0.001, 0.01, 0.1}
Ω_T scan: {0.05, 0.1, 0.2, 0.5, 1.0}
N=96, L=15 (single oscillon), t=500

Compile: `gcc -O3 -fopenmp -Wall -o deformed src/deformed.c -lm`
