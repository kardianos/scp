# V24-MF: Z₂ → U(1) Promotion via Complex Phase

## Background

V24-180 discovered a Z₂ symmetry: P → -P (equivalently φ₃ → -φ₃) leaves the
Lagrangian invariant. This discrete symmetry is the remnant of a continuous
U(1) if the triple product P is made complex.

Promote the real triple product to complex by allowing all three fields to
have complex phases. The modulus |P| stays fixed while the phase θ = arg(P)
rotates. If θ is gauged, the photon emerges from the Higgs mechanism.
If θ is NOT gauged, it's a massless Goldstone boson mediating 1/r.

## Setup

### Option 1: Ungauged (Goldstone)

Keep three REAL fields but add an angular variable θ(x,t) that represents
the "phase" of the triple product:

    P_eff = |P|·e^{iθ} where |P| = |φ₁φ₂φ₃|

The Lagrangian gets a kinetic term for θ:
    L_θ = ½f²[(∂_t θ)² - (∂_x θ)²]

where f is the "decay constant." The potential V depends only on |P|²
(already the case), so θ is an exact flat direction → massless Goldstone.

The Goldstone θ mediates a 1/r interaction (in 3D) or linear-in-|x|
(in 1D) between oscillons that carry "θ-charge."

### Option 2: Gauged (Higgs mechanism)

Promote θ to a U(1) gauge symmetry: θ → θ + α(x), A_μ → A_μ + ∂_μα.
The covariant derivative: D_μθ = ∂_μθ - eA_μ.

If the triple product has a VEV (⟨|P|⟩ = v > 0): the Goldstone θ is eaten
by A_μ, giving a massive photon with mass m_γ = ef.

If the triple product has NO VEV (⟨P⟩ = 0, which is the vacuum state):
the photon remains massless → Coulomb interaction.

For the oscillon: P ≠ 0 inside the core, P = 0 outside. The photon mass
varies with position: massive inside the oscillon (Meissner effect),
massless outside → long-range Coulomb.

### Implementation

The simplest version: add a REAL scalar θ field coupled to the oscillon
through a derivative coupling:

    L = L_oscillon + ½(∂θ)² - g·θ·(∂_t P)

where g·(∂_t P) acts as a current sourcing θ. This is NOT a gauge theory
but captures the essential physics: the time-varying triple product radiates
θ waves, which propagate at c (massless) and create 1/r interactions.

## What to Compute

### Phase 1: Goldstone Mode

1. Add a massless scalar θ field to the 1D triad code
2. Couple θ to the oscillon via: □θ = g·∂_t(P²)
   (derivative coupling ensures θ is sourced by oscillating P, not static P)
3. The oscillating oscillon (P ~ cos³ωt) sources θ at frequencies ω, 3ω
4. Measure: θ profile at large |x|. Does it fall as 1/|x| (1D: linear)?
5. Scan coupling g ∈ {0.01, 0.1, 1.0}

### Phase 2: Two-Oscillon Interaction

6. Two oscillons at separation D, both sourcing θ
7. The θ field mediates an interaction: attractive or repulsive?
8. Measure force vs D. Is it 1/D² (3D prediction) or constant (1D)?

### Phase 3: Backreaction

9. Does the θ field affect the oscillon? (Does it drain energy?)
10. At what coupling g does the θ radiation become significant?

## Reference Code

- v21 1D: `/home/d/code/scp/v21/src/triad1d.c`
- v22 gravity: `/home/d/code/scp/v22/src/two_oscillon.c` (massless mediator)

## Output

- `src/maxwell_f.c` — solver with Goldstone θ field
- `data/` — TSV output
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
g scan: {0.01, 0.1, 1.0}
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o maxwell_f src/maxwell_f.c -lm`
