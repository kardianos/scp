# V24-S2: Mediator Field Characterization

## Thesis

When two pairwise-coupled oscillons interact, what field configuration
exists between them? Is it purely the antisymmetric mode? Does it have
spatial structure (nodes, polarization)? How does it connect to the
oscillons at each end?

This is a diagnostic/visualization study, not a dynamics test.

## Method

### Step 1: Equilibrate + place two oscillons

1. Use the equilibrated profile from V24-S1 (or re-equilibrate at λ=0.99)
2. Place two oscillons at D=30 (within the Proca range ~10 at λ=0.99)
3. Evolve for t=2000

### Step 2: Decompose the field at each timestep

At each diagnostic output (every Δt=10):

Compute the symmetric and antisymmetric combinations:
    S(x) = (φ₁+φ₂+φ₃)/√3        (symmetric mode)
    A₁(x) = (φ₁-φ₂)/√2          (antisymmetric mode 1)
    A₂(x) = (φ₁+φ₂-2φ₃)/√6     (antisymmetric mode 2)

Output: S(x), A₁(x), A₂(x) as profiles at selected times.

### Step 3: Identify the mediating field

Between the two oscillons (D/4 < |x| < 3D/4):
- Is S nonzero? (symmetric mode leaked between oscillons)
- Is A₁ or A₂ nonzero? (antisymmetric mode = Proca mediator)
- What is the spatial profile? Exponential decay from each oscillon
  (Yukawa) or something else?

### Step 4: Frequency content of the mediator

At a point x_mid = 0 (midpoint between oscillons):
- DFT of S(x_mid, t), A₁(x_mid, t), A₂(x_mid, t)
- What frequencies are present? At the breathing frequency ω?
  At the Proca frequency √(k² + m_A²)? At DC (static component)?

### Step 5: Compare with and without pairwise coupling

Run the same two-oscillon setup at λ=0 (no pairwise) and λ=0.99 (strong
pairwise). Compare the field between the oscillons.

At λ=0: only the massive field tails overlap (range ~1)
At λ=0.99: the antisymmetric mode should extend much further (range ~10)

## Reference Code

- v24/proca_force (V24-S1) for equilibrated profiles
- v21/src/triad1d.c for base solver

## Output

- `src/proca_field.c` — field decomposition code
- `data/field_decomp_t{T}.tsv` — S, A₁, A₂ profiles at selected times
- `data/midpoint_spectrum.tsv` — DFT at midpoint
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ=0.99
D=30, Nx=8000, xmax=200, t_run=2000
Diagnostic output every Δt=10

Compile: `gcc -O3 -Wall -o proca_field src/proca_field.c -lm`
