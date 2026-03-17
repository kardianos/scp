# V24-S1b: Asymmetric Oscillon Proca Force Measurement

## Thesis

V24-S1 (clean force measurement) was inconclusive because the symmetric
oscillon (φ₁=φ₂=φ₃) does NOT couple to the antisymmetric Proca channel
(confirmed by V24-S2). The Proca mediator exists in the spectrum but is
never sourced.

This redo uses ASYMMETRIC oscillons that naturally couple to the Proca
channel. Two options:

### Option A: 180° Oscillon (φ, φ, -φ)

The 180° state has maximal antisymmetric content. The decomposition:
    S = (φ₁+φ₂+φ₃)/√3 = (f+f-f)/√3 = f/√3
    A₂ = (φ₁+φ₂-2φ₃)/√6 = (f+f+2f)/√6 = 4f/√6

The antisymmetric mode A₂ = 4f/√6 ≈ 1.63f — LARGE. This oscillon
strongly sources the Proca channel.

### Option B: UUD Oscillon (m_U ≠ m_D)

With per-field masses: m₁=m₂=1.0, m₃=0.95. The three fields equilibrate
to different amplitudes, creating nonzero antisymmetric content.

### Option C: Seeded Asymmetry

Start with a symmetric oscillon but add a permanent small antisymmetric
perturbation (e.g., scale φ₃ by 0.9). Track whether the asymmetry persists
under the dynamics with pairwise coupling.

## Method

### Phase 1: Equilibrate asymmetric oscillons

For each option (A, B, C) at λ=0.99:
1. Initialize the asymmetric oscillon
2. Evolve t=10000 for equilibration
3. Save the profile at a breathing peak
4. Measure: S, A₁, A₂ content. How much antisymmetric mode exists?

### Phase 2: Force measurement

For each option at D = {20, 30, 40, 60, 80}:
1. Place two equilibrated asymmetric oscillons
2. Evolve t=2000
3. Measure initial acceleration (quadratic fit, t∈[0,500])
4. Fit F(D) to Yukawa: extract range
5. Compare range with 1/m_A = 1/√(m²-λ) = 10 (at λ=0.99)

### Phase 3: Compare with symmetric

6. Run the same D values with SYMMETRIC oscillons at λ=0.99
7. The force from symmetric oscillons is pure tail overlap (short range)
8. The force from asymmetric oscillons should include Proca (long range)
9. The DIFFERENCE is the Proca-mediated force

### Key Prediction

At D=60 (6× the Proca range):
- Symmetric: F should be negligible (tail overlap decays as e^{-mr})
- Asymmetric: F should be detectable (Proca decays as e^{-m_A·D} = e^{-6})

If F(D=60, asym) >> F(D=60, sym): the Proca mediator IS active.

## The 180° Twist

The 180° oscillon has P = -f³ (same binding as 0° but opposite sign).
Two 180° oscillons interact through:
1. Symmetric tail overlap (same as 0°, short range)
2. Antisymmetric Proca exchange (NEW, range 1/m_A)

The Proca force between two 180° oscillons depends on their RELATIVE
orientation. If both are (φ,φ,-φ) with the SAME anti-phase field (φ₃),
the antisymmetric content is correlated → attractive Proca force.
If opposite anti-phase fields (one has φ₃ anti-phase, other has φ₁
anti-phase): the antisymmetric content is orthogonal → weaker coupling.

This orientation dependence is like a dipole interaction — the Proca
force has angular structure.

## Reference Code

- v24/proca_force: `/home/d/code/scp/v24/proca_force/src/proca_force.c`
  (equilibration + force measurement framework)
- v24/180char: `/home/d/code/scp/v24/180char/src/char180.c`
  (180° initialization)

## Output

- `src/proca_asym.c` — asymmetric oscillon force measurement
- `data/` — force tables, profiles, time series
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, λ=0.99
Options: A (180°), B (UUD, m_D=0.95), C (seeded, φ₃ scaled by 0.9)
D scan: {20, 30, 40, 60, 80}
Nx=8000, xmax=300, t_equil=10000, t_run=2000

Compile: `gcc -O3 -Wall -o proca_asym src/proca_asym.c -lm`
