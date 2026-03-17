# V24-S5: What Determines λ? — Origin of the Pairwise Coupling

## Thesis

The pairwise coupling λ is currently a FREE PARAMETER. For the model to
be predictive, λ must either:
(a) Be determined by the other parameters (μ, κ, m) through some
    self-consistency condition
(b) Arise dynamically from a deeper mechanism
(c) Be promoted to a dynamical field (modulus) that evolves to an
    equilibrium value

This investigation explores three mechanisms for the origin of λ.

## Mechanism A: Self-Consistency from Oscillon Existence

The oscillon exists only when ω < m_S = √(m²+2λ). The oscillon frequency
ω depends on μ, κ, and the field amplitude, which in turn depends on λ
(because the pairwise coupling changes the effective potential).

At what λ does the system self-consistently satisfy ω(λ) < m_S(λ)?

The V24-S4 results show ω tracks m_S with a gap margin of 3-5%:
    ω ≈ 0.95-0.97 × m_S

This relationship might FIX λ: the oscillon adjusts its frequency to
maintain a fixed gap margin, which determines ω(λ), which determines
the amplitude, which determines λ self-consistently.

**Test**: Is the gap margin ω/m_S independent of λ? If it's CONSTANT
(say, always 0.97), then this is a universal relationship, not a
coincidence of the parameters.

## Mechanism B: λ from Vacuum Energy Minimization

The pairwise coupling changes the vacuum energy. With V_pw = λΣφ_aφ_b:
- At φ=0: V_pw = 0 (vacuum unaffected)
- At φ≠0: V_pw modifies the potential landscape

If the system seeks the LOWEST total energy (vacuum + oscillon), then λ
is determined by minimizing E_total(λ). Since the oscillon energy E_osc
increases with λ (more pairwise energy) but the binding energy also
increases (deeper well), there could be an optimal λ*.

**Test**: Compute E_total(λ) = E_oscillon(λ) + E_vacuum(λ) across the
full range of λ. Is there a minimum?

## Mechanism C: λ as a Dynamical Field (Modulus)

Promote λ to a scalar field Λ(x,t):
    V_pw = Λ(x,t) · (φ₁φ₂+φ₂φ₃+φ₃φ₁)

Add kinetic and potential terms for Λ:
    L_Λ = ½(∂Λ)² - U(Λ)

where U(Λ) has a minimum at some Λ₀. The field Λ evolves to Λ₀ everywhere
in vacuum, but near the oscillon it could take a different value due to
the coupling.

This is like a dilaton or modulus field in string theory — it sets the
coupling constants dynamically.

**Test**: With U(Λ) = ½m_Λ²(Λ-Λ₀)², evolve the combined system
(φ_a + Λ). Does Λ settle to a definite value near the oscillon? Is this
value different from Λ₀ (vacuum value)?

## What to Compute

### Test A: Gap Margin Universality

1. At each λ ∈ {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99}:
   equilibrate the oscillon, measure ω and m_S = √(m²+2λ)
2. Compute gap margin δ = 1 - ω/m_S
3. Plot δ vs λ. Is it constant?

### Test B: Energy Minimum

4. At each λ: measure the total oscillon energy E_osc(λ)
5. Plot E_osc vs λ. Is there a minimum?
6. Also compute E_osc/m_S (energy normalized by effective mass gap)

### Test C: Dynamical Modulus

7. Add a scalar field Λ(x,t) with mass m_Λ and equilibrium Λ₀
8. The pairwise force becomes: -Λ(x)·(φ_b+φ_c) for field a
9. Λ evolves: ∂²Λ/∂t² = ∂²Λ/∂x² - m_Λ²(Λ-Λ₀) + (source from φ)
10. Measure: Λ(x=0) vs Λ₀ at late times. Does the oscillon shift Λ?

## Reference Code

- v24/maxwell_e: pairwise coupling implementation
- v24/proca_push: λ scanning code
- v21/src/triad1d.c: base solver

## Output

- `src/proca_origin.c` — tests A, B, C
- `data/` — gap margins, energies, modulus profiles
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
Test A+B: λ scan {0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99}
Test C: Λ₀=0.5, m_Λ=0.1, source coupling g_Λ=0.1
Nx=4000, xmax=100, t=10000

Compile: `gcc -O3 -Wall -o proca_origin src/proca_origin.c -lm`
