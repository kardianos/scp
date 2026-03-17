# Test D: Local Breathing Frequency → Emergent Gauge Field

## Thesis

The oscillon breathes at frequency ω globally. Promote ω to a local
field ω(x,t). The requirement that physics is invariant under local
reparametrization of the breathing phase requires a gauge connection.

## Implementation

Instead of promoting ω directly (which isn't a field in the Lagrangian),
use a simpler approach: add a scalar field Ω(x,t) that modulates the
effective mass:

    m²_eff(x) = m² · (1 + g_Ω · Ω(x))

The field Ω evolves freely (massless or with small mass):
    ∂²Ω/∂t² = ∂²Ω/∂x² - m_Ω²·Ω

The oscillon's breathing frequency depends on m_eff:
    ω(x) ≈ 0.87 · m_eff(x)

If Ω varies in space → ω varies → different clock rates → gravitational
time dilation. This is NOT a gauge field but captures the PHYSICS of
local clock rates.

For a TRUE gauge field: parametrize each field as φ_a = A_a(x)·cos(θ_a(x,t))
where θ is the breathing phase. The gauge transformation: θ → θ + α(x).
This requires A_μ to maintain covariance. But this is complex — start
with the simpler Ω modulation.

## What to Test

1. Equilibrate oscillon (standard parameters)
2. Add the Ω field with g_Ω = 0.1, m_Ω = 0 (massless modulator)
3. Source: ∂²Ω/∂t² = ∂²Ω/∂x² + g_source · ρ(x,t) (sourced by energy density)
4. Does Ω develop a 1/|x| profile (in 1D: linear growth)?
5. Does the oscillon's ω shift in response to Ω?
6. Two oscillons: does Ω mediate a force?

## Reference Code

- v21/src/triad1d.c
- v24/maxwell_f (similar: massless scalar sourced by oscillon)

## Output

- `src/gauge.c`, `data/`, `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0
g_Ω = 0.1, m_Ω = 0, g_source = 0.01
Nx=4000, xmax=100, t=10000

Compile: `gcc -O3 -Wall -o gauge src/gauge.c -lm`
