# V24-180B: UUD/UDD in the 180° Anti-Phase State

## Thesis

The 180° state (φ, φ, -φ) naturally provides a two-flavor structure:
two "positive" fields (up) and one "negative" field (down). Test the
UUD and UDD mass ordering in the 180° regime where binding IS present
(P = -f³), unlike the 120° regime where P≈0.

The V24-UDD result at 120° showed correct neutron > proton ordering.
Does this ordering persist in the 180° state where the oscillon is
actually bound?

## Setup

### UUD (proton): two up fields + one anti-phase down field
    φ₁ = +A_U·g(x), m₁ = m_U   [Up]
    φ₂ = +A_U·g(x), m₂ = m_U   [Up]
    φ₃ = -A_D·g(x), m₃ = m_D   [Down, anti-phase]

### UDD (neutron): one up field + two anti-phase down fields
    φ₁ = +A_U·g(x), m₁ = m_U   [Up]
    φ₂ = -A_D·g(x), m₂ = m_D   [Down, anti-phase]
    φ₃ = -A_D·g(x), m₃ = m_D   [Down, anti-phase]

Note: UDD has TWO anti-phase fields. The triple product:
P_UDD = A_U·(-A_D)·(-A_D)·g³ = +A_U·A_D²·g³ (POSITIVE P)
P_UUD = A_U·A_U·(-A_D)·g³ = -A_U²·A_D·g³ (NEGATIVE P)

The signs differ! UUD has P < 0, UDD has P > 0. With the saturating
potential V = (μ/2)P²/(1+κP²), the sign of P doesn't matter for the
potential energy (it depends on P²). But the FORCE ∂V/∂φ_a depends on
P linearly, so the dynamics differ.

## Method

1. Fix m_U = 1.0, A_U = A_D = 0.8, σ = 3.0
2. Scan m_D ∈ {1.0, 0.95, 0.90, 0.85, 0.80, 0.70}
3. For each m_D: run BOTH UUD and UDD configurations for t=10000
4. Measure: E_final, ω per field, fc, stability
5. Compare E_UUD vs E_UDD at each m_D → mass ordering

## Key Questions

1. Does the UDD > UUD mass ordering persist in 180°?
2. Are the 180° UUD/UDD oscillons more stable than the 120° versions?
3. Does the sign of P (negative for UUD, positive for UDD) affect dynamics?
4. At what m_D/m_U ratio does the mass ordering flip?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/uud180.c` — solver with per-field masses and 180° initialization
- `data/uud180_mD{val}_ts.tsv`, `data/udd180_mD{val}_ts.tsv`
- `data/mass_ordering.tsv` — E_UUD vs E_UDD per m_D
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m_U=1.0
m_D scan: {1.0, 0.95, 0.90, 0.85, 0.80, 0.70}
Nx=4000, xmax=100, tfinal=10000

Compile: `gcc -O3 -Wall -o uud180 src/uud180.c -lm`
