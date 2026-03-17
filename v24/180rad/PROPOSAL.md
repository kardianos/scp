# V24-180C: Radiation Comparison — 0° vs 180° Oscillon

## Thesis

The 0° oscillon has P = +f³cos³(ωt), which radiates at harmonics ω, 2ω, 3ω.
The 180° oscillon has P = -f³cos²(ωt)·cos(ωt) ... actually:

With φ₁=φ₂=+f·cos(ωt), φ₃=-f·cos(ωt):
P = f·cos(ωt) · f·cos(ωt) · (-f·cos(ωt)) = -f³·cos³(ωt)

So P₁₈₀ = -P₀. The TIME DEPENDENCE is identical (cos³ωt has the same
harmonic content regardless of sign). The radiation pattern should be
IDENTICAL in amplitude.

BUT: the SPATIAL structure could differ. The 180° state has φ₃ pointing
opposite to φ₁,φ₂. If the three fields couple to spatial gradients
differently (as in the elastic interpretation), the radiation DIRECTION
could be different.

Also: the force on each field depends on P·(∂P/∂φ_a). For the anti-phase
field (φ₃ = -f):
    ∂P/∂φ₃ = φ₁φ₂ = +f² (positive)
    P = -f³ (negative)
    P·(∂P/∂φ₃) = -f⁵ (negative)
    force = -μ·(-f⁵)/(1+κf⁶)² = +μf⁵/... = -|μ|f⁵/... (attractive toward -f)

For the in-phase fields (φ₁ = +f):
    ∂P/∂φ₁ = φ₂φ₃ = f·(-f) = -f² (negative!)
    P = -f³ (negative)
    P·(∂P/∂φ₁) = (-f³)·(-f²) = +f⁵ (positive)
    force = -μ·(+f⁵)/... = -μf⁵/... = +|μ|f⁵/... (attractive toward +f)

So BOTH the in-phase and anti-phase fields feel ATTRACTIVE forces toward
their respective values. The binding is symmetric.

**The radiation rate should be IDENTICAL to the 0° oscillon** because P²
(which determines the potential and the radiation source) is the same.

This test verifies this prediction numerically and checks for any
subtle differences from the sign asymmetry.

## Method

1. Evolve 0° and 180° oscillons for t=20000 at μ=-20, κ=20, m=1.0
2. Measure dE/dt for both over t ∈ [10000, 20000] (late-time steady state)
3. Compare harmonic amplitudes: DFT of P(t) for both
4. Compare spatial radiation patterns: DFT of φ₁(x_far, t) at large |x|
5. Check if the 180° state spontaneously transitions to 0° (or vice versa)

## Also Test: 180° with Asymmetric Amplitudes

If φ₁,φ₂ have amplitude A₊ and φ₃ has amplitude A₋:
- Equal: A₊ = A₋ = A (symmetric 180°)
- Asymmetric: A₊ = 1.2A, A₋ = 0.8A (broken symmetry)
- Strongly asymmetric: A₊ = 1.5A, A₋ = 0.5A

Does the asymmetry persist? Does it create a different radiation pattern?
This connects to the UUD/UDD investigation (V24-180B).

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/rad180.c` — radiation comparison code
- `data/rad_0deg_ts.tsv`, `data/rad_180deg_ts.tsv` — time series
- `data/rad_spectrum_0deg.tsv`, `data/rad_spectrum_180deg.tsv` — DFT
- `data/rad_asymmetric_*.tsv` — asymmetric amplitude tests
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
Nx=4000, xmax=100, tfinal=20000

Compile: `gcc -O3 -Wall -o rad180 src/rad180.c -lm`
