# V24-180A: Characterize the 180° Anti-Phase Oscillon

## Thesis

The 180° state (φ₁, φ₂, -φ₃) emerged spontaneously in V24-SYM and V24-PW.
It has P = -f³ (full triple-product binding with opposite sign). This
investigation fully characterizes the 180° oscillon: its frequency, energy,
profile, lifetime, and radiation rate, compared side-by-side with the
standard 0° oscillon.

## Key Questions

1. What is the 180° oscillon's breathing frequency ω₁₈₀?
2. What is its energy (mass) compared to the 0° state?
3. What is its radiation rate dE/dt? Better or worse than 0°?
4. Is the 180° state truly stable, or does it eventually collapse to 0°?
5. What does the harmonic spectrum of P(t) look like?
   (0° has P=+f³cos³ωt; 180° has P=-f³cos²ωt·cos(ωt+π) = +f³cos³ωt...
   Wait — if φ₃ = -f·cos(ωt), then P = f·f·(-f)·cos³ωt = -f³cos³ωt.
   The P spectrum should be IDENTICAL to 0° but with opposite sign.
   Does the sign of P matter for the dynamics?)
6. What is the profile shape? Same as 0° or different?

## Method

### Phase 1: Direct Initialization

Initialize the 180° state directly:
    φ₁(x,0) = +A·g(x)
    φ₂(x,0) = +A·g(x)
    φ₃(x,0) = -A·g(x)
    v₁ = v₂ = v₃ = 0

This gives φ₃ = -φ₁ at t=0. Since the EOM for φ₃ involves:
    ∂²φ₃/∂t² = ∇²φ₃ - m²φ₃ - μP·φ₁φ₂/(1+κP²)²

With P = φ₁φ₂φ₃ = -f³ (when φ₃=-f, φ₁=φ₂=+f):
    ∂P/∂φ₃ = φ₁φ₂ = +f²
So: force on φ₃ = -μ(-f³)(+f²)/(1+κf⁶)² = +μf⁵/(1+κf⁶)²

For the 0° state (all positive):
    force on φ₃ = -μ(+f³)(+f²)/(1+κf⁶)² = -μf⁵/(1+κf⁶)²

With μ < 0: the 0° force is +|μ|f⁵/... (attractive, pushes ω below gap).
The 180° force is -|μ|f⁵/... (REPULSIVE!).

**CRITICAL**: The 180° state has REPULSIVE triple-product force on the
anti-phase field! This means the binding is DIFFERENT from 0°.

Actually wait — let me reconsider. With φ₃ = -f (negative):
    m²φ₃ = -m²f (points TOWARD zero, restoring)
    The triple product force: -μP(∂P/∂φ₃)/(1+κP²)²

P = φ₁φ₂φ₃ = f·f·(-f) = -f³
∂P/∂φ₃ = φ₁φ₂ = f²

force_triple = -μ·(-f³)·f²/(1+κf⁶)² = +μf⁵/(1+κf⁶)²

With μ = -20: force_triple = -20f⁵/(1+κf⁶)² < 0

This points in the NEGATIVE direction (same direction as φ₃ = -f).
It pushes φ₃ MORE NEGATIVE — it's ATTRACTIVE toward the anti-phase state!

So the 180° state IS bound by the triple product, with the same force
magnitude as the 0° state. The sign of P doesn't matter for the binding
because the force depends on P·(∂P/∂φ_a), and the signs cancel.

### Phase 2: Side-by-Side Comparison

Run BOTH states to t=20000 with identical parameters:

0° state: φ₁=φ₂=φ₃ = +A·g(x)
180° state: φ₁=φ₂= +A·g(x), φ₃ = -A·g(x)

Measure for both:
- Breathing frequency ω (DFT of φ₁(0,t))
- Energy E(t) and loss rate dE/dt
- Peak amplitude A(t)
- Core fraction fc
- DFT of P(t) — harmonic content
- Profile f(x) at late times

### Phase 3: Stability Test

Perturb the 180° oscillon:
a. Flip perturbation: add +ε to φ₃ (push toward 0° state)
b. Amplitude perturbation: scale all fields by 1.2×
c. Asymmetric perturbation: different amplitudes for φ₁ vs φ₂

Does the 180° state resist these perturbations or collapse to 0°?

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/char180.c` — characterization code
- `data/char180_0deg_ts.tsv` — 0° control time series
- `data/char180_180deg_ts.tsv` — 180° time series
- `data/char180_stability_*.tsv` — perturbation tests
- `RESULTS.md`

## Parameters

μ=-20, κ=20, m=1.0, A=0.8, σ=3.0
Nx=4000, xmax=100, tfinal=20000

Compile: `gcc -O3 -Wall -o char180 src/char180.c -lm`
