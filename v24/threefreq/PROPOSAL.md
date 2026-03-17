# V24-D: Three-Frequency Non-Degenerate Oscillon

## Thesis

The standard oscillon has ω₁=ω₂=ω₃=ω (degenerate). The triple product
P = f³cos³(ωt) radiates at harmonics ω, 2ω, 3ω.

If the three fields oscillate at DIFFERENT frequencies, the radiation
pattern changes. Certain frequency combinations may radiate LESS than the
degenerate case.

This tests whether breaking the frequency degeneracy improves stability,
WITHOUT modifying the Lagrangian — only different initial conditions.

## Mathematical Setup

### Three-Frequency Triple Product

With φ_a = f_a cos(ω_a t + θ_a), the triple product is:

P = f₁f₂f₃ cos(ω₁t+θ₁) cos(ω₂t+θ₂) cos(ω₃t+θ₃)

Expanding with product-to-sum:

P = (f₁f₂f₃/4) [cos((ω₁+ω₂+ω₃)t + θ₁+θ₂+θ₃)
                + cos((ω₁+ω₂-ω₃)t + θ₁+θ₂-θ₃)
                + cos((ω₁-ω₂+ω₃)t + θ₁-θ₂+θ₃)
                + cos((-ω₁+ω₂+ω₃)t - θ₁+θ₂+θ₃)]

The four frequency components of P:
  Ω₁ = ω₁+ω₂+ω₃ (sum)
  Ω₂ = ω₁+ω₂-ω₃
  Ω₃ = ω₁-ω₂+ω₃
  Ω₄ = -ω₁+ω₂+ω₃

Radiation occurs when |Ω_k| > m (above the mass gap).

### Optimal Frequency Choice

To minimize radiation, we want ALL four |Ω_k| < m.

Constraint: each ω_a < m (otherwise the individual fields radiate).
Also ω_a > 0 (must oscillate).

The degenerate case ω₁=ω₂=ω₃=ω gives Ω₁=3ω, Ω₂=Ω₃=Ω₄=ω.
So Ω₁ > m when ω > m/3 ≈ 0.33.

For three DIFFERENT frequencies: can we push ALL Ω_k below m?

Ω₁ = ω₁+ω₂+ω₃ < m requires ω₁+ω₂+ω₃ < 1.0
Ω₂ = ω₁+ω₂-ω₃ < m requires ω₁+ω₂-ω₃ < 1.0
Similarly for Ω₃, Ω₄.

If all ω_a are equal at ω, then ω₁+ω₂+ω₃ = 3ω < 1 → ω < 0.33.
If frequencies differ: e.g., ω₁=0.4, ω₂=0.3, ω₃=0.2:
  Ω₁ = 0.9 < 1 ✓
  Ω₂ = 0.5 < 1 ✓
  Ω₃ = 0.3 < 1 ✓
  Ω₄ = 0.1 < 1 ✓

ALL below the gap! But can the triple-product binding support oscillons
with frequencies this low?

At v21 parameters (μ=-20, κ=20, m=1.0), the standard oscillon has ω=0.87.
Getting ω down to 0.4 requires MUCH stronger binding.

### Practical Test

Even at v21 parameters (ω₀≈0.87), breaking the degeneracy slightly may
help. Test with:
  ω₁ = ω₀ + δ, ω₂ = ω₀, ω₃ = ω₀ - δ  (symmetric splitting)

Then: Ω₁ = 3ω₀ (unchanged), Ω₂ = ω₀+2δ, Ω₃ = ω₀, Ω₄ = ω₀-2δ.

The splitting doesn't help Ω₁ (always 3ω₀ > m), but it SPLITS the
ω₀ component into ω₀±2δ. If this spreading reduces the resonant coupling
to radiation modes, the oscillon might live longer.

Also test ASYMMETRIC splitting: ω₁ = 0.9, ω₂ = 0.85, ω₃ = 0.80.
Then Ω₁ = 2.55, Ω₂ = 0.95, Ω₃ = 0.85, Ω₄ = 0.75.
The sum Ω₁ is still above gap, but Ω₂ is pushed right to the gap edge.

## What to Compute

### Phase 1: Direct Three-Frequency Initialization

1. For each frequency triplet (ω₁, ω₂, ω₃):
   Initialize three independent Gaussians with different velocities:
     φ_a(x,0) = A·g(x)
     v_a(x,0) = -ω_a·A·g(x)·sin(θ_a)  (with θ_a = 0 for simplicity)

   This gives each field its own oscillation frequency.

2. Evolve for t=10000. The triple product coupling will mix the frequencies,
   but the initial frequency content is set by the initialization.

### Phase 2: Scan

3. Test these frequency triplets (all with A=0.8, σ=3.0):

   a. Control: (0.87, 0.87, 0.87) — standard degenerate oscillon
   b. Small split: (0.90, 0.87, 0.84) — δ=0.03
   c. Medium split: (0.95, 0.87, 0.79) — δ=0.08
   d. Large split: (0.95, 0.80, 0.65) — asymmetric
   e. Low-frequency: (0.30, 0.30, 0.30) — all below m/3 (needs strong binding)
   f. Mixed low: (0.40, 0.30, 0.20) — all Ω_k < m

4. For cases (e) and (f): also try with μ=-60, -100 (stronger binding).

### Phase 3: Diagnostics

5. For each case: measure:
   a. Energy E(t) and dE/dt
   b. Per-field DFT to confirm the three frequencies
   c. DFT of P(t) to identify Ω_k components
   d. Oscillon lifetime and stability
   e. Whether the frequencies converge (to degenerate) or remain split

## Reference Code

- v21 1D solver: `/home/d/code/scp/v21/src/triad1d.c`

## Output

- `src/threefreq1d.c` — code
- `data/threefreq_{label}_ts.tsv` — time series per case
- `data/threefreq_summary.tsv` — comparison table
- `RESULTS.md` — analysis

## Parameters

μ=-20 (base), μ=-60,-100 (strong binding for low-frequency cases)
κ=20, m=1.0
Nx=4000, xmax=100, t_equil=0 (direct init), t_run=10000

Compile: `gcc -O3 -Wall -o threefreq1d src/threefreq1d.c -lm`
