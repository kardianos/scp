# V32 Binding-Weighted Gradient Coupling — Results

## GRAVITATIONAL ATTRACTION IN A SINGLE FIELD

| α | D(0) | D(300) | ΔD | E(300)/E(0) | Verdict |
|---|------|--------|-----|-------------|---------|
| 0.0 | 20.0 | 21.6 | +1.6 | 0.93× | Neutral (control) |
| **0.5** | **20.0** | **16.8** | **-3.2** | **2.4×** | **ATTRACTION** |
| 1.0 | 20.0 | 30.0 | +10.0 | 477× | Repulsion + energy growth |
| 2.0 | 20.0 | 30.5 | +10.5 | 460× | Repulsion + blowup |

## What Works

At α=0.5, two braids approach from D=20 to D=16.8 over T=300.
The braids are alive throughout (max_rho oscillates 2-3, fc≈0.14-0.21).
The binding weight w(P) suppresses the self-interaction at the core (w≈0.22)
while allowing the inter-braid coupling at the surface (w≈0.5).

## What Doesn't Work Yet

- α=1.0 (physical strength) causes repulsion + energy blowup (E grows 477×)
- α=0.5 has 2.4× energy growth (not catastrophic but not conservative)
- The "sweet spot" α=0.5 may be fine-tuned

## The Equation

    ∂²φ_a/∂t² = ∇²φ_a + w(P)×α×(∇ρ/ρ)·∇φ_a - m²φ_a - ∂V/∂φ_a

    w(P) = 1/(1 + |P|/P_thresh),  P_thresh = 10% of peak |P|
    ρ = energy density,  P = φ₀φ₁φ₂

Single field. No split. No smoothing. No c(ρ). Periodic BC.
Symplectic Velocity Verlet.

## Parameters

    N=128, L=30, dx=0.47, dt=0.056
    m²=2.25, μ=-41.3, κ=50, A_bg=0.1
    Two braids at (±10, 0, 0), D=20 initial

## Files

- src/v32_weighted_grad.c
- data/wgrad_a0.0/ — control (5 field snapshots + timeseries)
- data/wgrad_a0.5/ — ATTRACTION (5 field snapshots + timeseries)
- data/wgrad_a1.0/ — repulsion (5 snapshots)
- data/wgrad_a2.0/ — blowup (partial)
- data/wgrad_summary.tsv — all configs compared
