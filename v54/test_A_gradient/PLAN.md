# Test A: Gradient-Based Theta Self-Interaction

## Concept

Where neighboring theta vectors are aligned (reinforcing), theta is stable.
Where they oppose (destructive interference), the gradient energy |∇θ|² is high
and the energy converts to phi density.

Chirality bias: modulate conversion rate by theta helicity h = θ·curl(θ).
One handedness is protected (low conversion), the other self-destructs into density.

## New terms

```
theta_acc[a] -= sigma * grad_theta_sq * theta[a]
phi_acc[a]   += sigma * grad_theta_sq * curl_theta[a]

where:
  grad_theta_sq = Σ_a |∇θ_a|²
  sigma_eff = sigma * (1 + chi * theta·curl(theta) / (|theta|² + eps))
```

## Parameters

- sigma: gradient conversion strength (test: 0.1, 0.5, 1.0)
- chi: chirality bias (test: 0, 0.5, 1.0) — chi=0 is symmetric, chi>0 breaks chirality

## What to look for

- Density forming at theta domain walls
- One chirality surviving while the other converts
- Self-sustaining phi structures sourced by theta gradients

## Config

Base: m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
BC: absorb T=30 then periodic, T=600
