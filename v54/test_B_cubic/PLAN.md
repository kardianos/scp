# Test B: Higher-Order Theta Self-Interaction (Cubic)

## Concept

Theta has a nonlinear self-interaction: -σ·|θ|²·θ (cubic damping).
This saturates theta growth without an arbitrary threshold.
When |θ| is small, the term is negligible. When large, it dominates.

Combined with conversion: the drained theta energy goes to phi via curl(θ).

## New terms

```
theta_acc[a] -= sigma * theta_sq * theta[a]
phi_acc[a]   += sigma * theta_sq * curl_theta[a]

where:
  theta_sq = theta[0]² + theta[1]² + theta[2]²
```

## Parameters

- sigma: cubic self-interaction strength (test: 0.1, 0.5, 2.0)

## What to look for

- Theta amplitude stabilizing (not growing without bound)
- Energy cycling between phi and theta
- Whether the cubic term alone creates localized structures
- Compare with saturation threshold (tests E/F) — this is smoother

## Difference from saturation (tests E/F)

The saturation used a hard threshold (theta_sat). The cubic term is smooth —
always active but only significant at large |θ|. This is more physical:
there's no sudden switch, just increasing resistance to theta growth.

## Config

Base: m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
BC: absorb T=30 then periodic, T=600
