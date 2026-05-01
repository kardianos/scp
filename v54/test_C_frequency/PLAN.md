# Test C: Frequency-Matched Theta Interaction

## Concept

Theta at each voxel has a local angular frequency:
  ω_local = |θ_vel| / (|θ| + ε)

Neighboring voxels with matched frequencies reinforce (stable).
Mismatched frequencies → destructive interference → energy converts to phi.

This is a Kuramoto-like coupling: synchronized oscillators form coherent domains.
Domain boundaries (frequency mismatch) create density concentrations.

## New terms

```
For each neighbor j of voxel i:
  omega_i = |theta_vel_i| / (|theta_i| + eps)
  omega_j = |theta_vel_j| / (|theta_j| + eps)
  dw = |omega_i - omega_j|
  
  theta_acc[a] -= sigma * dw * theta[a]        (drain mismatched theta)
  phi_acc[a]   += sigma * dw * curl_theta[a]   (convert to phi)
```

The conversion strength is proportional to frequency mismatch (dw).
Zero conversion when perfectly synchronized.

## Parameters

- sigma: frequency-mismatch conversion strength (test: 0.05, 0.2, 0.5)
- eps: regularization to avoid divide-by-zero (fixed: 0.01)

## What to look for

- Frequency domains forming (regions of synchronized theta oscillation)
- Density concentrating at domain boundaries
- Different "particle" types = different natural frequencies
- Whether frequency-locked domains are truly stable

## Key difference from A and B

This doesn't depend on theta amplitude or gradient — it depends on temporal
behavior (how fast theta oscillates). Two voxels can have identical |θ| and
identical |∇θ| but different frequencies, and they'll interact differently.
This is genuinely new information not captured by spatial derivatives alone.

## Config

Base: m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
BC: absorb T=30 then periodic, T=600
