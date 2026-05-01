# Test E: Theta Self-Potential Only

## Lagrangian term
L_self = -(lambda/4) |θ|⁴

## Forces (energy-conserving by construction)
theta_acc[a] -= lambda * |θ|² * θ[a]

No phi modification — energy stays in theta sector.

## Parameters
sigma_cross = 0, lambda_self = 5

## Expected
- m_eff(phi) = 3.00 (unchanged)
- m_eff(theta) = 0.67 (weak confinement from self-interaction)
- Theta caps but doesn't co-localize with phi
- Like test B (cubic) but from a proper Lagrangian

## Config
m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
sigma_cross=0, lambda_self=5
BC: absorb T=0-30, periodic T=30-600
