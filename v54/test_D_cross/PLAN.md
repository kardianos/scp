# Test D: Cross-Potential Only

## Lagrangian term
L_cross = -(sigma/2) |θ|² |φ|²

## Forces (energy-conserving by construction)
phi_acc[a]   -= sigma * |θ|² * φ[a]
theta_acc[a] -= sigma * |φ|² * θ[a]

## Parameters
sigma_cross = 10, lambda_self = 0

## Expected
- m_eff(phi) = 3.15 (5% tighter than baseline)
- m_eff(theta) = 4.74 (from massless to confined!)
- Theta and phi co-localize: each confines the other
- Energy E_cross = (sigma/2)|θ|²|φ|² conserved

## Config
m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
sigma_cross=10, lambda_self=0
BC: absorb T=0-30, periodic T=30-600
