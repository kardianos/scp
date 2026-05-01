# Test F: Combined Cross + Self

## Lagrangian terms
L = -(sigma/2)|θ|²|φ|² - (lambda/4)|θ|⁴

## Forces (energy-conserving by construction)
phi_acc[a]   -= sigma * |θ|² * φ[a]
theta_acc[a] -= (lambda * |θ|² + sigma * |φ|²) * θ[a]

## Parameters
sigma_cross = 10, lambda_self = 5

## Expected
- m_eff(phi) = 3.15 (5% tighter)
- m_eff(theta) = 4.79 (strongly confined)
- Theta capped AND co-localized with phi
- Best of both: cap prevents runaway, cross-coupling prevents dispersal

## Config
m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
sigma_cross=10, lambda_self=5
BC: absorb T=0-30, periodic T=30-600
