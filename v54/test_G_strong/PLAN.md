# Test G: Strong Combined

## Lagrangian terms
L = -(sigma/2)|θ|²|φ|² - (lambda/4)|θ|⁴

## Forces (energy-conserving by construction)
phi_acc[a]   -= sigma * |θ|² * φ[a]
theta_acc[a] -= (lambda * |θ|² + sigma * |φ|²) * θ[a]

## Parameters
sigma_cross = 20, lambda_self = 10

## Expected
- m_eff(phi) = 3.29 (10% tighter)
- m_eff(theta) = 6.77 (very strongly confined)
- Tightest co-localization of all tests
- If any test shows persistent structure, it should be this one

## Config
m=3, mu=-90, kappa=50, eta=0.3, A=1.5, A_bg=0.1
sigma_cross=20, lambda_self=10
BC: absorb T=0-30, periodic T=30-600
