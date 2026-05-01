# V54 Vacuum Pair Extension — Postulates

## Motivation

Virtual particle-antiparticle pairs in QFT follow a creation-annihilation loop:

    vacuum → (e⁻, e⁺) → vacuum,   constrained by ΔE·Δt ≥ ℏ/2

In a density gradient, the pair separates before annihilating, creating a net
polarization that feeds back on the source field. This is the physical basis
for vacuum polarization, the Schwinger effect, and the Casimir force.

The SCP system already exhibits one-way theta drain: η×curl(φ) sources θ, but
θ disperses without feeding back. The vacuum pair model provides energy-conserving
coupling terms that:
- Create θ (pairs) when φ (source field) exceeds a critical amplitude
- Self-limit θ growth through pair self-interaction
- Conserve total energy (Lagrangian-derived terms)

## Postulates

### P1: Cross-Potential (Pair Creation)

The vacuum pair density couples to the source field through a quartic interaction:

    V_cross = -(σ_cross/2) |θ|² |φ|²

where σ_cross > 0. This generates Euler-Lagrange forces:

    F_φ_a = +σ_cross · |θ|² · φ_a    (vacuum screening of φ)
    F_θ_a = +σ_cross · |φ|² · θ_a    (pair creation from φ)

Physical interpretation: The φ field (source) reduces the effective mass² of θ
(pair field). When |φ|² is large enough, θ's effective mass² becomes negative —
the vacuum becomes unstable and pairs are spontaneously created.

### P2: Self-Potential (Uncertainty Cutoff)

Virtual pairs interact with each other, imposing a self-limiting quartic potential:

    V_self = +(λ_self/4) |θ|⁴       (λ_self > 0 for self-limiting)

Euler-Lagrange force:

    F_θ_a = -λ_self · |θ|² · θ_a    (pair saturation)

Physical interpretation: As pair density grows, self-interaction increases the
effective mass of θ, halting exponential growth. This implements the uncertainty
principle's constraint: pairs cannot accumulate without limit because their
energy cost rises with density.

### P3: Uncertainty Constraint (Derived)

From P1 + P2, the equilibrium pair density is:

    |θ|²_eq = (σ_cross · |φ|² - m_θ²) / λ_self

This is analogous to ΔE·Δt = ℏ/2: the pair density (related to ΔE) times
the effective lifetime (related to λ_self) equals a constant set by σ_cross.

## Complete Equations

    d²φ_a/dt² = ∇²φ_a - m²φ_a - ∂V/∂φ_a + η·curl(θ)_a + σ_cross·|θ|²·φ_a

    d²θ_a/dt² = ∇²θ_a - m_θ²θ_a + η·curl(φ)_a + σ_cross·|φ|²·θ_a - λ_self·|θ|²·θ_a

where V(P) = (μ/2)P²/(1+κP²) is the standard saturation potential.

## Linear Stability Analysis

### Setup

Background: φ_a = A·cos(k₀z + δ_a), θ_a = 0
Perturbation: δθ_a ~ exp(i(q·x - ωt))

### Dispersion Relation for θ Modes

    ω² = |q|² + m_θ² - σ_cross · ⟨|φ|²⟩

where ⟨|φ|²⟩ = (3/2)A² (time-averaged over the 3-field background).

### Stability Criterion

Vacuum is STABLE when ω² > 0 for all q:

    σ_cross · (3/2)A² < m_θ²

Vacuum becomes UNSTABLE (pair creation) when:

    A > A_crit = sqrt(2 m_θ² / (3 σ_cross))

### Special Cases

- m_θ = 0: A_crit = 0. ANY nonzero φ creates pairs. This explains the
  unconditional theta growth in all prior SCP simulations with massless θ.

- m_θ > 0, σ_cross = 0: Standard massive theta, no pair creation.
  Theta oscillates independently. (Matches V54 m_theta test result.)

### Growth Rate

Below cutoff (q < q_max): exponential growth rate:

    γ(q) = sqrt(σ_cross · (3/2)A² - m_θ² - |q|²)

Maximum growth rate (at q = 0):

    γ_max = sqrt(σ_cross · (3/2)A² - m_θ²)

Growth time: τ_growth = 1/γ_max

### Nonlinear Saturation

When θ grows, the λ_self term kicks in. Equilibrium:

    |θ|_sat = sqrt((σ_cross · ⟨|φ|²⟩ - m_θ²) / λ_self)

For m_θ = 0:

    |θ|_sat / |φ| = sqrt(σ_cross / λ_self) · sqrt(3/2)

The ratio σ_cross/λ_self controls the equilibrium θ/φ ratio.

## Test Plan

### Test 1: Stable Vacuum (A < A_crit)
- Set m_theta = 1.0 (m_θ² = 1.0), sigma_cross = 0.5, lambda_self = 1.0
- A_crit = sqrt(2·1.0 / (3·0.5)) = 1.155
- Use A = 0.8 (below threshold)
- Predict: theta stays near zero, no pair creation
- Measure: theta_rms vs time — should remain < 0.01

### Test 2: Unstable Vacuum (A > A_crit)
- Same parameters but A = 1.5 (above threshold)
- γ_max = sqrt(0.5·(3/2)·1.5² - 1.0) = sqrt(0.6875) = 0.829
- τ_growth = 1/0.829 = 1.21 time units
- |θ|_sat = sqrt((0.5·(3/2)·1.5² - 1.0)/1.0) = sqrt(0.6875) = 0.829
- Predict: theta grows exponentially at rate 0.83, saturates at ~0.83
- Measure: theta_rms vs time — exponential rise then plateau

### Test 3: Marginal Case (A ≈ A_crit)
- A = 1.15 (just below threshold)
- Predict: slow oscillation of theta, no sustained growth

## Parameters

| Parameter | Symbol | Default | Role |
|-----------|--------|---------|------|
| sigma_cross | σ | 0.0 | Cross-potential coupling |
| lambda_self | λ | 0.0 | Theta self-potential |
| m_theta | m_θ | 0.0 | Theta mass (pair creation threshold) |
| eta | η | 0.5 | Curl coupling (unchanged) |

## Energy Conservation

Both V_cross and V_self are Lagrangian-derived. The total energy:

    E = E_kin + E_grad + E_mass + E_pot + E_coupling
        + (σ_cross/2)|θ|²|φ|²    [cross energy, NEGATIVE contribution]
        + (λ_self/4)|θ|⁴          [self energy, POSITIVE contribution]

must be conserved by the Verlet integrator to machine precision.
