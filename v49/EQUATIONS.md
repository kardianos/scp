# V49 Field Equations — Unified Transfer Potential

## Fields

    phi_a(x,t)    a = 0,1,2    displacement (position sector)
    theta_a(x,t)  a = 0,1,2    rotation (angle/EM sector)

## Equations of Motion

### Phi equation

    d²phi_a/dt² = ∇²phi_a - m²phi_a - dVdP × dPda × f_phi + η × curl(θ)_a

### Theta equation

    d²theta_a/dt² = ∇²theta_a - m_θ²theta_a + theta_drive × theta_a + η × curl(φ)_a

## Definitions

### Triple product and potential

    P = phi_0 × phi_1 × phi_2

    dP/dphi_0 = phi_1 × phi_2                 (and cyclic)

    V_base(P) = (μ/2) × P² / (1 + κ × P²)    (NEGATIVE, μ < 0)

    |V_base|  = -V_base                        (POSITIVE)

    dVdP = μ × P / (1 + κ × P²)²              (NEGATIVE for P > 0)

### Transfer function

    Θ = theta_0² + theta_1² + theta_2²         (local theta energy density)

    f_transfer = ε + (1-ε) × Θ / (Θ + Θ_c)    (Michaelis-Menten proportion)

### Confinement function

    f_confine = γ × ln(1 + Θ / Θ_c)

### Phi force modulation

    f_phi = f_transfer - f_confine

The MINUS is because W_confine = |V_base| × f_confine has the opposite
P-derivative from W = V_base × f_transfer (since |V_base| = -V_base).

    dU/dP = dVdP × (f_transfer - f_confine) = dVdP × f_phi

Transfer ATTRACTS (strengthens binding). Confine REPELS (opposes binding).
This provides negative feedback: more confinement → weaker binding.

### Theta force

    df/dΘ = (1-ε) × Θ_c / (Θ + Θ_c)²         (transfer derivative)

    confine_deriv = γ / (Θ + Θ_c)              (confinement derivative)

    theta_drive = |V_base| × (df/dΘ - confine_deriv) × 2

Note: the theta force is +theta_drive × theta_a (PLUS sign).

    Small Θ: df/dΘ > confine_deriv → theta_drive > 0 → DRIVES theta growth
    Large Θ: confine_deriv > df/dΘ → theta_drive < 0 → CONFINES theta
    Equilibrium: df/dΘ = confine_deriv → theta_drive = 0

The equilibrium theta is:

    Θ_eq = Θ_c × ((1-ε)/γ - 1)

The initial growth rate ratio at Θ=0 is (1-ε)/γ. This must be close
to 1.0 to avoid tachyonic blowup (fast exponential growth of theta).

### Curl coupling

    curl(F)_0 = ∂F_2/∂y - ∂F_1/∂z            (and cyclic)

## Energy Functional

    E_total = E_phi_kin + E_theta_kin + E_grad_phi + E_grad_theta
            + E_mass_phi + E_mass_theta + E_transfer + E_confine + E_coupling

    E_phi_kin   = (1/2) Σ_a |∂phi_a/∂t|²
    E_theta_kin = (1/2) Σ_a |∂theta_a/∂t|²
    E_grad_phi  = (1/2) Σ_a |∇phi_a|²
    E_grad_theta= (1/2) Σ_a |∇theta_a|²
    E_mass_phi  = (1/2) m² Σ_a phi_a²
    E_mass_theta= (1/2) m_θ² Σ_a theta_a²
    E_transfer  = V_base(P) × f_transfer(Θ)                    (NEGATIVE, binding)
    E_confine   = |V_base(P)| × γ × ln(1 + Θ/Θ_c)             (POSITIVE, theta trap)
    E_coupling  = η × Σ_a phi_a × curl(theta)_a

## Constants

### Current best values (from sweep)

| Constant | Value | Type | Role |
|----------|-------|------|------|
| m² | 2.25 | forced | Phi mass (carrier frequency, confinement scale) |
| m_θ² | 0.0025 (m_θ=0.05) | **free** | Bare theta mass. Yukawa range 1/m_θ = 20. Gives theta a global decay — prevents box saturation with periodic BC. Acts as the "pion mass." |
| μ | -41.345 | free | Binding strength (NEGATIVE = attractive). From V28 CMA-ES. |
| κ | 50.0 | free | Potential saturation. Prevents P runaway. |
| η | 0.5 | free | Curl coupling (phi ↔ theta). EM coupling constant. |
| ε | 0.1 | free | Pilot light. 10% baseline binding when theta is absent. |
| Θ_c | 0.005 | **free** | Transfer activation threshold. Must be calibrated to match the natural theta amplitude at the braid core. Too large → transfer dormant. Too small → transfer always active. |
| γ | 0.85 | **free** | Confinement strength. Must be close to (1-ε) = 0.9 to avoid tachyonic blowup. The ratio (1-ε)/γ controls the initial theta growth rate. |
| δ | {0, 3.0005, 4.4325} | forced | Phase offsets (from V28 CMA-ES). |
| A_bg | 0.1 | forced | Background carrier amplitude. |

### Derived relationships

    Initial theta growth ratio = (1-ε)/γ = 0.9/0.85 = 1.06   (gentle)

    Theta equilibrium = Θ_c × ((1-ε)/γ - 1) = 0.005 × 0.059 = 0.000295
    → theta_rms_eq = sqrt(Θ_eq/3) ≈ 0.010

    Yukawa range (bare mass) = 1/m_θ = 20 code units
    Yukawa range (at core, from V_base) = depends on |V_base| and theta_drive

    Pilot light binding fraction = ε = 0.1 (10% of full V_base)
    Full binding fraction = 1.0 (when theta fully present)

### Boundary conditions

| Constant | Value | Role |
|----------|-------|------|
| bc_type | 0 | Absorbing sphere |
| damp_width | 3.0 | Boundary layer width |
| damp_rate | 0.01 | Damping coefficient |
| bc_switch_time | 50.0 | Switch absorbing → periodic at this time |

### Tuning notes

    γ must be close to (1-ε):
      γ = 0.1  → ratio 9.0  → instant tachyonic blowup
      γ = 0.5  → ratio 1.8  → fast growth, marginal
      γ = 0.85 → ratio 1.06 → gentle, stable
      γ = 0.9  → ratio 1.0  → no transfer (balanced)
      γ > 0.9  → ratio < 1  → theta suppressed (anti-transfer)

    Θ_c must match the natural theta level at the braid core:
      theta_rms at core ≈ 0.01-0.03 (from curl coupling + bare mass)
      Θ_core ≈ 3 × theta_rms² ≈ 0.0003 - 0.003
      Θ_c should be comparable: 0.001 - 0.01

      Θ_c = 1.0   → f_avg ≈ ε (transfer dormant)
      Θ_c = 0.05   → f_avg ≈ ε (still dormant, Θ_c >> Θ_core)
      Θ_c = 0.005  → f_avg ≈ 0.2 (transfer active)
      Θ_c = 0.001  → f_avg ≈ 0.3 (strongly active)
      Θ_c = 0.0001 → f_avg → 1-ε (saturated, too sensitive)

    m_θ controls theta spatial extent:
      m_θ = 0     → theta fills box (saturation)
      m_θ = 0.03  → range 33 (nearly fills L=15 box)
      m_θ = 0.05  → range 20 (moderate confinement)
      m_θ = 0.1   → range 10 (strong confinement, transfer suppressed)
      m_θ = 0.3   → range 3.3 (very tight, nuclear-only)

## What changed from V48

| Aspect | V48 | V49 |
|--------|-----|-----|
| Phi binding | V(P) = (μ/2)P²/(1+κP²) always at full strength | V_base × f_phi — modulated by theta presence |
| Theta mass | λ_θ P² (separate, additive) | Unified: from |V_base| × (df/dΘ - confine) + bare m_θ² |
| Phi back-reaction | -λ_θ P dPda |θ|² (separate term) | Implicit: -dVdP × dPda × f_phi (through f_phi) |
| Energy transfer | None (two independent mechanisms) | Unified: binding ↔ theta through single potential W |
| Parameters | μ, κ, η, λ_θ (4 free) | μ, κ, η, ε, Θ_c, γ, m_θ (7 free) |
| Theta in free space | Nearly massless (P²~10⁻⁸) | Mass m_θ² (bare, constant) |

## Open questions

1. Can the braid survive at full strength with the transfer active?
   Current runs show phi_max ~0.3-0.5 when transfer is on, vs 0.5-0.9
   in V48. The transfer weakens binding — is this too much?

2. What is the right Θ_c for a two-proton binding test? The single-proton
   theta level depends on N, L, and m_θ. Need to measure theta at the
   core vs background on the actual test grid before setting Θ_c.

3. Does the unified transfer produce a different inter-baryon force
   than the V48 λ_θ mechanism? The energy transfer (phi↔theta) should
   create an additional attractive channel between overlapping baryons.

4. The parameter count increased (4 → 7). Can any parameters be
   eliminated or derived from others? The relationships ε + γ ≈ 1
   and Θ_c ≈ Θ_core suggest constraints.
