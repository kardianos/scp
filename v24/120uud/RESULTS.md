# V24-UUD: Proton-Like 120-degree Oscillon Results

## Summary

**NEGATIVE RESULT.** The UUD 120-degree oscillon with per-field masses does NOT survive as a long-lived bound state at the tested parameters. The 120-degree phase structure is lost within a few hundred time units, and energy dissipates below 35% retention in all cases. The fundamental obstruction is that the pairwise coupling lambda reduces the effective mass gap for the 120-degree mode, pushing the oscillation frequency above the radiation threshold.

---

## Model

Three massive scalars with triple-product coupling and pairwise mass mixing:

    ddot(phi_a) = lapl(phi_a) - m_a^2 phi_a - lambda*(phi_b + phi_c)
                  - dV_triple/dphi_a

    V_triple = (mu/2) P^2 / (1 + kappa P^2),  P = phi_1 phi_2 phi_3

Fields 1,2 are "up" (mass m_U), field 3 is "down" (mass m_D).

The mass matrix M^2_{ab} = m_a^2 delta_{ab} + lambda (1-delta_{ab}) has eigenvalues:
- Symmetric mode (1,1,1): m^2 + 2*lambda  (heaviest)
- Antisymmetric 120-degree modes: m^2 - lambda  (lightest)

The 120-degree mode effective mass is m_eff = sqrt(m^2 - lambda). An oscillon at frequency omega survives only if omega < m_eff (subgap condition).

## Protocol

1. **Equilibrate** symmetric (0-degree) oscillon with triple-product coupling only (lambda=0, equal masses m_U). This establishes a stable localized oscillon with well-defined frequency omega_eq.
2. **Rotate** phases to 120 degrees (0, 2pi/3, 4pi/3), change field 3 mass to m_D, enable pairwise coupling lambda.
3. **Evolve** for t_run=5000, recording energy, phase differences, core fraction.

Parameters: Nx=4000, xmax=100, dx=0.05, t_equil=3000, t_run=5000.

## Key Physics Constraint

For the 120-degree oscillon to survive, the equilibrated frequency must satisfy:

    omega_eq < sqrt(m_D^2 - lambda)

This constrains the parameter space:

| lambda | m_eff(120) | omega_eq (mu=-10,k=10) | omega_eq (mu=-20,k=20) | Stable? |
|--------|-----------|------------------------|------------------------|---------|
| 0.85   | 0.387     | 0.822                  | 0.858                  | NO (omega >> m_eff) |
| 0.30   | 0.837     | 0.822                  | 0.858                  | MARGINAL (omega ~ m_eff) |
| 0.10   | 0.949     | 0.822                  | 0.858                  | YES in principle |
| 0.05   | 0.975     | 0.822                  | 0.858                  | YES in principle |

The proposal parameter lambda=0.85 is incompatible with the equilibrated omega. The 120-degree effective mass gap (0.387) is far below the oscillation frequency (0.82-0.86).

## Results: lambda=0.85 Scan

With the proposal parameters (mu=-20, kappa=20, lambda=0.85):

The pairwise coupling makes the total energy unbounded from below for non-localized 120-degree configurations. After phase rotation, the fields spread across the domain, the triple-product potential (which saturates at |mu|/(2*kappa) per unit length) becomes dominant, and the total energy goes deeply negative (E ~ -36 to -48).

**Conclusion**: lambda=0.85 is incompatible with stable 120-degree oscillons.

## Results: lambda=0.05 Scan (mu=-10, kappa=10)

Equilibrated at omega_eq=0.822 < m_eff(120)=0.975 (subgap condition satisfied).

| m_D  | m_D/m_U | E_init | E_final | E_retained | fc_final | omega_1 | omega_3 | Survived |
|------|---------|--------|---------|------------|----------|---------|---------|----------|
| 1.00 | 1.000   | 3.402  | 0.666   | 19.6%      | 0.238    | 0.978   | 0.978   | NO       |
| 0.98 | 0.980   | 3.397  | 0.690   | 20.3%      | 0.223    | 0.960   | 0.960   | NO       |
| 0.96 | 0.960   | 3.392  | 0.865   | 25.5%      | 0.220    | 0.978   | 0.942   | NO       |
| 0.94 | 0.940   | 3.387  | 1.030   | 30.4%      | 0.212    | 0.978   | 0.924   | NO       |
| 0.92 | 0.920   | 3.382  | 1.129   | 33.4%      | 0.189    | 0.978   | 0.906   | NO       |
| 0.90 | 0.900   | 3.378  | 1.054   | 31.2%      | 0.162    | 0.978   | 0.888   | NO       |
| 0.85 | 0.850   | 3.366  | 0.880   | 26.1%      | 0.149    | 0.978   | 0.840   | NO       |

### Observations

1. **No survivors at any m_D.** Even the equal-mass case (m_D=1.0) fails because the post-rotation transient pushes the oscillation frequency (0.978) near the mass gap (m=1.0).

2. **Frequency splitting confirmed.** omega_1 (up quark) stays at 0.978 while omega_3 (down quark) decreases linearly with m_D: omega_3 approximately equals m_D. This is the expected mass-splitting signature.

3. **Optimal m_D ~ 0.92.** Energy retention peaks at m_D=0.92 (33.4%). The heavier down quark (smaller m_D) creates a lower radiation threshold for field 3, but the lighter mass also reduces the binding.

4. **120-degree phase NOT maintained.** Phase differences drift away from +/-120 degrees within ~200 time units after rotation. The triple-product coupling is not strong enough to enforce the phase lock.

## Root Cause Analysis

The 120-degree oscillon faces a triple obstruction:

1. **Weak binding in 1D.** The triple-product binding P = phi_1*phi_2*phi_3 for 120-degree phases produces P oscillating at 3*omega (not omega). Since 3*omega=2.47 >> m=1.0, the coupling is averaged away over one mass period.

2. **Phase rotation transient.** Rotating from 0-degree to 120-degree injects energy into the antisymmetric modes. This transient causes ~50% energy loss in the first 500 time units.

3. **Pairwise coupling conflict.** The pairwise coupling favors the 0-degree (symmetric) configuration energetically. The 120-degree state has HIGHER pairwise energy than the 0-degree state for the same amplitude. The system naturally relaxes toward 0-degree.

## Comparison with v24/phase120

The v24/phase120 code (without pairwise coupling, equal masses) showed the same result at mu=-20, kappa=20: the 120-degree state retains only 17% of energy vs 99% for the 0-degree control. The 120-degree configuration is inherently unstable at these parameters because 3*omega >> m.

For 120-degree oscillons to work, need 3*omega < m, i.e., omega < m/3 = 0.333. This requires much stronger coupling than mu=-10 to -30.

## Files

| File | Description |
|------|-------------|
| `src/uud120.c` | Solver: per-field masses + pairwise coupling + triple product |
| `data/uud_summary.tsv` | Energy vs m_D scan results |
| `data/uud_mD{X}_ts.tsv` | Time series for each m_D value |
| `data/uud_mD{X}_phi{N}_spectrum.tsv` | DFT power spectra per field |

## Recommendations

To achieve a stable UUD 120-degree oscillon, would need:
- omega < m/3 = 0.333, requiring much stronger triple-product coupling (mu << -30)
- OR: a different coupling mechanism that directly enforces 120-degree phases
- OR: higher spatial dimension (3D), where the nonlinear shift is different
- The pairwise lambda must satisfy lambda < m_D^2 - omega^2 for the lightest field
