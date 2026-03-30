# V49 Variant B (Final): Unified Transfer Potential with Pilot Light

## The Potential

    W(P, Θ) = (mu/2) × P² × [epsilon + (1 - epsilon) × Θ / (Θ + Θ_c)]

Where:
    P  = phi_0 × phi_1 × phi_2       (triple product)
    Θ  = |theta|² = theta_0² + theta_1² + theta_2²
    Θ_c = theta_c²                    (transfer midpoint, new parameter)
    epsilon = 0.01                    (pilot light, hardcoded or config)

Note: using Θ (not Θ²) in the sigmoid to avoid the ghost energy problem.
With Θ/(Θ + Θ_c), the force on theta is:

    dW/dtheta_a = (mu/2) × P² × (1 - epsilon) × Θ_c × 2 theta_a / (Θ + Θ_c)²

This never goes to zero — at large Θ it decays as 1/Θ, not as 1/Θ².
Energy can always flow BACK from theta to the binding sector.

## Behavior

    Θ = 0:     W = (mu/2) × P² × epsilon ≈ 0.01 × V_max
               Pilot light: 1% binding active, enough to seed braids.

    Θ = Θ_c:   W = (mu/2) × P² × [epsilon + (1-epsilon)/2] ≈ 0.5 × V_max
               Half-transfer point. Binding at 50% strength.
               Force on theta is at maximum: energy flowing in.

    Θ >> Θ_c:  W → (mu/2) × P² × 1.0 = V_max
               Full binding. Theta has absorbed its energy share.
               Force on theta decays as 1/Θ (not zero — no dead zone).

The S-curve: binding_fraction = epsilon + (1 - epsilon) × Θ / (Θ + Θ_c)

## Equations of Motion

### Phi equation (from dL/dphi_a):

    d²phi_a/dt² = nabla²phi_a - m²phi_a
                  - [mu × P / (binding_sigmoid)] × dP/dphi_a
                  + eta × curl(theta)_a

Where binding_sigmoid = epsilon + (1 - epsilon) × Θ / (Θ + Θ_c)

The binding force dW/dP = mu × P × binding_sigmoid is:
- At Θ=0: mu × P × 0.01 (weak pilot light)
- At Θ=Θ_c: mu × P × 0.5 (half strength)
- At Θ→∞: mu × P × 1.0 (full strength, same as current V(P) at low P)

NOTE: The kappa saturation is REMOVED. The transfer sigmoid replaces it.
At high P, the binding force grows linearly (mu × P), but it's gated by
how much theta is present. The "saturation" is now physically meaningful:
binding can only be as strong as the theta excitation allows.

### Theta equation (from dL/dtheta_a):

    d²theta_a/dt² = nabla²theta_a
                    - [mu × P² × (1 - epsilon) × Θ_c / (Θ + Θ_c)²] × theta_a
                    + eta × curl(phi)_a

The theta "mass" term is now DERIVED from the unified potential, not
a separate lambda_theta parameter. It depends on P² (only active at
braid cores) AND on Θ (weakens as theta grows, preventing runaway).

At braid core (P ~ 0.08), this gives an effective theta mass of:

    m_theta_eff² = mu × P² × (1 - epsilon) × Θ_c / (Θ + Θ_c)²

The mass is largest at Θ = 0 (m_eff² = mu × P² × (1-eps) / Θ_c)
and decreases as theta is excited. This is the "transfer" — as theta
absorbs energy, the mass drops, allowing theta to propagate more freely.

### NO separate V(P) saturation. NO separate lambda_theta.

The old equations had:
    V(P) = (mu/2) P² / (1 + kappa P²)     ← kappa saturation
    m_theta² = lambda_theta × P²           ← additive coupling

The new equation has:
    W(P, Θ) = (mu/2) P² × sigmoid(Θ)      ← ONE unified term

Parameters eliminated: kappa, lambda_theta
Parameters added: Θ_c, epsilon
Net change: 3 params → 2 params (reduction)

## Parameter Calibration

### Θ_c (transfer midpoint)

At equilibrium, the braid should have Θ ≈ Θ_c (operating at the
sigmoid midpoint for maximum energy transfer). From V48 data:

    theta_rms at braid core ~ 0.05 (V34 hires)
    |theta|² at core ~ 3 × 0.05² = 0.0075

So Θ_c ≈ 0.01 is a starting estimate.

### mu (binding strength)

With the kappa saturation removed, mu needs recalibration. The current
mu = -41.345 was tuned with kappa = 50 to give a peak force at P ~ 0.08.
Without kappa, the force is mu × P × sigmoid — linear in P.

For the sigmoid to give similar effective binding at the equilibrium P:

    Old peak force: mu × P / (1+kappa P²)² at P=1/sqrt(3*kappa) ≈ 0.082
                  = -41.345 × 0.082 / (1 + 50×0.0067)² = -2.54

    New force: mu_new × P × sigmoid(Θ_c) = mu_new × 0.082 × 0.505 = -2.54
    mu_new = -2.54 / 0.041 ≈ -61

But this should be determined by running the simulation, not by matching
to the old theory. Start with mu = -41.345, Θ_c = 0.01, and adjust.

### epsilon (pilot light)

0.01 is the default. This gives 1% baseline binding — enough to seed
structure formation from pure background, weak enough to not affect
established braids (which have sigmoid ≈ 0.5-1.0).

## Energy Functional

    E_total = (1/2)|dphi/dt|² + (1/2)|nabla phi|² + (1/2)m²|phi|²
            + (1/2)|dtheta/dt|² + (1/2)|nabla theta|²
            + W(P, Θ)
            + eta × phi · curl(theta)

W(P, Θ) replaces BOTH E_pot and E_tmass from the old decomposition.
In diagnostics, report:

    E_transfer = W(P, Θ)    (the unified transfer energy)
    sigmoid_avg = <binding_sigmoid> over grid  (measures transfer state)

## Implementation

In compute_forces, replace the V(P) + lambda_theta computation:

```c
double P = p0*p1*p2;
double t0=theta[0][idx], t1=theta[1][idx], t2=theta[2][idx];
double Theta = t0*t0 + t1*t1 + t2*t2;

// Unified transfer sigmoid
double sigmoid = EPSILON + (1.0 - EPSILON) * Theta / (Theta + THETA_C);

// Force on phi: binding gated by sigmoid
double dWdP = MU * P * sigmoid;
for (int a = 0; a < 3; a++) {
    double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
    double lap = laplacian(phi[a], ...);
    double ct = curl(theta, a, ...);
    phi_acc[a] = lap - MASS2*phi[a] - dWdP*dPda + ETA*ct;
}

// Force on theta: energy transfer from binding
double P2 = P*P;
double denom = (Theta + THETA_C) * (Theta + THETA_C);
double transfer_coeff = MU * P2 * (1.0 - EPSILON) * THETA_C / denom;
for (int a = 0; a < 3; a++) {
    double lapt = laplacian(theta[a], ...);
    double cp = curl(phi, a, ...);
    theta_acc[a] = lapt - transfer_coeff * theta[a] + ETA*cp;
}
```

Config parameters:
    theta_c = 0.01       (transfer midpoint)
    epsilon = 0.01       (pilot light)
    mu = -41.345         (binding strength, may need recalibration)
    kappa: REMOVED
    lambda_theta: REMOVED

## Test Plan

1. Single proton (N=64, T=100): verify stability, measure equilibrium Θ
2. Compare energy partition vs old V(P): does theta absorb ~10% as before?
3. Two protons (D=14, N=384, T=1000): does unified transfer give binding?
4. Parameter sweep on Θ_c: {0.001, 0.005, 0.01, 0.05, 0.1}
5. Check: does removing kappa cause blowup at high P? The sigmoid
   gates the force but doesn't cap P itself. May need P-dependent
   epsilon or a weak residual kappa as safety valve.

## Risk: Unbounded P

Without kappa, V(P) = (mu/2) P² × sigmoid can grow without bound.
The sigmoid caps at 1.0, so V → (mu/2) P² at high theta. This is
the unsaturated potential that caused blowups in the original theory.

Mitigation options:
a) Keep a weak kappa as safety: W = (mu/2) P²/(1+kappa_safety P²) × sigmoid
   with kappa_safety small (e.g., 5 instead of 50). The sigmoid does the
   physics; kappa_safety prevents numerical blowup.
b) Trust the sigmoid: at high P, the force drives theta excitation,
   which increases the sigmoid, which increases the binding... this
   could still run away. Need to verify numerically.
c) Use the existing absorbing BC + bc_switch_time to damp transients,
   then check if the system is self-regulating.

Recommend starting with option (a) — keep kappa=5 as a safety valve.
If the sigmoid provides sufficient self-regulation, reduce kappa toward 0.
