# V49 Variant B v2: Unified Transfer with Linear Proportion

**STATUS: SIGN ERROR FIXED in v2.1 (see below)**

## Corrections from Review (v2)

1. Safety kappa KEPT (prevents P-collapse)
2. Theta force fixed (ln form for confinement)
3. No sigmoid — two options: linear proportion or hard transition

## Corrections from Review (v2.1)

4. CRITICAL: V_base is NEGATIVE (mu < 0). All theta forces must use
   |V_base| = -V_base, not V_base directly. Using V_base produces
   anti-confinement (theta blows up instantly). Fixed throughout.
5. Missing phi back-reaction from W_confine. The ln term also depends
   on P through V_base, so phi must feel it too. Fixed in implementation.
6. Θ_c increased to avoid parametric oscillation of effective mass.

## The Base Potential (unchanged, keeps braids alive)

    V_base(P) = (mu/2) × P² / (1 + kappa × P²)

This is the existing V(P). It's always present. Braids survive from this
alone, with or without theta. kappa prevents P runaway.

## The Transfer Function

The transfer function f(Θ) controls how much of the binding is "active."
It replaces the old lambda_theta mechanism. Two options:

### Option 1: Linear Proportion (recommended)

    f(Θ) = epsilon + (1 - epsilon) × min(Θ / Θ_c, 1.0)

    Θ < Θ_c:  f grows linearly with theta.
              Each unit of theta energy adds proportional binding.
              f = epsilon at Θ=0, f = 0.5+epsilon/2 at Θ=Θ_c/2.

    Θ >= Θ_c: f = 1.0. Full binding. Theta has done its job.

This is a simple ramp. "If something twists to a degree, something
else gives" — the binding is directly proportional to how much theta
is present, up to a cap.

The kink at Θ = Θ_c needs smoothing for the Verlet integrator.
Use a small blend zone of width delta:

    For Θ in [Θ_c - delta, Θ_c + delta]:
        f = smoothstep blend from ramp to flat

Or simply use the C¹ approximation:

    f(Θ) = epsilon + (1 - epsilon) × Θ / (Θ + Θ_c)

This is NOT a sigmoid — it's the Michaelis-Menten / Monod function,
which IS a linear proportion at small Θ (f ≈ Θ/Θ_c) that smoothly
caps at 1.0. It's the continuous version of the ramp with no kink.
Every engineer who's modeled saturation kinetics uses this form.

At Θ << Θ_c: f ≈ epsilon + Θ/Θ_c  (linear proportion)
At Θ = Θ_c:  f = epsilon + (1-epsilon)/2 ≈ 0.5
At Θ >> Θ_c: f → 1.0 (cap)

### Option 2: Hard Transition

    f(Θ) = epsilon            if Θ < Θ_threshold
    f(Θ) = 1.0               if Θ >= Θ_threshold

Binary: binding is either pilot-light or full. The transition is sharp.

For numerical stability, smooth over a narrow band (width w):

    f(Θ) = epsilon + (1-epsilon) × smoothstep((Θ - Θ_threshold) / w)

    smoothstep(x) = 0 for x<0, 3x²-2x³ for 0<=x<=1, 1 for x>1

This gives a C¹ transition over width w. Set w = 0.1 × Θ_threshold
for a fast but smooth switch.

The hard transition is cleaner conceptually: either the braid has
excited enough theta to "turn on" full binding, or it hasn't. No
partial credit. But it creates a binary world — braids are either
fully bound or barely bound, with little in between.

## The Unified Potential

    W(P, Θ) = V_base(P) × f(Θ)

    V_base(P) = (mu/2) × P² / (1 + kappa × P²)

    f(Θ) = option 1 or option 2 above

## Equations of Motion

### Phi equation

    d²phi_a/dt² = nabla²phi_a - m²phi_a
                  - dU/dphi_a
                  + eta × curl(theta)_a

The total potential is U = W + W_confine, and BOTH contribute to the
phi force. CRITICAL: W and W_confine have OPPOSITE P-derivatives because
W uses V_base (negative) while W_confine uses |V_base| (positive).

    dU/dP = dV_base/dP × f_transfer + d|V_base|/dP × gamma × ln(1+Θ/Θ_c)

Since |V_base| = -V_base: d|V_base|/dP = -dV_base/dP = -dVdP.

    dU/dP = dVdP × f_transfer + (-dVdP) × gamma × ln(...)
          = dVdP × (f_transfer - f_confine)

    dU/dphi_a = dVdP × (f_transfer - f_confine) × dP/dphi_a

    F_phi_a = -dU/dphi_a = -dVdP × (f_transfer - f_confine) × dP/dphi_a

    dV_base/dP = mu × P / (1 + kappa × P²)²    (dVdP, NEGATIVE for P>0)

SIGN CHECK: dVdP < 0, so -dVdP > 0.
  If f_transfer > f_confine: force is attractive (toward higher P) ✓
  If f_transfer < f_confine: force is repulsive (confinement cost dominates) ✓

The confinement term OPPOSES binding in the phi sector. This provides
negative feedback: more confinement → weaker phi binding → self-limiting.
WITHOUT this sign, the confinement ADDS to binding, creating runaway.

### Theta equation

**CRITICAL SIGN NOTE:** V_base(P) is NEGATIVE (mu = -41.345 < 0).
The theta "mass" must be POSITIVE to be confining. Therefore we use
|V_base| = -V_base in all theta force terms.

    d²theta_a/dt² = nabla²theta_a
                    - |V_base(P)| × [df/dΘ + gamma/(Θ+Θ_c)] × 2 × theta_a
                    + eta × curl(phi)_a

The effective theta mass squared at each voxel is:

    m_theta_local² = |V_base(P)| × [(1-ε)Θ_c/(Θ+Θ_c)² + gamma/(Θ+Θ_c)] × 2

This is ALWAYS POSITIVE (|V_base| >= 0, all other factors >= 0).
It provides a restoring force: F = -m_local² × theta_a (toward zero).

The two contributions:

The theta force from U = W + W_confine is F = -dU/dtheta_a:

    dU/dΘ = V_base × df/dΘ + |V_base| × gamma/(Θ+Θ_c)
           = -|V_base| × df/dΘ + |V_base| × gamma/(Θ+Θ_c)
           = |V_base| × [gamma/(Θ+Θ_c) - df/dΘ]

    F_theta_a = -dU/dtheta_a = -dU/dΘ × 2θ_a
              = |V_base| × [df/dΘ - gamma/(Θ+Θ_c)] × 2θ_a

SIGN CHECK: this is +|V_base| × (df - confine) × 2θ_a. Note the PLUS.

    Small Θ: df/dΘ > gamma/(Θ+Θ_c)
              → force in SAME direction as theta → DRIVES theta growth
              This IS the energy transfer: binding energy → theta excitation.

    Large Θ: gamma/(Θ+Θ_c) > df/dΘ
              → force OPPOSITE to theta → CONFINES theta
              Prevents escape to infinity.

    Equilibrium: df/dΘ = gamma/(Θ+Θ_c) → force is zero.
              Theta settles at this value. Stable balance.

The unified potential creates a POTENTIAL WELL AT FINITE THETA, not at
theta=0. It drives theta toward an equilibrium and holds it there.

### Theta Confinement (ln potential)

The confinement comes from:

    W_confine(P, Θ) = |V_base(P)| × gamma × ln(1 + Θ / Θ_c)

Note: |V_base|, NOT V_base. Since |V_base| >= 0, W_confine >= 0.
The force -dW_confine/dθ_a points TOWARD zero (confining).

## Full Energy Functional

    E_total = (1/2)|dphi/dt|² + (1/2)|nabla phi|² + (1/2)m²|phi|²
            + (1/2)|dtheta/dt|² + (1/2)|nabla theta|²
            + W(P, Θ)
            + W_confine(P, Θ)
            + eta × phi · curl(theta)

    W(P, Θ) = V_base(P) × f(Θ)                            (transfer, negative)
    W_confine(P, Θ) = |V_base(P)| × gamma × ln(1 + Θ/Θ_c) (confinement, positive)

Note: W is negative (attractive binding). W_confine is positive (theta
mass penalty). Their sum W + W_confine is the total potential.

For diagnostics, report:
    E_transfer = W(P, Θ)         (replaces old E_pot)
    E_confine = W_confine(P, Θ)  (new, measures theta trapping)
    f_avg = <f(Θ)> over grid     (measures how "transferred" the system is)

## Parameters

| Parameter | Value | Role | Replaces |
|-----------|-------|------|----------|
| mu | -41.345 | Binding strength (NEGATIVE) | same |
| kappa | 50 (or lower) | P saturation (safety) | same |
| Θ_c | 0.05 | Transfer midpoint (increased from 0.01 to avoid parametric noise) | lambda_theta |
| epsilon | 0.1 | Pilot light (baseline binding fraction) | new |
| gamma | 0.1 | Theta confinement strength (ln trap) | new |
| eta | 0.5 | Curl coupling | same |
| m² | 2.25 | Phi mass | same |
| m_theta² | 0.0 | Bare theta mass (can stay zero — confinement from W_confine) | same |

Net change: lambda_theta removed, kappa kept, added Θ_c, epsilon, gamma.
Total free physics params: was 6 (mu, kappa, eta, m_theta², lambda_theta, eta1),
now 6 (mu, kappa, eta, Θ_c, epsilon, gamma). Same count but cleaner physics.

SIGN CONVENTION: V_base < 0 always (mu < 0). The code must use
|V_base| = -V_base for all theta-related forces to ensure confinement.
The phi force uses dVdP directly (which is also negative, and the
double-negative gives the correct attractive force).

## Implementation

```c
// In compute_forces, per voxel:

double P = p0*p1*p2, P2 = P*P;
double den = 1.0 + KAPPA*P2;
double V_base = (MU/2.0) * P2 / den;          // existing V(P), NEGATIVE (mu<0)
double V_abs = -V_base;                         // |V_base|, POSITIVE
double dVdP = MU * P / (den*den);              // existing dV/dP, NEGATIVE

double t0=theta[0][idx], t1=theta[1][idx], t2=theta[2][idx];
double Theta = t0*t0 + t1*t1 + t2*t2;
double Theta_sum = Theta + THETA_C;

// Transfer function (linear proportion, smooth):
double f_transfer = EPSILON + (1.0-EPSILON) * Theta / Theta_sum;

// Confinement function (ln):
double f_confine = GAMMA * log(1.0 + Theta / THETA_C);

// Phi force modulation: transfer ATTRACTS, confinement REPELS.
// W = V_base * f_transfer  →  dW/dP = dVdP * f_transfer
// W_confine = |V_base| * f_confine  →  dWc/dP = -dVdP * f_confine
// Total: dU/dP = dVdP * (f_transfer - f_confine)
// Force: -dVdP * (f_transfer - f_confine) * dPda
//
// The MINUS on f_confine is because |V_base| = -V_base has the
// opposite P-derivative. Confinement COSTS energy when P grows,
// providing negative feedback that prevents runaway.
double f_phi = f_transfer - f_confine;

for (int a = 0; a < 3; a++) {
    double dPda = (a==0)?p1*p2:(a==1)?p0*p2:p0*p1;
    double lap = laplacian(phi[a], ...);
    double ct = curl(theta, a, ...);
    phi_acc[a] = lap - MASS2*phi[a] - dVdP*dPda*f_phi + ETA*ct;
}

// Theta force: F = +|V_base| × (df/dΘ - confine) × 2 × theta_a
// Small Θ: df > confine → DRIVES theta growth (transfer)
// Large Θ: confine > df → CONFINES theta (pulls back)
// Equilibrium where df = confine.
double df_dTheta = (1.0-EPSILON) * THETA_C / (Theta_sum*Theta_sum);
double confine_deriv = GAMMA / Theta_sum;
double theta_drive = V_abs * (df_dTheta - confine_deriv) * 2.0;
// theta_drive > 0 at small Θ (growth), < 0 at large Θ (confinement)

for (int a = 0; a < 3; a++) {
    double lapt = laplacian(theta[a], ...);
    double cp = curl(phi, a, ...);
    theta_acc[a] = lapt + theta_drive * theta[a] + ETA*cp;
}
```

### Sign verification:

    V_base < 0, V_abs > 0, dVdP < 0 (for P > 0, mu < 0)

    Phi force: -dVdP*(f_transfer - f_confine)*dPda
        If f_transfer > f_confine: = |dVdP| × positive × dPda → ATTRACTIVE ✓
        If f_confine > f_transfer: = |dVdP| × negative × dPda → REPULSIVE ✓
        Confinement opposes binding → negative feedback ✓

    Theta force: +theta_drive*theta_a = +V_abs*(df-confine)*theta_a
        Small Θ (df > confine): force SAME direction as theta → DRIVES growth ✓
        Large Θ (confine > df): force OPPOSITE to theta → CONFINES ✓
        Equilibrium at df = confine → stable balance ✓

    If theta escapes: confine_deriv dominates → pulls back → CORRECT ✓
    If theta collapses: df_dTheta dominates → drives back up → CORRECT ✓
    If theta escapes: ln confine term grows → stronger pull back → CORRECT

## What This Changes Physically

### Old (V48): Two separate mechanisms
- V(P) provides binding, saturates at high P (dead zone)
- lambda_theta × P² provides theta mass (additive, independent)
- Energy can pool in V(P) dead zone with no transfer to theta

### New (V49 B v2): One mechanism with proportional transfer
- V_base(P) provides binding, still saturates (safety)
- BUT the binding strength is proportional to theta excitation
- Theta is driven by the binding (and confined by the ln term)
- No dead zone: at any P, energy flows to theta proportionally
- The phi and theta sectors are COUPLED through one potential

### What stays the same
- The curl coupling (eta × curl) is unchanged
- The phi mass term (m²phi) is unchanged
- The Verlet integrator is unchanged
- All boundary conditions work the same
- The SFA output format is unchanged

## Test Plan

1. Single proton (N=64, T=100): does it survive? what's the equilibrium f_avg?
2. Compare f_avg across Θ_c values: {0.001, 0.005, 0.01, 0.05, 0.1}
3. Two protons (D=14, N=384): does the transfer improve binding vs V48?
4. Energy partition: does E_transfer + E_confine replace E_pot + E_tmass cleanly?
5. Check: can a braid form from pure background (no template seed)?
   This tests the pilot light epsilon.
