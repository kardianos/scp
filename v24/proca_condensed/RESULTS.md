# V24-P2: Condensed Phase Goldstone — RESULTS

## Summary

**NO true Goldstone mode exists** when the triple-product oscillon coupling is present.
The mode that would be Goldstone is a **pseudo-Goldstone** with mass that grows with
the condensate VEV. This is because the triple product P = phi_1 phi_2 phi_3 explicitly
breaks the SO(2) symmetry in the antisymmetric subspace.

The oscillon does NOT survive well in the condensed vacuum — the symmetric mode
frequency sits at or above the effective mass gap, so the oscillon slowly radiates.

## Setup

Three massive scalars with:
- Pairwise coupling: lambda(phi_1 phi_2 + phi_2 phi_3 + phi_3 phi_1)
- Triple-product (oscillon binding): (mu/2) P^2 / (1 + kappa P^2)
- O(3)-symmetric quartic (condensate stabilization): (g4/4)(sum phi_a^2)^2
- Parameters: mu=-20, kappa=20, m=1.0, g4=1.0

When lambda > m^2 = 1.0: antisymmetric modes are tachyonic in phi=0 vacuum.
The quartic g4 stabilizes a condensate at finite VEV.

**Key modification**: The original proposal used only the triple-product for
stabilization. This fails because mu < 0 makes the triple product attractive
(unbounded below). An O(3)-symmetric quartic (g4/4)(sum phi^2)^2 was added
to stabilize the condensate.

## Phase 1: Condensed Vacuum Structure

### A. Without triple product (mu=0) — Exact Goldstone

The O(3)-symmetric quartic + pairwise coupling has exact SO(2) symmetry in the
antisymmetric subspace. VEV = sqrt(2(lambda-m^2)/(3g4)).

| lambda | |VEV| | V_min | m^2_Gold | m^2_rad | m^2_sym |
|--------|--------|---------|----------|---------|---------|
| 1.01 | 0.100 | -2.5e-5 | 0.000000 | 0.020 | 3.030 |
| 1.05 | 0.224 | -6.3e-4 | 0.000000 | 0.100 | 3.150 |
| 1.10 | 0.316 | -2.5e-3 | 0.000000 | 0.200 | 3.300 |
| 1.20 | 0.447 | -1.0e-2 | 0.000000 | 0.400 | 3.600 |
| 1.50 | 0.707 | -6.3e-2 | 0.000000 | 1.000 | 4.500 |

All eigenvalues match analytical predictions exactly:
- Goldstone mass^2 = 0 (exact to ~1e-8)
- Radial mass^2 = 2(lambda - m^2) (Higgs mode)
- Symmetric mass^2 = m^2 + 2*lambda + 2*g4*|VEV|^2

### B. With triple product (mu=-20, kappa=20) — Pseudo-Goldstone

The triple product P = phi_1 phi_2 phi_3 explicitly breaks SO(2) -> discrete.
The "Goldstone" becomes a pseudo-Goldstone with mass growing as VEV^3.

| lambda | |VEV| | m^2_pseudo | m_pseudo | m^2_rad | m^2_sym |
|--------|--------|------------|----------|---------|---------|
| 1.01 | 0.101 | 0.000341 | 0.018 | 0.020 | 3.030 |
| 1.05 | 0.231 | 0.009428 | 0.097 | 0.094 | 3.148 |
| 1.10 | 0.339 | 0.044161 | 0.210 | 0.170 | 3.294 |
| 1.20 | 0.564 | 0.149722 | 0.387 | 0.343 | 3.611 |
| 1.50 | 1.059 | 1.680546 | 1.296 | 2.793 | 5.706 |

**Critical observation**: At lambda=1.5, the pseudo-Goldstone mass (1.30) exceeds the
bare mass (1.0). The triple product makes ALL modes heavy — no light mediator.

## Phase 2: Fluctuation Spectrum

### mu=0, lambda=1.01 (exact Goldstone)
- Goldstone wavefront propagates at v_group = 1.001 (= c)
- Linear dispersion omega = c*k confirmed (massless)
- Wavefront arrival at x=80: t=72.5, v=1.00

### mu=-20, lambda=1.01 (lightest pseudo-Goldstone)
- Pseudo-Goldstone mass = 0.018
- Wavefront still propagates at v ~= 1.0 (mass so small it's nearly massless)
- Peak frequency omega = 0.035 (close to but above mass)

### mu=-20, lambda=1.1
- Pseudo-Goldstone mass = 0.21, clearly massive
- Peak frequency omega = 0.21 at all probe points (standing oscillation, not traveling)
- Propagation speed 0.99 (slightly below c due to mass)
- This is Yukawa, NOT 1/r

## Phase 3: Oscillon in Condensed Vacuum

Oscillons initialized as symmetric Gaussian (A=0.8, sigma=3) on top of condensate.

| lambda | fc(t=0) | fc(t=2500) | fc(t=5000) | omega_osc | m_eff_sym | Stable? |
|--------|---------|------------|------------|-----------|-----------|---------|
| 1.01 | 1.000 | ~0.19 | ~0.19 | 1.744 | 1.741 | MARGINAL |
| 1.10 | 1.000 | ~0.27 | ~0.09 | 1.824 | 1.815 | NO |
| 1.50 | 1.000 | ~0.15 | ~0.71* | 2.128 | 2.121 | MARGINAL |

*lambda=1.5 with mu=0 shows better confinement; with mu=-20, fc drops faster.

The oscillon frequency sits at or just above the effective symmetric mass gap.
No clear below-gap oscillon regime exists in the condensed vacuum.

## Key Physics

### Why no true Goldstone?
The continuous SO(2) symmetry in the antisymmetric plane (rotations preserving
phi_1+phi_2+phi_3=0 and |phi|^2=const) is explicitly broken by the triple product
P = phi_1 phi_2 phi_3. Under SO(2) rotation by angle theta:
- phi_1' = phi_1 (unchanged, it's the radial direction)
- (phi_2', phi_3') rotate in the Goldstone plane

But P = phi_1 phi_2 phi_3 = phi_1 * (phi_2 cos(theta) - phi_3 sin(theta)) *
(phi_2 sin(theta) + phi_3 cos(theta)), which is NOT invariant under SO(2).

### When is the pseudo-Goldstone approximately massless?
Near threshold (lambda -> m^2+), the VEV -> 0, and the triple product contribution
scales as VEV^6, so the pseudo-Goldstone mass^2 ~ VEV^4 -> 0. At lambda=1.01,
m_pseudo = 0.018 << 1. This gives a long-range Yukawa tail with range ~ 1/0.018 = 56
code lengths. But it's never truly massless.

### Can the oscillon survive?
The symmetric mode mass^2_S = m^2 + 2*lambda + 2*g4*|VEV|^2 grows with lambda.
The oscillon frequency must be below sqrt(m^2_S) to be trapped. With the quartic
contribution, the gap is large enough that the Gaussian initial condition places
the oscillation frequency right at the gap edge — not safely below it.

## Conclusion

1. **Exact Goldstone (mu=0)**: EXISTS when the only quartic is O(3)-symmetric.
   Confirmed by eigenvalue analysis (m^2 = 0 to 1e-8) and wavefront propagation
   at v = c.

2. **With triple product (mu=-20)**: The Goldstone becomes a PSEUDO-GOLDSTONE
   with mass proportional to VEV^2. Only approximately massless near threshold.
   This is Yukawa, not 1/r.

3. **Oscillon survival**: MARGINAL. The symmetric mode frequency sits near the
   mass gap edge. No robust oscillon regime found in the condensed vacuum.

4. **For 1/r mediation**: Would need to remove or modify the triple-product coupling
   to restore the exact SO(2) symmetry. But without the triple product, there is no
   oscillon binding mechanism.

**Bottom line**: The mechanisms for oscillon binding (triple product) and Goldstone
mode (continuous symmetry breaking) are incompatible. The triple product that binds
the oscillon also explicitly breaks the symmetry that would give a massless Goldstone.
