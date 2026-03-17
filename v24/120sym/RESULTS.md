# V24-SYM Results: Symmetric Penalty (phi1+phi2+phi3)^2

## Parameters

mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0, Nx=4000, xmax=100, tfinal=10000

## Summary Table

| lambda_Q | Test | omega_peak | cos(1,2) | cos(2,3) | cos(3,1) | fc(end) | E_total | sum_suppressed | State |
|----------|------|-----------|----------|----------|----------|---------|---------|----------------|-------|
| 0.0  | 0deg  | 0.876 | +1.000 | +1.000 | +1.000 | 1.000 | 1.264 | NO  | Stable 0deg oscillon |
| 0.0  | 120deg| 0.876 | -1.000 | +1.000 | -1.000 | 0.994 | 1.270 | NO  | 180deg (1 vs 2,3) |
| 0.1  | 0deg  | 1.086 | +1.000 | +1.000 | +1.000 | 1.000 | 5.680 | NO  | Stable 0deg, omega>m |
| 0.1  | 120deg| 0.900 | -1.000 | +1.000 | -1.000 | 0.996 | 1.392 | weak | 180deg oscillon |
| 0.5  | 0deg  | 1.896 | +1.000 | +1.000 | +1.000 | 1.000 | 18.06 | NO  | Stable 0deg, omega>m |
| 0.5  | 120deg| 1.002 | -0.981 | +0.931 | -0.981 | 0.49  | 1.323 | YES (9e-3) | DISPERSING |
| 1.0  | 0deg  | 2.562 | +1.000 | +1.000 | +1.000 | 1.000 | 33.78 | NO  | Stable 0deg, omega>m |
| 1.0  | 120deg| 1.002 | -0.670 | -0.123 | -0.580 | 0.22  | 0.368 | YES (7e-2) | DISPERSED (scrambled) |
| 2.0  | 0deg  | 2.982 | +1.000 | +1.000 | +1.000 | 1.000 | 64.33 | NO  | Stable 0deg, omega>m |
| 2.0  | 120deg| 0.936 | -0.993 | +1.000 | -0.993 | 0.991 | 1.772 | YES (1e-3) | Stable 180deg oscillon |
| 5.0  | 0deg  | 2.970 | +1.000 | +1.000 | +1.000 | 1.000 | 155.8 | NO  | Stable 0deg, omega>m |
| 5.0  | 120deg| 0.918 | -1.000 | +1.000 | -1.000 | 0.994 | 1.738 | YES (4e-4) | Stable 180deg oscillon |

## Key Finding: NO Spontaneous 0deg -> 120deg Transition

The 0-degree state **never** breaks symmetry. At all lambda_Q values tested (0 to 5),
the three fields remain perfectly locked in phase (cos = +1.000). The penalty
raises the effective mass of the symmetric mode to m_S^2 = m^2 + 6*lambda_Q,
which increases the oscillation frequency above the mass gap (omega > m for
lambda_Q >= 0.1), but the 0-degree configuration is a local minimum of V_pen
(S=0 at the zero crossings, maximum at the peaks). The fields oscillate together
through zero, paying the penalty cost periodically but never desynchronizing.

**Why it fails**: The symmetric penalty is quadratic in S = phi1+phi2+phi3.
The 0-degree state already has a saddle structure in S(t): S oscillates between
+3A and -3A. But at each instant, the force -2*lambda_Q*S acts identically on
all three fields (it's the SAME force on each), so it cannot break the symmetry.
A perturbation delta_a with sum(delta_a)=0 feels NO penalty force at all.
The instability would require noise to seed it, but velocity-Verlet is symplectic
and introduces no noise.

## Phase Structure: 180deg Not 120deg

The "120-degree" initialization (phases 0, 2pi/3, 4pi/3) does NOT maintain
120-degree phase separation. Instead it evolves into a **180-degree** (anti-phase)
configuration: phi1 oscillates opposite to (phi2, phi3), which are locked together.
This is the Z_2 anti-symmetric mode, not the Z_3 mode.

The cross-correlation confirms: cos(phi1,phi2) = -1.0, cos(phi2,phi3) = +1.0.
True 120-degree would give cos = -0.5 for all three pairs.

This happens because the 120-degree initial condition has two fields starting
negative (phi2, phi3 at phases 120, 240 deg both have cos < 0 at t=0), and
the nonlinear dynamics locks them together.

## Sum Suppression Works for Anti-Phase

At lambda_Q >= 2.0, the anti-phase (180deg) oscillon achieves excellent sum
suppression: power_sum/power_phi1 < 10^-3. The sum S = phi1 + phi2 + phi3
nearly cancels because phi1 ~ -2*phi2 ~ -2*phi3 with phi1 having twice the
amplitude. But this is 180deg cancellation, not 120deg.

## Stability Window

- lambda_Q = 0.5: 120deg init DISPERSES (fc drops to 0.05, amplitudes < 0.1)
- lambda_Q = 1.0: 120deg init FULLY DISPERSED (phases scrambled, E drops to 0.37)
- lambda_Q = 2.0, 5.0: 120deg init survives as 180deg oscillon

The dispersal at intermediate lambda_Q is because the penalty cost exceeds the
triple-product binding. At higher lambda_Q, the 180deg mode (small S) pays
minimal penalty and remains bound.

## Conclusion

The symmetric penalty V_pen = lambda_Q*(phi1+phi2+phi3)^2 does NOT achieve the
desired 0deg -> 120deg transition. The mechanism fails because:

1. The penalty force is identical on all three fields (cannot break degeneracy)
2. The natural anti-symmetric mode is 180deg (Z_2), not 120deg (Z_3)
3. True 120deg requires a Z_3-symmetric coupling (e.g., cos(phi1-phi2-...))

To get true 120-degree phase locking, one would need either:
- Explicit Z_3 phase coupling (circular permutation symmetry)
- Noise/thermal perturbation to break the 0deg saddle
- A different potential that penalizes pairwise alignment
