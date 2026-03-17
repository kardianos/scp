# V26-DynB: Rotating Braid --- Results

## Parameters
- mu=-20, kappa=20, m=1.0, A0=0.8, R_tube=3.0
- k = 2*pi/(2L) (one full twist across domain)
- N=128, L=20, dx=0.315, dt=0.063, tfinal=200
- Periodic BC in z, absorbing in x,y
- Omega scan: {0.05, 0.1, 0.2, 0.5}

## Initialization
```
phi_a = A(r_perp) * cos(theta + k*z + 2*pi*a/3)
vel_a = Omega * A(r_perp) * sin(theta + k*z + 2*pi*a/3)
```
where theta = atan2(y,x) gives azimuthal winding.

## Summary Table

| Omega | fc_i   | fc_f   | |P|_i   | |P|_f     | Lz_i     | Lz_f    | l2_frac | survived |
|-------|--------|--------|---------|-----------|----------|---------|---------|----------|
| 0.05  | 0.9990 | 0.7127 | 0.1269  | 2.97e-05  | -54.14   | -6.93   | 0.070   | YES      |
| 0.10  | 0.9990 | 0.7124 | 0.1269  | 3.00e-05  | -108.28  | -13.85  | 0.073   | YES      |
| 0.20  | 0.9990 | 0.7113 | 0.1269  | 3.14e-05  | -216.57  | -27.70  | 0.081   | YES      |
| 0.50  | 0.9991 | 0.7048 | 0.1269  | 4.37e-05  | -541.41  | -69.33  | 0.123   | YES      |

## Key Findings

### 1. All configurations survive (fc > 0.70)
Every Omega value produces a surviving localized structure at t=200.
Core fraction fc drops from ~1.0 to ~0.71, nearly independent of Omega.
The transient dynamics are also Omega-independent: fc dips to ~0.18 at t=40,
recovers to ~0.69 at t=100, stabilizes near 0.71.

### 2. Triple product collapses universally
|P|_final ~ 3e-05 for all Omega (down from 0.127). The triple-product
coupling is NOT maintaining the 120-degree phase structure. The surviving
energy is in the mass term (individual phi_a oscillations), not the braid.

### 3. Angular momentum Lz scales linearly with Omega
- Lz_init = -1082.8 * Omega (exactly proportional)
- Lz_final ~ -138.5 * Omega (also proportional)
- Retention ratio: Lz_f/Lz_i ~ 0.128 for all Omega

Angular momentum is NOT conserved --- it is lost to the absorbing boundary.
About 87% of initial Lz is radiated away by t=200. The retention ratio is
Omega-independent, meaning angular momentum dynamics decouple from Omega.

### 4. l=2 multipole increases with Omega
| Omega | l=0 frac | l=2 frac |
|-------|----------|----------|
| 0.05  | 0.929    | 0.070    |
| 0.10  | 0.927    | 0.073    |
| 0.20  | 0.918    | 0.081    |
| 0.50  | 0.877    | 0.123    |

Higher Omega shifts strain from l=0 (monopole) toward l=2 (quadrupole).
At Omega=0.5, l=2 content is 1.75x the Omega=0.05 value. This is the
expected signature of rotation: centrifugal deformation produces quadrupolar
strain on the measurement shell.

Initial state has l=2 fraction = 0.44 (strong quadrupole from azimuthal winding).
This decays to 0.07-0.12 as the braid disperses, retaining only the rotation signal.

### 5. DFT peak omega = 0.07 for all runs
The center density oscillates at omega_peak = 0.07 regardless of Omega.
This is NOT the rotation frequency --- it is likely the mass oscillation
frequency (m=1.0 gives omega ~ 1.0 for free field, but nonlinear effects
and the confined geometry reduce the effective frequency).

### 6. Energy: Omega determines kinetic contribution only
| Omega | Et_init | Et_final | Ek_init | Ek_final |
|-------|---------|----------|---------|----------|
| 0.05  | 913.9   | 84.9     | 1.4     | 61.9     |
| 0.10  | 918.0   | 85.3     | 5.5     | 62.0     |
| 0.20  | 934.4   | 87.3     | 21.9    | 62.4     |
| 0.50  | 1049.3  | 101.1    | 136.8   | 65.4     |

Initial Ek scales as Omega^2 (as expected). Final Ek is nearly Omega-independent
(~62-65), meaning the rotation energy is radiated away. The surviving structure
has roughly the same kinetic content regardless of initial Omega.

## Conclusion

Rotation does NOT stabilize the braid. The fc, |P|, and energy evolution are
nearly identical across all Omega values. The 120-degree phase structure (measured
by |P|) collapses by a factor of 4000x regardless of rotation. What survives is
a localized energy concentration (fc~0.71) that is NOT a rotating braid but rather
individual field oscillations trapped by the mass term.

The one clear Omega-dependent signature is the l=2 quadrupole content, which
increases monotonically with Omega. This confirms the rotation IS generating
quadrupolar strain, but without topological protection the braid unwinds anyway.

Gyroscopic stabilization requires a conserved angular momentum, but the absorbing
boundary radiates 87% of Lz by t=200. Even so, the 13% retained is proportional
to Omega, indicating angular momentum IS carried by the structure --- just not enough
to prevent unwinding.
