# V27-M5: Topological Confinement via Chern-Simons Term

## Configuration
- Grid: N=128, L=20, periodic z, absorbing xy boundaries
- Fields: m=0, mu=-50, kappa=50 (V27-M4 optimal m=0 propagating braid)
- A0=0.8, R_tube=3.0, k_wave=2*pi/(2L)

## Test M5a: Topological Charges of the Propagating Braid

Evolved m=0 propagating braid for t=500. Computed at t=0,100,...,500:

| t   | Q_jac        | Q_CS_xy      | phase_wind | peak_P | fc    | E_total |
|-----|-------------|-------------|-----------|--------|-------|---------|
| 0   | -9.04e-33   | -2.28e-33   | -1.000    | 0.127  | 0.983 | 25.8    |
| 100 | -4.01e-16   | 2.27e-16    | -1.000    | 1.515  | 0.364 | -8688   |
| 200 | 6.79e-18    | 7.24e-17    | -1.000    | 1.796  | 0.365 | -8830   |
| 300 | 8.92e-17    | -2.21e-17   | -1.000    | 0.911  | 0.367 | -8901   |
| 400 | -2.02e-17   | 8.54e-17    | -1.000    | 1.102  | 0.363 | -8950   |
| 500 | 1.87e-18    | -3.70e-17   | -1.000    | 2.506  | 0.365 | -8970   |

### Key findings:

1. **Q_jac (Jacobian determinant integral) = 0**: The integral of det(d_i phi_a)
   is identically zero to machine precision (~10^-17) at all times. This is because
   the helical braid ansatz maps R^3 -> R^3 with zero topological degree (the fields
   are localized Gaussians times cosines, not a surjective map).

2. **Q_CS_xy = 0**: The Chern-Simons density eps_{abc} phi_a (dx phi_b)(dy phi_c)
   also integrates to zero. Same reason: the Gaussian envelope kills any net winding
   in the xy-plane.

3. **Phase winding is EXACTLY conserved**: wind = -1.000 at every single diagnostic
   step from t=0 to t=500, including all 50 fine-grained measurements. The helical
   phase structure kz + 2*pi*a/3 is perfectly maintained by the dynamics.

4. **|P| (triple product) is robustly maintained**: oscillates between 0.9 and 3.5
   throughout the evolution, never dropping below 0.1. The braid binding is stable.

5. **Helicity integrals H_12, H_23, H_31**: All ~10^-18 (zero). No net helicity
   in the field configuration.

## Test M5c: Static Braid Comparison

Same parameters but with zero initial velocity (no propagation).

| t   | Q_jac        | Q_CS_xy      | phase_wind | peak_P | fc    | E_total |
|-----|-------------|-------------|-----------|--------|-------|---------|
| 0   | -9.07e-33   | -2.25e-33   | -1.000    | 0.127  | 0.966 | 12.3    |
| 100 | 1.69e-17    | -5.48e-17   | -1.000    | 1.378  | 0.365 | -8697   |
| 200 | 2.21e-17    | 3.08e-17    | -1.000    | 1.683  | 0.367 | -8846   |
| 300 | -2.05e-18   | 4.43e-17    | -1.000    | 0.975  | 0.367 | -8905   |
| 400 | -2.09e-17   | 1.94e-17    | -1.000    | 1.077  | 0.364 | -8958   |
| 500 | 3.70e-17    | 2.00e-17    | -1.000    | 2.298  | 0.366 | -8970   |

### Key findings:

1. **Static braid ALSO maintains topology**: wind = -1.000, |P| > 0.9 at all times.
   This is DIFFERENT from V26 where the static braid unwound (|P| decayed 0.13 -> 0.001).

2. **The difference from V26**: The mu=-50, kappa=50 strong binding (vs V26's mu=-20,
   kappa=20) is what prevents unwinding. The stronger triple-product potential creates
   a deeper well that the braid cannot escape from.

3. **Propagation is NOT required for topology preservation**: Both static and propagating
   braids maintain identical topological measures. The binding potential alone is sufficient.

4. **Energetics are nearly identical**: Both converge to E ~ -8970 by t=500.

## Test M5b: Chern-Simons Potential Scan

Added global CS potential: V_CS = lambda_CS * (Q_CS_xy - Q0)^2 where Q0 is initial value.

| lambda_CS | Q_jac(0)   | Q_jac(500) | Q_CS_xy(0) | Q_CS_xy(500) | wind(0/500) | peak_P(500) |
|-----------|-----------|-----------|-----------|-------------|------------|------------|
| 0.01      | -9.06e-33 | 2.23e-17  | -2.27e-33 | 1.57e-17    | -1/-1      | 2.506       |
| 0.10      | -9.06e-33 | -1.14e-17 | -2.25e-33 | 1.57e-17    | -1/-1      | 2.506       |
| 1.00      | -9.04e-33 | -2.96e-17 | -2.27e-33 | 6.74e-17    | -1/-1      | 2.506       |

### Key findings:

1. **CS penalty has NO EFFECT**: All three lambda_CS values produce identical evolution
   to the lambda_CS=0 case. The reason: the topological charges Q_jac and Q_CS_xy
   are identically zero (within machine precision), so Q_CS - Q0 ~ 10^-17 and the
   penalty force is negligible.

2. **The Jacobian density is a total divergence**: The analysis in the code proves that
   d_m[cofactor_{a,m}] = 0 identically (second derivatives cancel by symmetry of
   mixed partials + antisymmetry of epsilon). So the Jacobian integral is EXACTLY
   conserved by ANY smooth dynamics, making the penalty redundant.

3. **The xy-plane CS density also integrates to zero**: For the Gaussian-envelope helical
   ansatz, the net CS charge vanishes due to the radial symmetry of the envelope.

## Summary and Implications

### The topology of the m=0 braid is AUTOMATICALLY conserved

The phase winding number w = -1 is maintained exactly throughout all 500 time units
for both propagating and static braids at mu=-50, kappa=50. No Chern-Simons term
is needed.

### Why the topology survives:

1. **Phase winding is protected by periodicity**: With periodic z-boundary conditions,
   the winding number (integral of d(theta)/dz around the z-cycle) can only change
   by integers. Since the dynamics are smooth and continuous, the winding number
   is trapped at -1.

2. **Strong binding prevents amplitude collapse**: At mu=-50, kappa=50, the triple-product
   potential well is deep enough to prevent the field amplitudes from decaying to zero
   (which would make the phase undefined and allow topology change).

3. **The Jacobian integral is a total divergence**: Q_jac = integral of det(J) is
   exactly conserved by any smooth dynamics. This is a mathematical identity, not
   a dynamical statement.

### What this means for confinement:

The braid is confined by TWO mechanisms working together:
- **Topological**: periodic BC + smooth dynamics = exact winding conservation
- **Energetic**: strong triple-product binding prevents amplitude death

The Chern-Simons term is unnecessary because the topology is already exactly conserved.
The relevant question for confinement shifts from "how to conserve linking number" to
"how to prevent amplitude collapse" -- which is solved by the strong binding potential.

### Comparison with V26:

| Property         | V26 (mu=-20,k=20) | V27-M5 (mu=-50,k=50) |
|-----------------|-------|---------|
| Static |P| at t=500 | 0.001 (dead) | 2.30 (alive) |
| Phase winding   | lost  | -1.000 (conserved) |
| Mechanism       | weak binding | strong binding |

The V26 "unwinding" was actually amplitude death (|P| -> 0), not topological change.
With sufficient binding strength, the topology is automatically protected.

## Files

- `data/m5a_topo_propagating.tsv` — full topology diagnostics, propagating braid
- `data/m5a_topo_fine_prop.tsv` — fine-grained topology (every ~10 time units)
- `data/m5c_topo_static.tsv` — full topology diagnostics, static braid
- `data/m5c_topo_fine_static.tsv` — fine-grained topology, static
- `data/m5b_topo_lcsX.XXXX.tsv` — topology with CS potential at each lambda_CS
- `data/m5b_summary.tsv` — CS scan summary
