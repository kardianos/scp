# T7: Two-Braid Interaction — Results

## Setup
- Grid: N=192, L=40, dx=0.419, dt=0.0838
- Two braids at x=+/-15, initial separation D=30
- BIMODAL sweet spot params, mass=1.50 (m^2=2.25), mass2_override=-1
- Absorbing xy BC, periodic z
- T=500, diagnostics every T=25 (20 snapshots per config)
- Wall time: ~46 min total (22 min config 1, 13 min config 2, 12 min config 3)

## Configurations
| Config | Description | D(0) | D(500) | deltaD | Direction |
|--------|-------------|-------|--------|--------|-----------|
| Same twist | Both W=-1 | 30.00 | 33.16 | +3.16 | REPULSION |
| Opposite twist | W=-1, W=+1 (phi_1 flipped) | 30.00 | 32.51 | +2.51 | REPULSION |
| Control | Single braid at center | -- | -- | -- | -- |

## Key Findings

### 1. Both configurations show REPULSION
Both same-twist and opposite-twist braids repel, with similar magnitude
(deltaD ~ +2.5 to +3.2 over T=500). The separation grows monotonically
after an initial transient oscillation (T ~ 0-150).

### 2. Repulsion is weakly twist-dependent
- Same twist: deltaD = +3.16 (slightly stronger repulsion)
- Opposite twist: deltaD = +2.51

The difference (0.65 code units, ~20% of signal) suggests a small
twist-dependent component on top of a dominant twist-independent repulsion.
This is closer to "universal repulsion" than "EM-like" behavior.

### 3. Quadratic acceleration fits
- Same twist: a = 1.15e-5 (D units/T^2)
- Opposite twist: a = 7.77e-6 (D units/T^2)

These are tiny accelerations — the dominant effect is radiation-driven drift,
not a clean force law.

### 4. Energy loss is identical across configs
All three configurations lose ~38.5% of energy over T=500:
- Same twist: dE/E = -38.5%
- Opposite twist: dE/E = -38.4%
- Control (single braid): dE/E = -38.5%

The energy loss per braid is the same whether alone or paired, indicating
the braids are too far apart (D=30) for significant energy exchange.
The loss is from radiation into the damping layer.

### 5. Transient dynamics
Both two-braid configs show a characteristic oscillation in D:
- T=0-75: D drops from 30 to ~28.4 (initial compression from overlapping tails)
- T=75-200: D recovers to ~30 (relaxation)
- T=200-500: D grows monotonically to ~33 (sustained repulsion)

The fc (core fraction) tracks this: it dips at T~50 (radiation emission),
recovers by T~125, then remains stable at ~0.92.

### 6. Control braid behavior
The single centered braid has fc ~ 0.01-0.08 (very low), because the
com_half diagnostic splits a centered braid across the x=0 boundary.
Its D fluctuates around 5-10 (noise from asymmetric radiation).
Energy loss matches the two-braid configs per-braid (E_single ~ E_double/2).

## Interpretation

The repulsion is likely dominated by **radiation pressure** from the
oscillating field tails, not a long-range force. Evidence:
1. Both twist orientations repel (not EM-like)
2. The magnitude is small (deltaD/D ~ 10% over T=500)
3. Energy loss rate is identical with or without a neighbor
4. The acceleration is not constant — it builds after the transient

The twist-dependent component (20% difference) may reflect interference
between the phi_1 fields: same-twist braids have constructive interference
in the overlap region, producing slightly more radiation pressure.

## Data
- `data/t7_results.tsv` — full time series (config, T, D, E, fc_left, fc_right)
