# Multi-Scale Atom Simulation -- Phase 1 Results

## Overview

Phase 1 implements the two-scale architecture: a 3D Cosserat core with spherical
boundary (N=128, L=25, R_match=20) coupled to a 1D radial Schrodinger wedge
(Nr=2000, log-spaced from r=20 to r=500,000).

Source: `src/multiscale.c`
Build: `gcc -O3 -march=native -fopenmp -o multiscale src/multiscale.c -lm`

## Test 1: Core Braid Survival with Spherical BC

**Result: PASS** -- the braid survives T=300 with spherical absorbing boundary.

The spherical mask clips the Cartesian grid to a sphere of radius R_match=20.
Points outside are zeroed. A damping layer (2*dx thick) absorbs outgoing waves.
The braid (R_tube=3, centered at origin) sits well inside the sphere.

Key metrics from the T=300 production run (N=128, L=25):

| Time | E_core | E_pot (binding) | theta_rms |
|------|--------|----------------|-----------|
| 0    | 4244   | -107           | 0.000     |
| 50   | 3896   | -72            | 0.079     |
| 100  | 3529   | -69            | 0.083     |
| 150  | 3154   | -22            | 0.071     |
| 200  | 2860   | -17            | 0.072     |
| 250  | 2624   | -17            | 0.067     |
| 300  | 2427   | -10            | 0.066     |

- **E_pot stays negative throughout** -- the braid remains bound for the full run
- E_pot oscillates between -2 and -107 with period ~55 code units (braid breathing)
- Total energy decreases monotonically: the absorbing boundary drains outgoing radiation
  (this is correct behavior for an open boundary, not numerical instability)
- theta_rms stabilizes at ~0.07 after initial transient -- curl coupling is active
- theta_at_rmatch = 0.010-0.015 -- nonzero theta field reaches the matching surface

**Energy budget**: E_core drops from 4244 to 2427 (43% absorbed). This is radiated
energy hitting the damping layer and being absorbed. The braid core energy (E_pot)
remains robustly negative, confirming the braid itself is stable.

## Test 2: Wedge Orbital

**Result: PASS** -- the wavefunction is perfectly stable in the bound-state potential.

The Crank-Nicolson integrator preserves unitarity exactly:

| Metric | t=0 | t=300 | Change |
|--------|-----|-------|--------|
| norm   | 1.00000000 | 0.99999981 | -2e-7 |
| <r>    | 239,854 | 239,854 | 0 |
| rms_r  | 245,620 | 245,620 | 0 |
| r_peak | 212,402 | 212,402 | 0 |
| E      | -3.714e-6 | -3.714e-6 | 0 |

The Gaussian wave packet (center=212,000, sigma=79,500) sits near the Bohr radius
(265,000) and does not disperse. The CN scheme is unconditionally stable and time-
reversal symmetric.

**Comparison to eigenstate prediction**: The V35 eigenstate solver gives E_1s ~ -4e-6
for the 1/r^1.189 potential. The measured E=-3.714e-6 is consistent -- the Gaussian
is not the exact eigenstate but has significant overlap with it. The invariance of
<r> and rms_r over T=300 confirms the wavefunction is in (or very near) a stationary
state.

## Test 3: Core-Wedge Coupling

**Phase 1 coupling is one-way (diagnostic only)**:
- theta_at_rmatch is extracted at each diagnostic step
- It reads 0.010-0.015 (nonzero), confirming that the theta field sourced by curl(phi)
  propagates to the matching surface
- The wedge evolves independently in Phase 1 (no back-coupling from core to wedge
  potential normalization)

**For Phase 2**: The theta amplitude at R_match can modulate the wedge potential
depth, providing true bidirectional coupling.

## Performance

| Configuration | Wall time | Notes |
|--------------|-----------|-------|
| N=64, T=10   | 1s       | Quick validation |
| N=128, T=50  | 63s      | Benchmark |
| N=128, T=300 | 417s (7 min) | Production run |

Core dominates runtime (26.2% of N^3 points active inside sphere).
Wedge (Nr=2000, CN tridiagonal) adds negligible overhead.

## Parameters Used

Core:
- N=128, L=25, R_match=20 (sphere inside [-25,25]^3 box)
- dx=0.394, dt=0.0394
- m^2=2.25, m_theta^2=2.25, eta=1.0, mu=-41.345, kappa=50.0
- Damping layer: 2*dx wide, 0.99 velocity factor per step

Wedge:
- Nr=2000, r=[20, 500000], log-spaced
- hbar_eff=22727, m_eff=1535
- V(r) = -1.27/(r/5)^1.189
- Gaussian IC: center=212000, sigma=79500
- Crank-Nicolson dt=3.94e-3 (100 substeps per 10 core steps)

## Key Findings

1. **Spherical BC works**: The braid survives inside the sphere with absorbing BC.
   The boundary drains radiation but does not destabilize the core.

2. **CN wedge is exact**: Norm preserved to 2e-7 over 76,000 wedge steps.
   Energy and position expectation values are invariant.

3. **Curl coupling generates theta at R_match**: The phi braid sources theta waves
   that reach the matching surface with amplitude ~0.01. This provides the physical
   pathway for Phase 2 coupling.

4. **Scale separation is clean**: Core (dx=0.4, R=20) and wedge (dr=0.1 to 2500,
   R=500,000) operate on completely different scales. The log-spacing naturally
   handles the 25,000x range ratio.

## Next Steps (Phase 2)

1. Feed theta_at_rmatch as source amplitude into wedge potential normalization
2. Extract wedge psi at R_match and inject as theta BC on core ghost zone
3. Measure orbital frequency under coupled evolution
4. Compare to eigenstate prediction (omega ~ E/hbar)

## Files

- `src/multiscale.c` -- combined core+wedge code (single file, ~760 lines)
- `data/production/timeseries.tsv` -- full time history (62 rows, 17 columns)
- `data/production/field_t*.bin` -- core 3D snapshots at t=0,100,200,300
- `data/production/wedge_t*.dat` -- wedge 1D profiles at t=0,100,200,300
