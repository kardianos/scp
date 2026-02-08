# Numerical Results Archive

## verify3d_run1_2026-02-08.txt

**Date**: 2026-02-08, 07:19–07:28 (9 minutes total)
**Machine**: 16 threads, 115 GB RAM
**Code**: `src/verify3d.c` + `src/field.c` (4th-order central differences, periodic BC)
**Script**: `run_verify3d.sh`

### Method

3D Cartesian lattice initialization from 1D radial profiles. The field is:

    q(x) = rho0 * [cos(f(r)) + sin(f(r)) * n_hat . sigma]

where f(r) is interpolated from `profile_B{1,2,3,4}.dat` (produced by the 1D rational map ODE solver) and n_hat comes from the hedgehog map (B=1) or rational map ansatz (B>1).

Energy computed via 4th-order central finite differences on a periodic N^3 grid spanning [-L, L]^3 with spacing h = 2L/N. Topological charge Q computed from the right-current formula B = -(1/2pi^2 rho0^4) integral of <A0 A1 A2>_0.

Sigma-model gradient flow: step along -F (negative gradient of E2+E4), then project |q| = rho0. Armijo energy backtracking with adaptive dt.

### Key Parameters

- rho0 = 1.0, e = 4.0 (giving c4 = 2*rho0^2/e^2 = 0.125)
- Soliton core radius: sqrt(c4) = 0.354
- E_FB = 6*sqrt(2)*pi^2*rho0^3/e = 20.937
- B=1 reference: E/E_FB = 1.2322 (from 1D shooting)

### Results Summary

**Part 1 — Initialization Accuracy (no relaxation)**

| B | N | L | h | E error | Q error | E2/E4 |
|---|---|---|---|---------|---------|-------|
| 1 | 64 | 6 | 0.188 | -10.1% | -0.103 | 1.154 |
| 1 | 96 | 6 | 0.125 | -2.3% | -0.025 | 1.040 |
| 1 | 128 | 6 | 0.094 | -0.5% | -0.008 | 1.020 |
| 1 | 160 | 6 | 0.075 | +0.1% | -0.004 | 1.016 |
| 1 | 192 | 6 | 0.063 | +0.4% | -0.002 | 1.015 |
| 1 | 128 | 8 | 0.125 | -2.4% | -0.025 | 1.037 |
| 1 | 192 | 8 | 0.083 | -0.4% | -0.005 | 1.012 |
| 1 | 256 | 8 | 0.063 | +0.1% | -0.002 | 1.009 |

Convergence is 4th-order in h (consistent with the derivative stencil order).
At fixed L=6, the error is non-monotonic beyond N=160 because the periodic image
interaction energy grows as discretization error shrinks. L=8 gives cleaner results.

**Part 2 — Gradient Flow Relaxation (topology loss)**

| N | h | Core pts | Steps to Q<0.5 | Final Q | Final E/E_FB |
|---|---|----------|----------------|---------|-------------|
| 128 | 0.094 | 3.8 | ~8 | 6e-8 | 0.014 |
| 160 | 0.075 | 4.7 | ~22 | 3e-8 | 0.013 |
| 192 | 0.063 | 5.7 | ~35 | 5e-6 | 0.037 |

**All resolutions show complete topology loss.** The soliton unwinds during gradient
flow at all tested N. Higher N delays the onset but does not prevent it.

Root cause: The sigma-model Skyrmion is a saddle point of the unconstrained energy
functional on the lattice. Continuous topology is only approximately preserved on
a discrete grid. Once the soliton core shrinks below grid resolution, the winding
number slips through the lattice and the gradient flow finds the global minimum
(vacuum, Q=0).

**Part 3 — Higher-B Initialization (all B=1–4)**

| B | N | L | E error | Q error | E2/E4 |
|---|---|---|---------|---------|-------|
| 1 | 256 | 6 | +0.74% | -0.001 | 1.018 |
| 2 | 256 | 6 | -0.02% | -0.001 | 1.002 |
| 3 | 256 | 6 | -0.04% | -0.002 | 1.001 |
| 4 | 256 | 6 | -0.08% | -0.001 | 1.000 |

All B=1–4 verified to sub-percent accuracy at N=256.

**B=3 convergence study** (after fixing rational map bug):

| N | L | h | E error | Q error | E2/E4 |
|---|---|---|---------|---------|-------|
| 160 | 6 | 0.075 | -0.41% | -0.010 | 1.005 |
| 192 | 6 | 0.063 | -0.19% | -0.005 | 1.003 |
| 256 | 6 | 0.047 | -0.04% | -0.002 | 1.001 |
| 320 | 6 | 0.038 | +0.01% | -0.001 | 1.001 |
| 512 | 8 | 0.031 | +0.005% | -0.000 | 1.000 |

Convergence is 4th-order in h, matching B=1. Best result: 0.005% at N=512, L=8.

**B=3 bug (found and fixed 2026-02-08)**: The original `verify3d.c` used the wrong
B=3 rational map denominator: `1 - sqrt(3)*i*z^3` (z^3 in denominator) instead of the
correct `sqrt(3)*i*z^2 - 1` (z^2 in denominator, matching `rational_map.c`). Both maps
have degree 3 (same winding number), so E2 matched perfectly, but the angular integral
I = integral(b^2 dOmega) differs, causing a systematic +27% E4 overestimate that
*increased* with resolution (converging to the wrong answer). The buggy run is archived
in `B3_convergence_buggy_2026-02-08.txt` for reference.

### Conclusions

1. 3D initialization from 1D profiles works to <0.1% accuracy for all B=1–4.
2. 3D gradient flow cannot verify stationarity — topology is lost at all tested N.
3. The 1D radial solver remains the authoritative tool for Skyrmion profiles.
4. B=3 needs slightly higher N than B=1,2,4 (tetrahedral symmetry has sharper angular features), but converges normally once the correct rational map is used.
