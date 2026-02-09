# Hopfion Search — Numerical Skyrmion Solver

Numerical solver for soliton (Skyrmion) solutions of the CHPT Lagrangian in Cl+(3,0,1).

## Directory Structure

```
hopfion_search/
  src/        Source code (C11, headers)
  bin/        Compiled executables (built by make)
  build/      Object files (built by make)
  data/       Generated data: profiles, Maxwell fields, scatter output
  scripts/    Shell and Python scripts for batch runs and visualization
  results/    Analysis documents, archived run logs, figures
  Makefile    Build system (gcc + OpenMP)
```

## Building

```bash
make          # build all targets into bin/
make clean    # remove bin/ and build/
make scatter  # build only the scatter executable
```

Requires: `gcc` with OpenMP support, `libm`.

## Executables

| Binary | Source | Description |
|--------|--------|-------------|
| `bin/radial_solver` | `src/radial.c` | B=1 Skyrmion profile via shooting method |
| `bin/rational_map_solver` | `src/rational_map.c` | B=1-4 profiles via rational map ansatz |
| `bin/finite_lambda_solver` | `src/finite_lambda.c` | Coupled f-shooting + rho-BVP at finite lambda |
| `bin/scatter` | `src/scatter.c` | Soliton-soliton scattering (leapfrog time evolution) |
| `bin/verify3d` | `src/verify3d.c` | 3D grid initialization from 1D profiles, energy check |
| `bin/maxwell` | `src/maxwell.c` | Sourced Maxwell equations from soliton current |
| `bin/veff` | `src/veff.c` | Effective potential for degenerate modes |
| `bin/soliton_search` | `src/main.c` | 3D gradient flow (historical, topology loss on lattice) |
| `bin/verify_gradient` | `src/verify.c` | Gradient verification (finite-difference check) |

## Common Workflows

**Generate a radial profile:**
```bash
bin/radial_solver -e 1 -o data/profile_sigma_e1.dat
```

**Run a scattering simulation (repulsive B+B at v=0.5c):**
```bash
export OMP_NUM_THREADS=8
bin/scatter -v 0.5 -z 1.5 -N 192 -L 10 -e 1 -lambda 5000 \
  -dt 0.001 -T 5 -out 50 -isorot -profile data/profile_sigma_e1.dat
```

**Compute finite-lambda profile:**
```bash
bin/finite_lambda_solver -e 2 -lambda 5000 -o data/profile_finlam_e2_5000.dat
```

**Verify 3D initialization from 1D profile:**
```bash
export OMP_NUM_THREADS=8
bin/verify3d -B 1 -N 256 -L 8 -profile data/profile_B1.dat
```

## Key Parameters for Scattering

The critical parameter for lattice stability is grid points across the soliton core:
- Core radius = sqrt(2 * rho0^2 / e^2). For e=1: sqrt(2) = 1.414.
- Need >= 13 grid points across core for ~2.4 time units of stability.
- Use sigma-model profiles (generated with `radial_solver`) to eliminate breathing modes.
- No damping or settling -- conservative evolution only.

Working collision parameters: `e=1, N=192, L=10, sigma-model profile, v=0.5c, z0=1.5, lambda=5000, dt=0.001`.

## Data Files

Profile naming convention:
- `profile_B{N}.dat` — sigma-model profiles at e=4 for B=1-4
- `profile_sigma_e{N}.dat` — sigma-model profiles at specified e
- `profile_finlam_e{E}_{L}.dat` — finite-lambda profiles at specified e and lambda
- `profile_lam{X}.dat` — finite-lambda sweep profiles (historical)

Output data:
- `scatter_*.dat` — time series from scattering simulations
- `maxwell_*.dat` — electromagnetic field data from soliton currents
- `veff_*.dat` — effective potential data for degenerate modes
