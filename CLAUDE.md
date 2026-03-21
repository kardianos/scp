# SCP Project Conventions

## Data Format Policy
- **ALL simulation output MUST use SFA format** (`.sfa` files via `sfa.h`)
- **NO `.bin` file output** — remove `save_field()` functions that write raw binary
- SFA header: `/home/d/code/scp/sfa/format/sfa.h` — single-header C library
- Usage: `#define SFA_IMPLEMENTATION` in exactly ONE .c file before including
- Link with: `-lzstd`
- Column convention: phi_x/phi_y/phi_z (SFA_POSITION), theta_x/theta_y/theta_z (SFA_ANGLE)
- All analysis tools must READ SFA files using `sfa_open()` / `sfa_read_frame()`
- SFA viewer: `/home/d/code/scp/sfa/viewer/volview`

## Simulation Parameters (standard)
- Equation: 6-field Cosserat (3 phi + 3 theta)
- V(P) = (mu/2) P^2 / (1 + kappa P^2), P = phi_0 * phi_1 * phi_2
- Default: m^2=2.25, mu=-41.345, kappa=50, eta=0.5, A_bg=0.1
- Phase offsets: delta = {0, 3.0005, 4.4325} (from v28 optimization)

## Build Convention
- C with OpenMP: `gcc -O3 -march=native -fopenmp -o <binary> src/<file>.c -lzstd -lm`
- **ONE copy of sfa.h** at `/home/d/code/scp/sfa/format/sfa.h` — do NOT copy to other directories
- Include via relative path: `#include "../../sfa/format/sfa.h"` (adjust depth as needed)
- No external deps beyond zstd, OpenMP, math

## Diagnostics Requirements
- Every simulation MUST include fragmentation detection (connected component analysis)
- Track per-cluster mass, centroid, aspect ratio — not just global totals
- Time-averaged death check (rolling window, not instantaneous)
- Absorbing boundary damping for compact objects (not periodic-only)
