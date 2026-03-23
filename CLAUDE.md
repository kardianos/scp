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

## Remote GPU Policy
- **ALL remote simulations MUST run on the GPU, NOT the CPU**
- If running on Vast.ai or any remote GPU server, use CUDA or JAX — never OpenMP C
- Running C/OpenMP on a GPU server wastes money and time (same speed as local CPU)
- The CUDA simulator (`sfa/sim/scp_sim.cu`) runs 10-13x faster than CPU
- If a CUDA version doesn't exist for a tool, BUILD ONE before deploying remotely
- Monitor GPU utilization (`nvidia-smi`) every 10s to ensure GPU is actually being used

## Build Convention
- C with OpenMP: `gcc -O3 -march=native -fopenmp -o <binary> src/<file>.c -lzstd -lm`
- **ONE copy of sfa.h** at `/home/d/code/scp/sfa/format/sfa.h` — do NOT copy to other directories
- Include via relative path: `#include "../../sfa/format/sfa.h"` (adjust depth as needed)
- No external deps beyond zstd, OpenMP, math

## Physics Requirements — CRITICAL
- **ALL simulations MUST use the FULL 6-field Cosserat equation (3 phi + 3 theta)**
- The equation: d²φ/dt² = ∇²φ - m²φ - V'(P) + η×curl(θ), d²θ/dt² = ∇²θ - m_θ²θ + η×curl(φ)
- NEVER build a 3-field-only simulator. The theta coupling is ESSENTIAL physics.
- Default: m²=2.25, m_θ²=0 (massless theta), η=0.5, μ=-41.345, κ=50
- Verify: 18 arrays (not 9), curl terms in forces, theta_rms grows from zero
- Reference implementation: `v34/torsion_coupling/src/v33_cosserat.c`

## Simulation Kernel Policy — CRITICAL
- **The unified sim kernel (`sfa/sim/scp_sim.c` and `sfa/sim/scp_sim.cu`) MUST NOT be modified unless the user EXPLICITLY requests changes to the simulation code.**
- "Explicitly" means the user says something like "modify scp_sim", "change the kernel", "update the simulator code", "add X to scp_sim". General requests like "run a simulation", "test this parameter", "try a new experiment" do NOT authorize modifying the kernel.
- To run experiments: write a config file and/or a seed generator. The kernel reads config — never modify it for a specific experiment.
- To add new initialization patterns: write a new seed generator in `sfa/seed/`, NOT by editing the kernel.
- The same protection applies to `sfa/format/sfa.h` — do NOT modify the SFA format without explicit user authorization.
- If a simulation requires physics not supported by the current kernel, ASK the user before modifying it. Describe what change is needed and why.

## Simulation Kernel Location
- CPU kernel: `sfa/sim/scp_sim.c` — build with `gcc -O3 -march=native -fopenmp -o scp_sim scp_sim.c -lzstd -lm`
- GPU kernel: `sfa/sim/scp_sim.cu` — build with `nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm`
- Seed generators: `sfa/seed/gen_braid.c`, `sfa/seed/gen_oscillon.c`
- Config files: `sfa/sim/*.cfg`
- SFA header: `sfa/format/sfa.h` (single copy, include via relative path `../format/sfa.h`)
- Reference implementation (historical): `v34/torsion_coupling/src/v33_cosserat.c`

## Diagnostics Requirements
- Every simulation MUST include fragmentation detection (connected component analysis)
- Track per-cluster mass, centroid, aspect ratio — not just global totals
- Time-averaged death check (rolling window, not instantaneous)
- Absorbing boundary damping for compact objects (not periodic-only)
