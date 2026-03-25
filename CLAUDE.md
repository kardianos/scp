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
- **Vast.ai GPU search**: use `gpu_name=Tesla_V100` (NOT `V100`).
  Two variants: Tesla_V100 16 GB and Tesla_V100 32 GB (SXM2). Both sm_70.
  Multi-GPU: available as 4× and 8× for very large grids (N=512+).
  Search: `vastai search offers 'gpu_name=Tesla_V100 num_gpus=1 rentable=true disk_space>=20' -o 'dph'`
  Multi: `vastai search offers 'gpu_name=Tesla_V100 num_gpus>=4 rentable=true' -o 'dph'`
  Pricing: $0.12-0.20/hr per GPU. 4×V100 ~$0.50-0.80/hr.
  RTX 4090 ($0.30/hr, sm_89, 24 GB VRAM) is faster for single-GPU but no multi-GPU.
  Compile for Tesla_V100: `nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm`
  Compile for RTX 4090: `nvcc -O3 -arch=sm_89 -o scp_sim_cuda scp_sim.cu -lzstd -lm`

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

## Remote Simulation Data Transfer Policy
- **Use rsync with `--append-verify` to incrementally download SFA files DURING simulation**
  ```bash
  rsync -avz --partial --append-verify \
      -e 'ssh -o StrictHostKeyChecking=no -p PORT' \
      root@HOST:/root/output.sfa ./local_results/
  ```
  This downloads new frames as they are written, with automatic resume on connection drop.
  Run every 10 minutes via a monitor agent or background loop.

- **NEVER destroy a GPU instance before verifying ALL downloads are complete**
  Verification protocol:
  1. Run analysis and f16 conversion ON THE REMOTE before downloading
  2. Download all files (SFA, JSON, TSV)
  3. For EACH file, compare local size to remote size:
     ```bash
     remote_size=$(ssh remote "stat -c%s /root/file")
     local_size=$(stat -c%s ./file)
     [ "$remote_size" = "$local_size" ] && echo "VERIFIED" || echo "MISMATCH"
     ```
  4. Run `sfa_info` on downloaded SFA to check for truncation warnings
  5. ONLY after ALL files verified → destroy instance
  6. If any download fails, RETRY up to 3 times. Do NOT destroy on failure.
  See `sfa/sim/REMOTE_PROTOCOL.md` for the full protocol.

## SFA Archival (rclone)
- Completed SFA files can be archived to cloud storage via rclone:
  ```bash
  rclone copy local_file.sfa scpsfa:scpsfa/v42/
  rclone ls scpsfa:scpsfa/           # list archived files
  rclone copy scpsfa:scpsfa/v42/file.sfa ./  # retrieve when needed
  ```
- The `scpsfa` rclone remote is pre-configured locally for B2 cloud storage.
- Archive SFA files after analysis is complete to free local disk space.
- Keep diag.tsv, analysis.json, and freq.json locally (small files, always needed).
- Large f32 output SFAs (10-30 GB) should be archived and deleted locally.
  The f16 viewing copies (1-9 GB) can be kept locally or also archived.
- To check if a file is already archived: `rclone ls scpsfa:scpsfa/`

## Diagnostics Requirements
- Every simulation MUST include fragmentation detection (connected component analysis)
- Track per-cluster mass, centroid, aspect ratio — not just global totals
- Time-averaged death check (rolling window, not instantaneous)
- Absorbing boundary damping for compact objects (not periodic-only)
