# SCP Project Conventions

## Document Roles — CRITICAL
- **CONCEPT.md** is the THEORY DOCUMENT — it presents the current best understanding
  of the physics as a cohesive, replicable description. Write it like a textbook chapter,
  NOT a lab notebook. Someone reading only CONCEPT.md should be able to understand and
  replicate the theory without knowing about the trial and error that got there.
  - Do NOT include chronological experiment narratives or failed approaches
  - DO include null results that clarify the theory ("it is NOT X, it IS Y")
  - DO include key numerical data that confirms claims
  - Present the physics as understood NOW, not as it was discovered over time
  - Use correct terminology: "particle" or "proton" for gravitationally responsive
    objects, "braid" only for the z-aligned sub-component (quark analog)
- **DISCOVERIES.md** is the chronological lab notebook — records what was found, when,
  including failed approaches and process. This is where the history lives.
- **EM_THEORY.md** is the detailed electromagnetic sector theory document (same
  standards as CONCEPT.md — cohesive, not chronological).
- **FUTURE.md** tracks open questions and proposed experiments.
- **Version directories** (v28/, v34/, v41/, v42/, v43/) contain per-experiment
  plans, results, analysis, and generated data. Each version builds on the previous.

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

## Running Simulations — Use the MCP Runner

The `scp-runner` MCP server (source: `sfa/runner/`, binary: `bin/scp-runner`) manages ALL simulation
execution, both local and remote. It exposes `sim_*` tools via MCP that handle
instance management, building, running, monitoring, and downloading — with NO
sleep polling required. Use these tools instead of manual SSH/SCP/rsync.

**IMPORTANT: Only one executor is active at a time.** Calling `sim_setup` replaces
the current executor. To run local then remote, complete the local workflow first,
then call `sim_setup(executor="remote")` to switch.

**Local execution** (CPU, for quick tests):
```
sim_setup(executor="local")
sim_build(sources=["sfa/sim/scp_sim.c"])          ← auto-detects gcc
sim_run(config="N=64\nL=15\nT=5\n...", id="test_001")  ← config content, uses last-built binary
sim_run_status(id="test_001")   ← instant, no polling
```

**Remote execution** (GPU, for production runs):
```
sim_setup(executor="remote")     ← auto-provisions Vast.ai instance (reuses existing if running)
sim_build(sources=["sfa/sim/scp_sim.cu", "sfa/format/sfa.h"])  ← auto-detects nvcc
sim_run(config="N=384\nL=100\nT=200\n...", id="gradient_test")
sim_run_status(id="gradient_test")
sim_download(remote_path="output.sfa", local_path="/space/scp/results/")
sim_teardown()                   ← verifies downloads, destroys instance
```

**Remote with existing instance** (connect to already-running GPU):
```
sim_setup(executor="remote", host="ssh5.vast.ai", port=12345)
```

**Custom GPU filter** (default is V100 16GB):
```
sim_setup(executor="remote", gpu_filter="gpu_name=RTX_4090 num_gpus=1 rentable=true disk_space>=20")
```

**sim_run config**: Pass config file content (key=value lines). The runner writes it to
a `.cfg` file and invokes the last-built binary automatically. Alternatively, pass a
command string (e.g. `/path/to/binary /path/to/config.cfg`) for direct execution.

**sim_build**: Omit `cmd` for auto-detection (gcc for .c, nvcc for .cu). If you need
a custom build command, use `${OUTPUT}` as the output path placeholder.

**DO NOT** use manual SSH, SCP, rsync, or `sleep N` polling for simulation work.
The runner handles all of this internally with goroutines — every tool call
returns instantly with cached state.

**Key tools**: `sim_setup`, `sim_status`, `sim_build`, `sim_run`, `sim_run_status`,
`sim_run_cancel`, `sim_upload`, `sim_download`, `sim_download_status`,
`sim_list_files`, `sim_exec`, `sim_teardown`

## GPU Notes
- **V100-16GB**: Fits N=384 (10.3 GB). Use for standard proton/gradient tests.
- **V100-32GB**: Fits N=512 (19.3 GB). Use for large-grid or high-resolution runs.
- **RTX 4090**: 24 GB, sm_89, faster single-GPU but no multi-GPU.
- The CUDA kernel (`sfa/sim/scp_sim.cu`) uses async hooks for I/O — snapshots
  and diagnostics overlap with physics compute. GPU utilization should be ~100%.
- The `init=template` mode loads a small proton template (5 MB) and generates
  the background analytically — no large seed files needed.
- For gradient tests, use `bc_type=1` with `gradient_A_high`/`gradient_A_low`.

## Build Convention
- **`make -C sfa install`** builds everything and installs to `bin/` (18 binaries)
- **`make -C sfa runner`** builds just the MCP runner
- **`make -C sfa analysis`** builds just the analysis tools
- **`make -C sfa cuda`** builds CUDA kernels (requires nvcc)
- C with OpenMP: `gcc -O3 -march=native -fopenmp -o <binary> src/<file>.c -lzstd -lm`
- CUDA: `nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm -lpthread`
- **ONE copy of sfa.h** at `/home/d/code/scp/sfa/format/sfa.h` — do NOT copy to other directories
- Include via relative path: `#include "../../sfa/format/sfa.h"` (adjust depth as needed)
- All binaries go to `bin/` via `make install` — do NOT commit binaries

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
- Seed generators: `sfa/seed/gen_braid.c`, `sfa/seed/gen_oscillon.c`,
  `sfa/seed/gen_phase_confined.c` (UUD/UDD composites),
  `sfa/seed/gen_proton_analytical.c` (algebraic proton/neutron),
  `sfa/seed/gen_composite.c` (stamp templates into grids for composite nuclei)
- Config files: `sfa/sim/*.cfg`
- SFA header: `sfa/format/sfa.h` (single copy, include via relative path `../format/sfa.h`)
- Reference implementation (historical): `v34/torsion_coupling/src/v33_cosserat.c`

## Data Transfer and Storage
- The `scp-runner` handles all remote file transfers via `sim_download`.
  It uses rsync with `--append-verify` internally — no manual rsync needed.
- `sim_teardown` verifies downloads before destroying instances.
- Large SFA files should go to `/space/scp/` (separate disk, 600+ GB free).
- Local working disk (`/home/d/code/scp/`) has limited space — keep only
  small files (diag.tsv, templates, analysis results) on the local disk.

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
