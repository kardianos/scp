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
  - Use correct terminology (U(1) era, v66+): "ball"/"baryon" for the stable
    charged Q-ball, "component"/"quark" for the Φ_a sub-fields, "flavor" for the
    frequency/charge partition (ω_a, Q_a), "nucleus" for fused multi-ball droplets.
    "Braid", "oscillon", "UUD/UDD proton/neutron" are real-field-era (≤v53) terms —
    use only in historical context.
- **DISCOVERIES.md** is the chronological lab notebook — records what was found, when,
  including failed approaches and process. This is where the history lives.
- **EM_THEORY.md** is the detailed electromagnetic sector theory document (same
  standards as CONCEPT.md — cohesive, not chronological). NOTE: the photon analog is
  the gauge field A (v68+); the old θ-polariton mapping is historical.
- **FUTURE.md** tracks open questions and proposed experiments.
- **Version directories** (v28/ … v71/) contain per-experiment plans, results,
  analysis, and generated data. Current era: v66 (complexified U(1) Q-balls),
  v67–v68 (characterization + gauge design), v69 (gauged kernel-v3), v70 (skeptical
  verification + existence dynamics), v71 (quark substructure, nuclei, flavor).

## Data Format Policy
- **ALL simulation output MUST use SFA format** (`.sfa` files via `sfa.h`)
- **NO `.bin` file output** — remove `save_field()` functions that write raw binary
- SFA header: `/home/d/code/scp/sfa/format/sfa.h` — single-header C library
- Usage: `#define SFA_IMPLEMENTATION` in exactly ONE .c file before including
- Link with: `-lzstd`
- Column convention: phi_x/phi_y/phi_z (SFA_POSITION), theta_x/theta_y/theta_z (SFA_ANGLE)
- All analysis tools must READ SFA files using `sfa_open()` / `sfa_read_frame()`
- SFA viewer: `sfa/volview` (binary `bin/volview`) — handles 12/24/30-column files;
  views: 4 field / 5 velocity / 6 accel / 8 U(1) gauge (|E|,|A|) / 9 charge (±ρ_Q) /
  0 flavor (per-component RGB + inline clock analysis) / C local clock;
  headless export: `volview -snapshot N -view K -out f.webp file.sfa`

## Simulation Parameters (standard)
- **Current standard (U(1) era, v66+)**: complexified 12-field Cosserat, optionally
  gauged. Config: `complex_phi=1`, `complex_gauge=1`, `g_gauge=0.05`, `m_theta=1.6`
  (closes the θ-drain), `eta=0` for particle experiments, absorbing BC.
  Vt(s) = (mu/2) s/(1+kappa s), s = Π|Φ_a|²; m²=2.25, mu=-41.345, kappa=50.
  Q-ball window ω ∈ (1.3087, 1.5) ungauged; (1.406, 1.5) and Q_max=921 at g=0.05.
  Seed via radial profiles (init=qball on CPU; gen_qball_* + init=sfa on GPU —
  the GPU kernel has no init=qball path; the init Gauss projection builds E).
- Real-kernel legacy defaults (≤v65): m^2=2.25, mu=-41.345, kappa=50, eta=0.5,
  A_bg=0.1, delta = {0, 3.0005, 4.4325} (from v28 optimization)

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
- **ALL simulations MUST use the unified kernel with the full field content** —
  never a reduced/3-field-only simulator. The theta sector and (in complex mode)
  the imaginary sector are essential physics even when seeded to zero.
- Real mode: d²φ/dt² = ∇²φ - m²φ - V'(P) + η×curl(θ), d²θ/dt² = ∇²θ - m_θ²θ + η×curl(φ)
- Complex/gauged mode (v66/THEORY.md, v68/GAUGE_DESIGN.md, v69/SPEC.md): all matter
  derivatives link-covariant; charge conservation and the discrete Gauss law are
  exact by construction — verify `gauss_max` stays at the 1e-13 floor (a drift is
  implementation-bug tripwire #1).
- Verify (complex): 36 arrays (+ gauge blocks when complex_gauge=1), Q_phi conserved
  to the integrator floor at eta=0, g=0 byte-identical to the ungauged path.
- Reference docs: `v66/THEORY.md` (complexification), `v69/SPEC.md` (gauged lattice);
  historical real-field reference: `v34/torsion_coupling/src/v33_cosserat.c`

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
- Seed generators (U(1) era): `gen_qball_pair.c` (two balls, phases/signs),
  `gen_qball_boost.c` (boosted single ball; v=0 → single static ball),
  `gen_qball_bath.c` (ball + θ bath), `gen_qball_multi.c` (N balls — nuclei),
  `gen_qball_quark.c` (single-COMPONENT lumps, mask + per-component centers),
  `gen_qball_flavored.c` (per-component profiles/frequencies; multi-ball;
  accepts 2-col symmetric or 4-col flavored profiles). Profiles from
  `radial_qball` (ungauged) or the v69 gauged shooter. NOTE: seed writers MUST
  call `sfa_finalize_header()` before `sfa_write_frame()`.
- Seed generators (real-field era, historical): `gen_braid.c`, `gen_oscillon.c`,
  `gen_phase_confined.c` (UUD/UDD), `gen_proton_analytical.c`, `gen_composite.c`
- **Analytical seed warning**: `gen_proton_analytical` (Level 2) produces baryons
  with ~4× the equilibrium binding density (P_int ~1270 vs ~320 per baryon in V42).
  These seeds need T=2000+ to relax and are NOT suitable for energy comparisons.
  For energy-sensitive experiments (mass defect, binding), use template seeding
  (`init=template` with `proton_template.sfa`) or the `gen_deuterium.c` generator
  which produces better-equilibrated initial conditions.
- Pre-converged templates: `v43/proton_formation/proton_template.sfa` (64³, UUD),
  `v43/proton_formation/neutron_template.sfa` (192³, UDD from V41)
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
