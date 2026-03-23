# SFA — SCP Field Archive

Single-file compressed container for volumetric field simulation data.

## Directory Structure

```
sfa/
  format/
    sfa.h           — Single-header C library (write + read + KVMD metadata)
    sfa_test.c      — Roundtrip test
  sim/
    scp_sim.c       — Unified CPU simulation kernel (OpenMP)
    scp_sim.cu      — Unified GPU simulation kernel (CUDA)
    braid_default.cfg — Standard braid config
    verify_parity.cfg — CPU/GPU parity test config
  seed/
    gen_braid.c     — Braid seed generator → single-frame SFA
    gen_oscillon.c  — Oscillon seed generator → single-frame SFA
  viewer/
    volview.c       — Interactive SFA volume viewer (SDL2, raymarcher)
    volview_bin.c   — Legacy .bin file viewer (single snapshot)
```

## Build

```bash
# Simulation kernel (CPU)
cd sim && gcc -O3 -march=native -fopenmp -o scp_sim scp_sim.c -lzstd -lm

# Simulation kernel (GPU)
cd sim && nvcc -O3 -arch=sm_70 -o scp_sim_cuda scp_sim.cu -lzstd -lm

# Seed generators
cd seed && gcc -O3 -o gen_braid gen_braid.c -lzstd -lm
cd seed && gcc -O3 -o gen_oscillon gen_oscillon.c -lzstd -lm

# SFA volume viewer
cd viewer && gcc -O3 -fopenmp -o volview volview.c \
    -I../format $(sdl2-config --cflags --libs) -lzstd -lm

# Library test
cd format && gcc -O2 -o sfa_test sfa_test.c -lzstd -lm && ./sfa_test
```

## Simulation Workflow

### 1. Generate a seed

```bash
# Standard braid at N=128
sfa/seed/gen_braid -N 128 -L 15 -o seed.sfa

# Braid at custom position
sfa/seed/gen_braid -N 128 -L 15 -cx 10 -cy 0 -cz 0 -o offset_braid.sfa

# Oscillon blob
sfa/seed/gen_oscillon -N 128 -L 10 -A 0.8 -sigma 3.0 -o oscillon.sfa
```

Seed generators write single-frame 12-column SFA files (6 fields + 6 velocities)
with physics parameters embedded as KVMD metadata.

### 2. Run simulation

```bash
# From config file
OMP_NUM_THREADS=8 sfa/sim/scp_sim sim/braid_default.cfg

# From seed SFA (parameters loaded from KVMD metadata, no config needed)
sfa/sim/scp_sim seed.sfa -T 200 -output run.sfa

# With CLI overrides
sfa/sim/scp_sim sim/braid_default.cfg -N 256 -T 500 -precision f64

# On GPU
sfa/sim/scp_sim_cuda sim/braid_default.cfg -output gpu_run.sfa
```

### 3. Resume / restart

```bash
# Continue from last frame (parameters auto-loaded from KVMD)
sfa/sim/scp_sim run.sfa -T 100 -output continued.sfa

# Continue from specific frame
sfa/sim/scp_sim run.sfa -init_frame 10 -T 50 -output branch.sfa

# CPU output → GPU restart (cross-platform, verified identical physics)
sfa/sim/scp_sim_cuda run.sfa -T 500 -output gpu_continued.sfa
```

### 4. View results

```bash
sfa/viewer/volview run.sfa

# Controls: mouse=rotate, scroll=zoom, left/right=frames, space=play
# 1/2/3=toggle channels, O/P=opacity, +/-=brightness, S=screenshot, Q=quit
```

## SFA File Format

12-column layout (restartable):

| Column | Name | Semantic | Component | Description |
|--------|------|----------|-----------|-------------|
| 0-2 | phi_x/y/z | POSITION | 0-2 | Displacement fields |
| 3-5 | theta_x/y/z | ANGLE | 0-2 | Rotation fields |
| 6-8 | phi_vx/vy/vz | VELOCITY | 0-2 | Displacement velocities |
| 9-11 | theta_vx/vy/vz | VELOCITY | 3-5 | Rotation velocities |

Output precision: f16, f32 (default), or f64. Compression: BSS + zstd (lossless on stored dtype).

### KVMD Metadata

SFA files embed simulation parameters as key-value metadata (KVMD chunks).
This makes files self-describing — `scp_sim file.sfa` loads all parameters
from the file itself. Priority: CLI args > config file > KVMD metadata > defaults.

Multiple KVMD sets per file (with frame ranges) support multi-parameter chains.

## Config File Format

Plain key=value text. Lines starting with `#` are comments.

```ini
# Grid
N          = 128
L          = 15.0
T          = 200.0
dt_factor  = 0.025

# Physics (6-field Cosserat)
m          = 1.5          # phi mass (m^2 = 2.25)
m_theta    = 0.0          # theta mass (massless)
eta        = 0.5          # phi-theta curl coupling
mu         = -41.345      # V(P) coefficient
kappa      = 50.0         # V(P) saturation

# Mass coupling mode: 0=constant, 1=inverse, 3=density-kappa
mode       = 0
kappa_gamma = 2.0         # mode 3 only

# Absorbing boundary (spherical)
damp_width = 3.0
damp_rate  = 0.01

# Initialization: braid | oscillon | sfa | exec
init       = braid
A          = 0.8
A_bg       = 0.1
R_tube     = 3.0
ellip      = 0.3325
delta      = 0,3.0005,4.4325

# SFA restart (when init=sfa)
init_sfa   = prev_run.sfa
init_frame = -1           # -1 = last frame

# Output
output     = output.sfa
precision  = f32          # f16 | f32 | f64
snap_dt    = 5.0
diag_dt    = 2.0
diag_file  = diag.tsv
```

## CPU/GPU Parity

The CPU (OpenMP) and GPU (CUDA) kernels produce identical physics.
Diagnostics match byte-for-byte. SFA frame data matches to machine
epsilon (~10⁻¹⁵) due to FMA instruction differences between x86 and CUDA.
Verified on RTX 4090 (SM 8.9) vs 16-core x86.

## Dependencies

- **libzstd** — compression (required for all tools)
- **SDL2** — viewer only
- **CUDA toolkit** — GPU kernel only
- No other external dependencies
