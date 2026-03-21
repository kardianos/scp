# V36: GPU Implementation + Multi-Atom Simulation

## IMPORTANT: Data Format Policy

**All simulation output from V36 onward MUST use the SFA format.**
No more raw .bin files. This ensures:
- Compression (BSS+zstd, 5× on active data)
- Seekable frame access
- Self-describing schema
- Streaming support for remote GPU runs
- Future AMR multi-grid compatibility

The SFA library is at `/home/d/code/scp/sfa/format/sfa.h`.

---

## Phase 1: CUDA Port (single V100)

Port the 6-field Cosserat force computation to CUDA.

### What gets ported
- `compute_forces()`: one CUDA thread per grid point
  - Laplacian (6-neighbor stencil)
  - Triple product V'(P)
  - Curl coupling η×curl(θ)
- `verlet_step()`: trivially parallel (per-point update)

### What stays on CPU
- Initialization (one-time)
- Diagnostics (infrequent, small)
- SFA compression + writing (async thread)

### Architecture
```
GPU memory:
  phi[3][N³]      — position fields
  vel_phi[3][N³]  — velocities
  acc_phi[3][N³]  — accelerations
  theta[3][N³]    — angle fields
  vel_theta[3][N³]
  acc_theta[3][N³]
  Total: 18 × N³ × 8 bytes

CPU memory:
  SFA write buffer (one frame)
  Async compression thread
```

### Kernel design
```cuda
__global__ void compute_forces_kernel(
    double *phi0, double *phi1, double *phi2,
    double *theta0, double *theta1, double *theta2,
    double *acc_phi0, ... , double *acc_theta2,
    int N, double idx2, double idx1,
    double MASS2, double MTHETA2, double MU, double KAPPA, double ETA) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N*N*N) return;

    int i = idx / (N*N), j = (idx/N) % N, k = idx % N;
    // Periodic neighbors
    int ip = (i+1)%N, im = (i-1+N)%N;
    // ... same force computation as CPU ...
}
```

### File structure
```
v36/
  src/
    cosserat_cuda.cu     — CUDA simulation (GPU kernel + host code)
    cosserat_cpu.c       — CPU fallback (for local testing)
    sfa_writer.c         — async SFA write thread
  scripts/
    vast_deploy.sh       — spin up, upload, run, download, destroy
    vast_search.sh       — find best V100 offer
  PLAN.md               — this file
```

### Build (on Vast instance)
```bash
nvcc -O3 -arch=sm_70 -o cosserat_gpu src/cosserat_cuda.cu -lzstd -lm -lpthread
```

### Disk space estimate (per run)
- Source: < 1 MB
- Binary: < 1 MB
- SFA output (N=256, 100 frames, f32 lossy): ~300 MB
- SFA output (N=256, 100 frames, lossless): ~1 GB
- Wedge data: < 10 MB
- **Total per run: ~1-2 GB**
- **Recommended Vast disk: 20 GB** (room for multiple runs)

---

## Phase 2: Validation

Reproduce V34 results on GPU:
1. Single braid stability (N=128, T=300): compare E_pot, θ_rms to CPU
2. Two-braid force (N=128, D=15): compare ΔD to CPU
3. Winding reversal: confirm charge-dependent force on GPU
4. Verify bitwise reproducibility is NOT required (GPU float arithmetic
   differs from CPU), but physics must match within 1%

---

## Phase 3: Multi-Scale on GPU

Combined core + wedge on a single V100:
- Core: N=256 (2.4 GB) — double the resolution of CPU runs
- Wedge: Nr=500 (CPU, negligible)
- SFA output: streaming to B2 bucket
- Target: 10× faster than current CPU → T=3000 in 30 minutes

---

## Phase 4: Multi-Atom (x8 V100 campaign)

Two hydrogen atoms interacting:
- 2× core braids (N=256 each, 2.4 GB each)
- 2× electron wedges
- Shared coarse interaction grid (N=128, 300 MB)
- Total GPU memory: ~6 GB (fits on one V100, use 8 for parameter sweep)

Physics targets:
- Do same-winding braids form a molecular bond through θ overlap?
- Does the bond length match molecular hydrogen?
- Can we observe orbital hybridization (sp mixing)?

---

## Vast.ai Deployment

### Search for V100
```bash
vastai search offers 'gpu_name=V100 num_gpus=1 disk_space>=20 inet_down>=200' \
  --order 'dph_total'
```

### Create instance
```bash
vastai create instance <OFFER_ID> \
  --image nvidia/cuda:12.2.0-devel-ubuntu22.04 \
  --disk 20 \
  --ssh
```

### Upload + compile + run
```bash
# Upload source (< 100 KB)
scp -P PORT src/*.cu src/*.c sfa.h scripts/run.sh user@HOST:~/

# SSH in
ssh -p PORT user@HOST
apt-get update && apt-get install -y libzstd-dev
nvcc -O3 -arch=sm_70 -o sim src/cosserat_cuda.cu -lzstd -lm
./sim -N 256 -T 300 -o output.sfa

# Download results
scp -P PORT user@HOST:~/output.sfa ./
# Or sync to B2:
# b2 sync ~/results/ b2://bucket/v36/
```

### Destroy when done
```bash
vastai destroy instance <INSTANCE_ID>
```

### Cost estimate
- V100 × 1: ~$0.20/hr
- Typical session: 15-30 min = $0.05-0.10
- Full V36 development: ~$5-10
- x8 campaign: ~$2-5
- **Total budget needed: ~$15 of $25 available**
