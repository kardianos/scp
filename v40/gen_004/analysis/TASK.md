# Gen 4 Structure Characterization Task

## Goal

Numerically characterize WHAT makes stable vs unstable regions of the T=200
survivor structures. For each surviving structure (S20, UDD_R4, CB15), analyze
the field properties in regions that persist vs regions that fragment/disperse.

## Specific Questions

1. **What does ρ (field energy density) look like in stable vs unstable regions?**
   - Radial profile of ρ around cluster centroids vs inter-cluster voids
   - Is there a critical ρ threshold below which regions fragment?

2. **What does θ look like in stable vs unstable clusters?**
   - Is θ_rms higher in surviving clusters than in dispersing ones?
   - Does the θ field have coherent structure (aligned curl) or random?
   - Does θ_rms correlate with binding longevity?

3. **What do velocities tell us?**
   - Are surviving clusters expanding or contracting (div(v) sign)?
   - Is there systematic drift (whole cluster moving) vs oscillation?
   - Kinetic energy concentration: is KE localized in cores or shells?

4. **What does the acceleration field reveal?**
   - Where are the strongest forces? Core, surface, or inter-cluster?
   - Is the V(P) force dominant or is the curl coupling force significant?
   - Do stable regions have balanced forces (near-zero net) vs unstable?

5. **Can we define a local stability metric?**
   - Per-voxel or per-cluster stability score based on (ρ, θ, v, a)
   - Predict from t=0 data which regions will survive to t=200

## Available Data

- `gen_004/S20_output.sfa` — 11 GB, 42 frames, N=192, 12 columns (phi/theta/vel)
- `gen_004/UDD_R4_output.sfa` — 6 GB, N=192 (frames after t=100 may be duplicated)
- `gen_004/CB15_output.sfa` — 11 GB, 42 frames, N=192, 12 columns
- Spatial analysis already done (cluster counts, centroids, aspect ratios)

## Tools Available

- `v40/tools/analyze_sfa` — frame-by-frame global metrics
- `v40/tools/spatial_analysis` — cluster detection, centroids, R_rms
- `sfa/format/sfa.h` — SFA reader API (sfa_open, sfa_read_frame)
- Custom C tools can be written in `v40/gen_004/analysis/`

## Constraints

- Do NOT modify `sfa/sim/scp_sim.c`, `sfa/sim/scp_sim.cu`, or `sfa/format/sfa.h`
- New analysis tools go in `v40/gen_004/analysis/`
- Build: `gcc -O3 -fopenmp -o tool tool.c -I/home/d/code/scp -lzstd -lm`
  (use `#define SFA_IMPLEMENTATION` and `#include "sfa/format/sfa.h"`)
