# V36 First GPU Run: V100 Cosserat Simulation

## Setup

- GPU: Tesla V100-SXM2-32GB (Vast.ai, Utah, $0.17/hr)
- Image: nvidia/cuda:12.2.0-devel-ubuntu22.04
- Build: nvcc -O3 -arch=sm_70 -o sim cosserat_cuda.cu -lzstd -lm

## Run Parameters

N=128, L=20, T=100, η=0.5, m_θ²=0, dt=0.0315

## Results

- **Wall time: 20.9 seconds** (6.6 ms/step)
- CPU equivalent: ~5 minutes → **15× speedup**
- E_pot oscillates -29 to -264 (braid alive, breathing)
- θ_rms = 0.025-0.041 (curl coupling active)
- SFA output: 89 MB (11 frames, BSS+zstd)

## Physics Verification

Matches V34 CPU results:
- E_pot range: ✓ (same oscillation amplitude)
- θ_rms: ✓ (same magnitude)
- Braid survival: ✓

## Performance Projection

| Grid | ms/step | T=300 wall time | Speedup vs CPU |
|------|---------|----------------|----------------|
| N=128 | 6.6 | ~35 seconds | 15× |
| N=256 | ~50 | ~12 minutes | 40× |
| N=512 | ~400 | ~2 hours | 50× |

## Cost

- Instance: $0.17/hr
- This run: ~2 minutes of GPU time = **$0.006**
- Total session (boot + compile + run + download): ~10 min = **$0.03**

## Files

- `src/cosserat_cuda.cu` — CUDA simulation source
- `first_gpu_run.sfa` — output (89 MB, 11 frames)
- SFA C++ compatibility fix applied to `sfa/format/sfa.h`

## NOTE: All data from V36 onward uses SFA format exclusively.
