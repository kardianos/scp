# V50 Hardware Performance Notes

## Current Development System

    CPU: 8-core Intel (presumed), DDR4-3200, 2 channels
    Practical memory bandwidth: ~38 GB/s
    AVX2 (256-bit), no AVX-512
    L3 cache: ~16 MB

## Performance Characteristics (from V50/C4 analysis)

The SCP simulation is **memory-bandwidth bound**. Arithmetic intensity is
0.36 FLOPs/byte — 21× below the machine's compute ridge point. At N=100,
the working set is 144 MB (9× L3). The stencil accesses 6 neighbors per
voxel across 12+ arrays, creating scattered memory access patterns that
cannot be SIMD-vectorized.

Key measurements (N=100, current system):
- 4 threads saturate memory bandwidth (31 ms/step)
- 8 threads give no improvement (31 ms/step)
- Hyperthreading (16T) is 6× SLOWER (cache thrashing)
- The two-pass C4 force computation adds 45% overhead at 8T (extra memory traffic)
- Auto-vectorization fails on the stencil loop (scattered access)

## AMD EPYC 9255 (Turin/Zen 5) — Recommended Target

    Cores: 24 (48 threads)
    Base: 3.2 GHz, Boost: 4.3 GHz, All-core: 4.0 GHz
    Memory: 12-channel DDR5-6400 (614 GB/s theoretical, ~400 GB/s practical)
    L3 cache: 128 MB
    AVX-512: Native 512-bit execution (Zen 5)
    TDP: 200W (configurable 200-240W)
    Socket: SP5

### Why this CPU changes the equation

    Current BW:  38 GB/s → saturates at 4 cores, 8 cores idle
    EPYC 9255:  400 GB/s → saturates at ~16 cores, 8 cores spare

The 10× bandwidth increase means:
1. Physics stencil scales to ~16 cores (vs 4 currently)
2. ~8 spare cores available for parallel COLZSTD frame compression
3. Physics and I/O can genuinely overlap without competing for bandwidth
4. The two-pass overhead (C4 Cosserat + hardening) drops because the
   temporary arrays partially fit in the 128 MB L3

Native 512-bit AVX-512 benefits:
- Verlet integration loops (contiguous v += dt*a) at 8 doubles/instruction
- BSS byte-shuffle encoding for SFA output
- zstd internal SIMD in parallel compression threads
- The main stencil loop still won't vectorize (scattered access), but
  the surrounding code runs faster

### Estimated performance

    N=100:  ~3-4 ms/step  (vs 31 ms/step local, ~8× faster)
    N=384:  ~15-20 ms/step (vs ~190 ms/step local, ~10× faster)

    T=400 two-proton at N=100:   ~3-4 minutes  (vs ~35 min local)
    T=1000 production at N=384:  ~30-40 minutes (competitive with V100 GPU)

### NUMA considerations

EPYC uses chiplets (CCDs) with separate memory controllers. For optimal
performance, pin threads and memory to the same NUMA domain:

    OMP_PLACES=cores OMP_PROC_BIND=close numactl --localalloc ./scp_sim_c4 ...

Cross-NUMA access halves effective bandwidth. Always verify NUMA topology
with `numactl --hardware` on the target system.

### Cost comparison (speculative, bare metal cloud)

A bare metal EPYC 9255 instance at ~$1-2/hr would be cost-competitive
with a V100 GPU instance at ~$0.15-0.30/hr for short runs, considering:
- No CUDA porting needed (CPU kernel runs as-is)
- Burst writes 8× faster (parallel compression on spare cores)
- The C4 kernel with Cosserat + hardening is not yet ported to CUDA
- For N≤384 the CPU performance is within 2-3× of V100

For N>512 or production runs >T=5000, GPU (V100/A100) remains faster
due to HBM bandwidth (~900-2000 GB/s). But for development iteration
and parameter sweeps at N=100, the EPYC is ideal.

## Other options considered

    EPYC 9124 (Genoa/Zen 4): 16 cores, DDR5-4800, split AVX-512
      - 25% less bandwidth than Turin
      - AVX-512 at 256-bit execution (half speed)
      - Fine if cheaper; ~20% slower than 9255

    EPYC 9135 (Turin/Zen 5): 8 cores, DDR5-6000
      - Same per-core performance as 9255
      - Only 8 cores: saturates bandwidth with no spare for I/O
      - No headroom for parallel compression during bursts

    GPU (V100/A100): HBM bandwidth dominates for large grids
      - Requires CUDA port of C4 kernel (not yet done)
      - Better for N≥512 production runs
      - Worse for rapid iteration (compile+test cycle slower)
