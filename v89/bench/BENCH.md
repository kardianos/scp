# The speedup ladder — measured (ratchet-gated at every rung)

Protocol: `./bench.sh <label>` — 3 repeats per box, median wall time;
boxes: bench_blob (L=24 g4-class: blob + footprint + quiet vacuum,
9741 cells / 85886 links) and bench_noise (L=20 hot noise: everything
active). Laws: the standing table. Every rung passes the FULL battery
before its numbers are recorded here.

| rung | change | blob ms/step | noise ms/step | cum. speedup (blob) |
|---|---|---|---|---|
| 0 | baseline (single thread, -O2 -march=native) | 18.30 | 6.25 | 1.0× |
| 1 | OpenMP, 8 threads (deterministic gather refactor + edge-colored hops + serial RNG draws) | 7.45 | 2.96 | **2.46×** |
| 1 | — same binary, 1 thread | 32.28 | 10.71 | 0.57× (the gather double-visits + snapshots cost ~1.8× serial; parallelism pays it back) |
| 1 | — same binary, 16 threads | 15.13 | 10.11 | 1.2× (hyperthread collapse: 8 physical cores; use OMP_NUM_THREADS=8) |

Rung-1 notes: results are **thread-count independent to the byte**
(verified: 1-thread vs 16-thread logs identical) — gathers and color
batches have fixed evaluation order, the tumble gaussian stream is
drawn serially in the legacy order, instrument/conversion prints stay
serial. Numerics differ from the pre-refactor kernel at roundoff
(gather summation order, Jacobi space limiter, colored hop order,
snapshot clocks) — full battery is the judge. Host: 8 physical cores /
16 HT. Battery runs now use jobs×threads = cores (threads =
nproc/jobs).
