# V38 Experiment Plan: Braid-as-Electron Characterization

## Goal

Characterize the helical braid from V33/V34 at high resolution on GPU.
Determine whether the theta field is truly 1/r (Maxwell-like) and whether
multi-braid bound states exist.

## What Was Missing in V37

V37 ran CPU simulations on GPU servers (wasting money). V38 fixes this
with proper CUDA implementations:
- `seedrun_cuda.cu`: GPU seed runner with absorbing BC, f32 SFA output
- `braid_analyze_cuda.cu`: GPU braid characterizer with radial theta profiling

## Budget

8 hours of V100 time at ~$0.12/hr = ~$1.00 total.

## Phase 1: Single Braid Characterization (N=256, est. 30 min)

**Objective**: Measure the theta field profile around the braid axis with
enough resolution to distinguish 1/r from exponential or other fall-offs.

**Parameters**:
- N=256, L=20, T=200 (dx=0.157, 12.7 pts across braid core)
- m^2=2.25, eta=0.5, m_theta=0 (massless theta)
- Standard braid init: A=0.8, delta={0, 3.0005, 4.4325}

**Measurements** (every 10 timesteps):
- Azimuthally-averaged |theta|(r) at z=0 plane
- Comparison with 1/r Biot-Savart prediction
- Winding number of phi field
- Full energy decomposition (9 components)
- Braid z-centroid velocity

**Expected output**:
- `output/braid_n256/sim.sfa` (~2 GB)
- `output/braid_n256/radial_profiles.tsv`
- `output/braid_n256/timeseries.tsv`

**Success criteria**:
- theta(r) follows 1/r for r > 2 R_tube (Biot-Savart regime)
- Winding number stable at +/-1
- Energy conservation < 0.1%

**Command**:
```
./braid_analyze -N 256 -L 20 -T 200 -eta 0.5 -mt 0 -m 1.5 -snap 10 -aevery 10 -o output/braid_n256
```

## Phase 2: Two-Braid Interaction (N=256, est. 60 min)

**Objective**: Determine whether two braids attract, repel, or bind.
Run two braids separated by D=12 (just outside overlap) and watch.

**Parameters**:
- N=256, L=20, T=300
- Two braids at x=-6, x=+6 (same winding)
- Same physics as Phase 1

**Measurements**:
- Same as Phase 1 plus inter-braid separation vs time
- Check for bound state oscillation or escape

**Expected outcomes**:
A. Same-winding: attractive force (V34 result). May form bound state.
B. If they merge, the resulting structure is interesting (meson analog?).

**Command**:
```
./braid_analyze -N 256 -L 20 -T 300 -eta 0.5 -mt 0 -m 1.5 -braids 2 -D 12 -snap 10 -aevery 10 -o output/two_braid_n256
```

## Phase 3: Evolutionary Search at N=48 on GPU (JAX)

**Objective**: Use the V37 JAX evolutionary pipeline but run ON GPU.
Search for compact structures that survive T=500.

**Note**: This phase uses the existing JAX `field_evolve.py` from V37,
not a new CUDA binary. The JAX code already runs on GPU if available.

**Parameters**:
- N=48, L=8 (compact search domain)
- Population=64, generations=200
- CMA-ES with survival fitness

**Estimated time**: 2-3 hours (batched on GPU, much faster than V37 CPU).

**This phase is lower priority** -- only run if Phase 1+2 give clear results.

## Phase 4: Validate Survivors at N=128/256 (seed runner)

**Objective**: Take any interesting structures from Phase 3 (or from V37
survivors) and run them at high resolution with `seedrun_cuda`.

**Parameters**:
- N=128 or N=256, L=15, T=500
- Absorbing BC: damp_width=4, damp_rate=0.005
- f32 SFA output every 5 time units

**Command**:
```
./seedrun_cuda -seed <file.bin> -N 256 -L 15 -T 500 -snap 5 -damp_width 4 -damp_rate 0.005 -o output/seed_validate
```

## Deployment

```bash
# From local machine:
cd /home/d/code/scp/v38

# Phase 1 (single braid):
./scripts/deploy.sh user@gpu-host 1

# Phase 2 (two braids):
./scripts/deploy.sh user@gpu-host 2

# Phase 4 (seed validation):
./scripts/deploy.sh user@gpu-host 4 /path/to/seed.bin
```

The deploy script handles: upload, compile, run, download results.
The monitor script runs in background checking GPU utilization every 10s.

## Key Questions to Answer

1. **Is theta ~ 1/r?** If yes, the braid generates a Coulomb-like field
   from the curl coupling. This is the electromagnetic field analog.

2. **What is the braid's "charge"?** The coefficient of the 1/r tail
   gives an effective charge Q_eff. Compare with electron charge.

3. **Do same-sign braids bind?** If two same-winding braids attract
   and form a bound state, this could be a meson or baryon analog.

4. **Is the theta response truly spin-1?** The theta field has 3 components.
   Check if the radial profile decomposes into l=1 (dipole) pattern.

5. **What sets the braid velocity?** The traveling wave reconstructs at
   each step. Is the velocity a function of amplitude, or is it quantized?

## Memory Budget (N=256 on V100 16GB)

- 18 arrays x 256^3 x 8 bytes = 2.42 GB (f64 fields)
- f32 buffer: 256^3 x 4 bytes = 67 MB
- Total GPU: ~2.5 GB (fits easily on V100)
- SFA output at snap_dt=10: ~200 frames x 6 x 256^3 x 4 = ~12 GB disk
  (with zstd compression: ~3-4 GB)

## N=512 Stretch Goal

If Phase 1 at N=256 shows interesting structure at the core scale,
run at N=512 (dx=0.078, 25 pts across core).

Memory: 18 x 512^3 x 8 = 19.3 GB. Requires V100 32GB or A100.
Time: ~4x longer than N=256, roughly 2 hours.

Only attempt if budget allows and Phase 1 results warrant it.
