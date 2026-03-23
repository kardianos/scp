# V38: Discover a Stable, Bounded, Non-Axis-Aligned Structure

## Goal

Find a self-sustaining, localized field configuration analogous to a proton
or neutron. Must be:
- **Bounded** (doesn't extend to infinity, no periodic BC dependence)
- **Non-axis-aligned** (3D structure, not a 1D tube)
- **Self-sustaining** (survives T≥200 with absorbing BC)
- **Carries charge** (winding number, chirality)

## Approach: Toroidal Braid

The helical braid self-sustains because the traveling wave continuously
reconstructs V(P) binding. It fails as a particle because it extends
forever along one axis.

**Solution**: Bend the braid into a torus. The wave travels around the
ring continuously — no endpoints, no corners, no redirection. Naturally
bounded and 3D.

Key insight: the V37 "ring" tests failed because they put the three fields
on SEPARATE rings (co-location violated). The correct approach: ALL three
fields on ONE ring with phase offsets along the circumference, exactly like
the working braid but closed.

### Closure condition
The braid oscillates as cos(k*s + delta_a) along arc length s.
For a torus of circumference C = 2*pi*R_major:
- Need k*C = 2*pi*n (integer n oscillations around the ring)
- k = n/R_major
- For n=3 (three oscillations, matching the 3-fold phase structure):
  R_major = 3/k

With the braid's natural k related to the mass gap (k ~ pi/L_arm ~ 0.2):
R_major ≈ 3/0.2 = 15 code units.

### Structure dimensions
- R_major (torus center to tube center): ~12-15 code units
- R_minor (tube radius, same as braid): 3 code units
- Total extent: 2*(R_major + R_minor) = ~36 code units
- The torus inner hole: R_major - R_minor = ~12 code units (well-defined hole)

### For a proton
One toroidal braid = one constituent. A proton could be:
- Three interlocking toroidal braids (each in a different plane)
- UUD chirality (winding direction around the torus)
- Held together by theta-mediated force between the tori
- Total extent: ~50-60 code units

Start with ONE toroidal braid first. If it survives, attempt multi-torus.

## Grid Requirements

| Configuration | Extent | L | N (min) | N (good) | Memory |
|---------------|--------|---|---------|----------|--------|
| Single torus | ~36 | 30 | 192 | 256 | 1.0/2.4 GB |
| Single torus (generous) | ~36 | 40 | 256 | 384 | 2.4/8.2 GB |
| Three interlocking tori | ~60 | 50 | 320 | 384 | 4.7/8.2 GB |

**Recommended**: Start with N=256, L=30 for single torus (2.4 GB, fits V100).
Scale to N=384, L=45 for multi-torus (8.2 GB, fits V100).

Resolution: at N=256, L=30: dx=0.235, tube radius = 12.8 cells. Adequate.

## Experiments

ALL on V100 GPU using CUDA. Log per-iteration time. Output SFA for evaluation.
Download SFA via rsync (supports resume + checksum).

### Phase 1: Single Toroidal Braid (1-2 hr V100)

Initialize ONE toroidal braid:
- R_major = 12, 15, 18 (sweep to find best)
- n_osc = 2, 3, 4 (oscillations around the ring)
- Tube carries all three fields with delta = {0, 3.0005, 4.4325}
- Traveling wave velocity: wave propagates AROUND the torus
- N=256, L=30, absorbing BC
- T=300

For each (R_major, n_osc) combination:
- Log: E_pot, P_int, theta_rms, cluster count every 5 time units
- Log: wall time per iteration (verify GPU is working)
- Output: SFA file (f32) for visualization
- Measure: survival time, binding retention, aspect ratio

**Success criterion**: E_pot retention > 30% at T=200, 1 cluster.

### Phase 2: Toroidal Braid Characterization (if Phase 1 succeeds, 1 hr)

For the best (R_major, n_osc) from Phase 1:
- Measure theta field around the torus: is it toroidal? Poloidal?
- Measure the "charge" (winding number) — is it well-defined?
- Compare with the linear braid's theta field
- Vary eta (curl coupling): is there a sweet spot?

### Phase 3: Multi-Torus Composite (if Phase 2 succeeds, 2-3 hr)

Three interlocking toroidal braids:
- Each in a different plane (xy, xz, yz)
- UUD chirality (proton analog) and UDD (neutron analog)
- N=384, L=45
- T=300 with absorbing BC
- Do the three tori attract/repel via theta? Do they form a bound state?

### Phase 4: Evolutionary Refinement (if Phase 3 shows promise, 2+ hr)

Use field_evolve.py (JAX, GPU) to optimize:
- Seed: 3-torus composite from Phase 3
- Fitness: E_pot retention at T=50 with absorbing BC
- Mutate: torus radii, crossing angles, phase relationships

## Tools

### Must build: `init_torus_braid.py`
Generate initial condition for a toroidal braid:
```python
# For each grid point (x,y,z):
#   1. Compute distance to torus center-circle: d_tube
#   2. Compute angle along torus: theta_tor = atan2(y, x)
#   3. Envelope: exp(-d_tube^2 / (2*R_minor^2))
#   4. Phase: k * R_major * theta_tor + delta[a]
#   5. Velocity: traveling wave around the torus
#
# For three strands (like braid3):
#   Each strand center spirals around the tube cross-section
#   strand_b at angle: n_twist * theta_tor + 2*pi*b/3
#   This gives a helical braid wound around the torus
```
Output: numpy .npz seed file, convertible to C binary via seed_to_c_init.py.

### Existing tools (from V37/V38)
| Tool | Status | Notes |
|------|--------|-------|
| `src/seedrun_cuda.cu` | Built | GPU sim, absorbing BC, SFA output |
| `src/braid_analyze_cuda.cu` | Built | Needs torus init mode added |
| `backprop/seed_to_c_init.py` | Built | numpy → C binary converter |
| `backprop/field_evolve.py` | Built | JAX evolutionary search |
| `scripts/deploy.sh` | Built | Needs update for new grids |
| `scripts/monitor.sh` | Built | 10s GPU monitoring |
| `v37/src/sfa_frag.c` | Built | Fragmentation analysis |
| `v37/src/sfa_structure.c` | Built | Structure analysis |
| `sfa/viewer/volview` | Built | Wireframe + depth-tested |

### Download protocol
Use rsync for all SFA transfers (supports resume + checksum):
```bash
rsync -avz --partial --progress \
  -e "ssh -o StrictHostKeyChecking=no -p PORT" \
  root@host:/path/file.sfa local/file.sfa
```

## Per-Iteration Timing Requirements

Every experiment must log:
- Wall time per Verlet step (ms)
- GPU utilization (from monitor.sh every 10s)
- Estimated time to completion

Expected V100 performance (from V36 benchmarks):
- N=128: ~1.9 ms/step → ~80 steps/s
- N=256: ~15 ms/step → ~67 steps/s (estimate, 8× more work)
- N=384: ~50 ms/step → ~20 steps/s (estimate)

For T=300 at N=256: ~300/(0.1*0.235) ≈ 12,766 steps × 15 ms = ~190s ≈ 3 min.
For T=300 at N=384: ~300/(0.1*0.157) ≈ 19,108 steps × 50 ms = ~955s ≈ 16 min.

## Budget
~6 hours V100 time, ~$0.80 total.
Must use GPU. Kill if idle.
