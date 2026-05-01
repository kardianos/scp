# V55 Stage 1a Results — Voronoi Foam Cosserat (Eulerian)

## Status (cell-native + 9.4× CPU speedup): production-ready prototype

The kernel now writes its own native cell-data SFA format (FMSH for the
mesh, FCEL for per-cell field values) — no voxel resampling at simulation
time, file ~50% smaller than the voxel SFA, and a separate
`cellsfa_to_voxel` converter feeds the existing volview pipeline.

CPU optimisations (`v55/foam/bench/results.tsv`) compounded to **9.4×**:

| Stage | ms/step | speedup | drift |
|-------|---------|---------|-------|
| baseline | 375.1 | 1.00× | -0.446% |
| +face AoS | 221.0 | 1.70× | -0.446% |
| +Morton reorder | 103.7 | 3.62× | -0.446% |
| +combined grad pass | 39.9 | **9.40×** | -0.446% |

L=40 (3.26M cells, voxel-output mode) wall time dropped from a projected
4 hours to 21 minutes.

### Cell-native SFA pipeline (FMSH + FCEL)

`sfa.h` extension (`SFA_CHUNK_FMSH`, `SFA_CHUNK_FCEL`,
`SFA_FRAME_MESH`/`CELL` codes; `sfa_write_mesh_frame`,
`sfa_write_cell_frame`, `sfa_read_mesh_frame`, `sfa_read_cell_frame`).

The mesh is stored as a frame, not a header chunk: subsequent mesh frames
supersede earlier ones, so a single SFA can hold high-res keyframes
followed by lower-res cell frames, or remeshed adaptive sequences. The
JMPF index distinguishes voxel/vec/mesh/cell frame types via the
`frame_type` field; readers without FMSH/FCEL knowledge skip them via the
JMPF jump table without breaking.

Pipeline now:

```
foam_sim (cell_native=1) → cell-native SFA  (FMSH + FCEL frames)
                                ↓
                 cellsfa_to_voxel  (one-time conversion)
                                ↓
                   voxel SFA  (FRMD frames, volview-compatible)
```

Smoke test (50k mesh, T=10):

| Metric | Voxel SFA (existing path) | Cell-native SFA |
|--------|---------------------------|-----------------|
| Size | 15.5 MB | 7.2 MB (**46%**) |
| Wall time | 5 sec + voxel index | 1 sec |
| Round-trip → voxel render | direct | via cellsfa_to_voxel |
| Drift at T=10 | -0.379% | -0.379% |

L=40 production run (3.26M cells, T=50, voxel_N=192 rendering):

| Metric | Voxel SFA | Cell-native SFA |
|--------|-----------|-----------------|
| Size | 3.36 GB | 1.89 GB (**56%**) |
| Wall time | 21 min | 22 min |
| Drift envelope | ±0.9% | ±0.9% (identical) |
| Final drift | +0.048% | +0.048% (identical) |
| Frames | 27 (FRMD) | 27 (1 FMSH + 26 FCEL) |
| Visual at v34 cam | reference | indistinguishable |

The cell-native size advantage scales with frame count: at fewer frames
the mesh overhead dominates and savings are smaller; at many frames the
mesh amortises and savings approach the per-frame ratio (~25 MB FCEL vs
129 MB FRMD ≈ 5× reduction).

### volview cell-native support

`sfa/volview/main.go` extended (~340 LOC additions): chunk-type constants,
auto-detection of cell-native files via JMPF scan, FMSH parser, FCEL
parser with on-the-fly cell→voxel resampling. The mesh is read once per
`FMSH` frame seen; subsequent FCEL frames reuse the cached voxel-to-cell
index. CLI flag `-cell-voxel-N` (default 192) controls the rendering
grid resolution independently of the mesh's cell count.

Installed at `$HOME/bin/volview`. Reads cell-native and voxel SFAs
transparently — passes through unchanged for FRMD/FRVD frames, decodes
via the cell-native path when FMSH/FCEL/FCEP are detected.

### Dense temporal capture: FCEL v2 + FCEP P-frames

For dense `snap_dt` (≤ 1 t-unit) the kernel now writes:

- **FCEL v2 keyframes** with embedded temporal model
  `value(t) = mean[c, k] + amp[c, k] · cos(ω·t + phase[c, k])` per cell,
  per column. Refit at every I-frame using sliding accumulator sums
  (Fourier projection over the previous interval).
- **FCEP P-frames** carrying only `(cell_id, residual)` pairs for cells
  whose actual value differs from the model prediction by more than
  `cell_delta_tol`. Sparse encoding when the model fits.

Config knobs (in `foam_sim` cfg):

```
cell_iframe_interval = 20    # I-frame every N snaps; should cover ≥ 1 period
cell_delta_tol       = 0.02  # absolute threshold; ~20% of A_bg keeps vacuum cells
cell_omega           = 1.508 # carrier-wave frequency, default for L=20 m=1.5
```

Reader algorithm (cellsfa_to_voxel, volview):
- On FCEL v2 → store cell values + load temporal model
- On FCEP → reconstruct `state[c, k] = mean[c, k] + amp[c, k]·cos(ω·t + phase[c, k]) + delta[c, k]`
  where `delta[c, k] = 0` for cells not in change list

Smoke test (50k mesh, T=10, snap_dt=0.25 → 41 snaps, tol=0.02):

| Phase | n_changed/N_cells | Comment |
|-------|-------------------|---------|
| Bootstrap (no fit) | 100% | First I-frame interval, model is constant |
| After 1st refit | 28.5% | Model captures vacuum cells; soliton + boundary stay |
| After 2nd refit | (tracked) | Model fluctuates with soliton breathing phase |

File size 46 MB at tol=0.02 vs 73 MB at the strict 0.005 — the model is
real but not perfect, mostly because the soliton interior has both
carrier (ω≈1.5) and breathing (ω≈2.86) frequencies and the single-ω
model can only capture one. For pure vacuum cells far from the braid
the residuals are <0.001.

Round-trip visualisation (rendered at v34 camera through `$HOME/bin/volview`)
matches the cell-native and voxel paths — the braid cores, θ halo, and
green curl bands are all preserved.

Tuning recommendations:

| Use case | iframe_interval | delta_tol | Expected sparsity |
|----------|-----------------|-----------|-------------------|
| Visualisation | 20–50 (≥ 1 period) | 0.02–0.05 | 30–60% per P-frame |
| Quality preview | 10 | 0.005 | 80–100% (more I-frames, less compression but tight) |
| Long-run capture | 50 | 0.05 | 5–30% (mostly soliton interior) |

## Status (initial post-fix): **PROTOTYPE VALIDATED — visually indistinguishable from voxel**

After fixing a missing 1/2 factor in the divergence-theorem gradient (the
discrete curl was 2× too big, pumping fake energy into θ via the curl
coupling) and switching to face-summed E_grad/E_coupling diagnostics, a
500k-cell foam run is *visually identical* to the voxel reference at the
same time points (T=10 and T=50, identical camera). Energy drift is
bounded around -0.8% with no secular trend.

| Metric | Voxel (128³ = 2.1M) | Foam initial | Foam post-fix (500k) |
|--------|---------------------|--------------|---------------------|
| Drift at t=50 | -0.11% | +5.20% (growing) | -1.08% (oscillating) |
| θ_rms at t=10 | 0.031 | 0.070 (2× too fast) | 0.031 (matches) |
| θ_rms at t=50 | 0.035 | 0.054 | 0.034 (matches) |
| phi_max range | 0.7-0.8 | 0.7-1.2 | 0.7-1.2 |
| Visual at v34 cam | reference | qualitatively similar | indistinguishable |

The fix isolates the issue: it was discrete-curl magnitude, not topology
of the unstructured mesh. With the correct operator, Voronoi geometry
gives the *same* physics as voxels for this seed and parameters.

## Status (initial): PROTOTYPE WORKS — Same physics, with caveats

The foam kernel runs the 6-field Cosserat dynamics on an unstructured Voronoi
mesh and produces braid behavior consistent with the voxel kernel:
z-aligned braid forms, θ excites via curl coupling, V(P) potential traps
the field into bright cores, periodic BC respected. Visual character at
t=50 matches v34's voxel run.

Key caveats: (1) ~5% energy drift over T=50 from non-Hamiltonian
discretization of the curl coupling, (2) θ excitation is ~2× faster than
voxel (related to (1)), (3) 45k cells is 47× lower resolution than the
2M-voxel reference run.

## Build artefacts

```
v55/foam/voro_src/         BSD-licensed voro++ source vendored from
                           github.com/chr1shr/voro
v55/foam/voro_src/libvoro++.a   562 KB static library
v55/foam/gen_foam_mesh.cpp      C++ mesh generator using voro++
v55/foam/foam_sim.c             C dynamics kernel (OpenMP)
v55/foam/foam_to_voxel.c        Resampler: per-cell .fsnp → voxel SFA
v55/foam/foam_mesh.bin formats   binary mesh (header + cells + faces + CSR)
v55/foam/foam_*.fsnp             per-cell snapshot binary (custom format)
```

Build:
```
cd voro_src && make libvoro++.a    # one-time
g++ -O3 -Wall -Ivoro_src -o gen_foam_mesh gen_foam_mesh.cpp voro_src/libvoro++.a
gcc -O3 -march=native -fopenmp -o foam_sim foam_sim.c -lm
gcc -O3 -fopenmp -o foam_to_voxel foam_to_voxel.c -lzstd -lm
```

## Stage 1a runs

### Mesh: 45235 cells in [-20, 20]^3, periodic

Built via `gen_foam_mesh -L 20 -N 50000 -o foam50k.bin` in 1 second.

```
N_cells = 45235  (target 50000; Bridson Poisson-disk tolerance)
N_faces = 339473  (each face stored once)
avg degree = 15.0  (canonical for 3D Voronoi)
V_cell mean = 1.415, range [0.876, 2.406]  (max/min ratio 2.7×)
V_total = 64000.000 (matches box exactly — periodic Voronoi closes)
dx_min = 0.941  (smallest face-pair distance)
```

### Run: T=50, dt_factor=0.025, dt=0.0235, 850 steps, 29s wall

| t | E_total | drift | E_pot | phi_max | θ_rms |
|---|---------|-------|-------|---------|-------|
| 0    | 5561 | 0.00% | -120.7 | 0.89 | 0.000 |
| 9.88 | 5725 | +2.94% | -178.5 | 1.18 | 0.070 |
| 19.75 | 5799 | +4.28% | -143.4 | 0.87 | 0.065 |
| 29.63 | 5850 | +5.20% | -204.1 | 0.74 | 0.054 |
| 39.50 | 5849 | +5.17% | -179.7 | 0.75 | 0.078 |
| 49.38 | 5755 | +3.48% | -98.1 | 0.72 | 0.055 |

E_pot oscillates -98 to -204 (voxel range -127 to -230, similar).
phi_max settles 0.7-0.8 (matches voxel).
θ_rms reaches 0.066 by t=10 (voxel reaches 0.062 only by t=300 — 30× faster).

### Voxel reference (v34_full at t=50)

```
E drift = -0.107%, E_pot = -169, phi_max = 0.78, θ_rms = 0.035
```

## What works

1. **Mesh generation** is robust. voro++ produces a clean periodic
   Voronoi tessellation in seconds for 50k cells. Cell volumes sum to
   the box volume to machine precision; no unbounded cells; average
   degree 15 is canonical.

2. **Operators are stable.** The two-point flux Laplacian and the
   divergence-theorem gradient produce a well-behaved time evolution
   with bounded drift (no exponential blowup, no NaNs).

3. **Physics qualitatively matches voxel.** Same z-aligned braid, same
   θ excitation channel, same V(P) binding cores, same breathing
   period in E_pot oscillations. Visual rendering (after nearest-
   neighbor resample to 128^3 voxels) shows the same structural
   character as the voxel reference.

4. **Periodic BC works.** Face deltas already encode the periodic wrap
   so the kernel is BC-agnostic. Energy across the periodic boundaries
   evolves consistently.

## What needs fixing

### Issue 1: Curl coupling is not strictly Hamiltonian

The continuum identity ∫ φ · curl(θ) dV = ∫ θ · curl(φ) dV (on a
periodic manifold) does not hold exactly for the discrete operators.
With:
  curl(F)_a(c) = ε_abc × (1/V_c) Σ_n A_cn × n̂_cn,b × (F_c(n) − F_c(c))

the two integrals differ by O(h²) per face. Verlet then conserves a
slightly-different discrete Hamiltonian, so the cell-based diagnostic
shows ~5% drift over T=50. The drift is bounded (oscillates rather
than grows monotonically) but it artificially pumps energy into the
θ field, which is why θ_rms saturates 30× faster than in voxel.

**Fix**: replace the cell-centered curl with a face-centered formulation
where the same face contribution appears in both φ and θ equations of
motion. Specifically, define a face quantity:
  Q_f = A_f × ε_abc × n̂_b × (φ_a(c2) + φ_a(c1))/2 × (θ_c(c2) − θ_c(c1))
  H_coupling = -η Σ_f Q_f
The φ and θ EOMs are then both gradients of the same H, and Verlet
becomes exactly symplectic.

### Issue 2: Energy diagnostic uses cell-based |∇f|² × V_c

The discrete conserved energy on a Voronoi mesh is face-based:
  E_grad = (1/2) Σ_f A_f × (φ(b) − φ(a))² / d_f
not the cell-summed |∇f|²V_c that I computed. Fixing the diagnostic
will show a smaller "drift" (closer to the genuine integrator error,
not the discretization error).

### Issue 3: Resolution

45k cells at L=20 means cell radius ~0.7 vs braid R_tube=3.0 — only
~4 cells across the braid. For a fair voxel-vs-foam comparison we
need ~250k cells (cell radius ~0.3, matching voxel dx=0.315). Mesh
build for 250k cells should still be under a minute.

## Side-by-side images

```
v55/snapshots/voxel_t50.webp      voxel (v34_full) at t=50
v55/snapshots/foam/frame_*.webp   foam (foam_50k) at t=0,10,20,30,40,50
```

Both rendered at the v34 camera (elevation=0.8, azimuth=0.5, distance=2.0,
brightness=2.5, opacity=4.0).

## Next: Stage 1b candidates

In rough order of usefulness for the hypothesis (does discretization
matter for v54's stability puzzle?):

1. **Fix the symmetric curl coupling** — replace the cell-centered curl
   with a face-centered Hamiltonian formulation. Estimated 200 LOC.
   This should drop drift to <0.5% over T=50 and make the θ excitation
   rate match voxel.

2. **Run a v54-unstable parameter** (m=1.5, μ=-41.345, α=18.4) on the
   foam mesh at 250k resolution. Does it dissolve like voxel, or
   survive due to isotropic stencils?

3. **Stage 2: integer field amplitudes.** Replace `double phi[N]` with
   `int32_t phi[N]` plus a fixed-point V(P), test whether the
   quantization floor stabilizes the v54 instability.
