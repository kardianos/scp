# V55 — Voronoi Foam Discretization (exploratory)

## Motivation

V54 concluded that the 6-field Cosserat equations under voxel discretization
produce no stable particles across 750+ configurations. Before declaring the
equations dead, test whether the discretization is the problem.

Two specific hypotheses to test in order:

1. **Voronoi foam mesh** (cell-based, density emerges from cell volume,
   isotropic stencils). Fixes voxel grid's cubic anisotropy and lets cells
   adapt to soliton shape.

2. **Integer field amplitudes** (quantum-like discrete representation). The
   v53 stability condition α<7 is about whether breathing frequency exceeds
   the radiation gap. If φ is integer-quantized, there is a natural floor
   below which oscillations are forbidden — the v54-unstable α=18.4 regime
   may stabilize because tiny continuous radiation leakage is structurally
   absent.

Both are NEW kernels (not modifications of `sfa/sim/scp_sim.*`). The current
unified kernel is preserved per CLAUDE.md kernel policy.

## Stage 0 — V34 Recreation (this stage)

Goal: reproduce the v34 image (`v34/torsion_coupling/screenshot.png`) — a
single braid with surviving structure and clear theta halo, under the
original V34 6-field Cosserat (no C4 hardening).

This validates that:
- The current kernel (with `alpha_cs=0`, `beta_h=0`) reproduces V34 physics
- The single-braid seed at A=0.8, R_tube=3.0, ellip=0.3325 still works
- We have a known-good baseline before changing the discretization

Original parameters (from `v34/torsion_coupling/RESULTS_cosserat.md`,
m_theta=0 stable case):
- N=128, L=20, T=300, dt = 0.10 × dx
- m²=2.25, m_θ²=0 (massless), η=0.5, μ=-41.345, κ=50
- Braid: A=0.8, R_tube=3.0, ellip=0.3325, δ={0, 3.0005, 4.4325}
- A_bg=0.1, periodic BC, no damping
- Result: braid survives T=300, θ_rms=0.062, E_pot=-44, drift=-0.6%

Success criterion: a snapshot at t≈300 visually matching the screenshot
(two stacked bright P-cores along z, blue θ halo, persistent structure).

## Stage 1 — Voronoi Foam Mesh

Conditional on Stage 0 working. Build a CPU prototype that runs the same
6-field Cosserat on an unstructured Voronoi cell complex:
- Each cell has position (centroid), volume V_c, neighbor list
- Field values stored at cell centroids; gradients/curl via Delaunay-dual
  finite-volume stencils
- Density implicit: ρ_local ~ N_quanta_per_cell / V_cell (cells with high
  field activity are smaller)
- Initial mesh: Poisson-disk sample, refined where |φ| or |∇φ| is large
- Test: same single-braid seed, see if the braid survives T=300 with
  comparable energy/drift behavior

Compare to Stage 0 voxel result. The interesting outcomes:
- (A) Voronoi result matches voxel: discretization not the issue → return to v54 conclusion
- (B) Voronoi result better: try hard cases (v54 m=1.5 α=18.4) under foam
- (C) Voronoi result worse: foam adds error not symmetry → abandon

## Stage 2 — Integer Quantum Representation

Conditional on Stage 1 producing a working foam kernel (any of A/B/C
where cells work mechanically). Replace double field amplitudes with
integers, fixed-point arithmetic for V(P)=μP/(1+κP²), discrete-step
velocity Verlet.

Test: a v54-unstable case (m=1.5, μ=-41.345, α=18.4) with integer
amplitudes. Hypothesis: radiation gap that floats over continuously in
v54 becomes a hard step at amplitude=1, suppressing slow leakage.

## Files

- `recreate_v34.cfg` — Stage 0 config
- `v34_screenshot_target.png` — for visual reference (copied from v34/)
