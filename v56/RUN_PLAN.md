# V56 — first SFA-producing run plan

## Target

Run a v56 stage-A simulation on the existing L=40 Voronoi mesh and
produce a cell-native SFA. The dynamics differ from v55 in three
specific ways; everything else (mesh, integrator, SFA format) is the
same.

## Kernel changes (foam_sim.c)

| Change | v55 | v56 stage A |
|--------|-----|-------------|
| Field count per cell | 6 (φ₀, φ₁, φ₂, θ₀, θ₁, θ₂) | 8 (relativistic-quaternion components M[0..7]) |
| Mass term | +m² φ + 0 × θ | g_{ii} m² M[i], g = diag(−,−,−,−,+,+,+,+) |
| Self-coupling | μ P²/(1+κP²), P = φ₀φ₁φ₂ | λ (\|M\|²_bulk − v²)² with \|M\|²_bulk = M[4]²+M[5]²+M[6]²+M[7]² − M[0]²−M[1]²−M[2]²−M[3]² |
| Curl coupling | η ∇×θ in φ EOM | none (multivector self-couples through bulk norm) |
| Init seed | braid (3-axis Gaussian carrier) | tunable; default Q-ball (localised perturbation in rotor part) |

Per-cell EOM:

```
□ M[i] = − λ (|M|²_bulk − v²) × g_{ii} × M[i]   (no sum over i)
```

Component-wise:

| i | basis | g_{ii} | EOM mass term sign |
|---|-------|--------|---------------------|
| 0 | e₄₁₀  | −1     | +λ(|M|²_bulk−v²)M[i]  (tachyonic — drives toward vacuum) |
| 1 | e₄₂₀  | −1     | +λ(...)M[i]                                            |
| 2 | e₄₃₀  | −1     | +λ(...)M[i]                                            |
| 3 | 𝟙     | −1     | +λ(...)M[i]                                            |
| 4 | e₁    | +1     | −λ(...)M[i]                                            |
| 5 | e₂    | +1     | −λ(...)M[i]                                            |
| 6 | e₃    | +1     | −λ(...)M[i]                                            |
| 7 | e₃₂₁  | +1     | −λ(...)M[i]                                            |

The negative `g_{ii}` on rotor components inverts the sign — for those
the same potential pushes them away from origin, generating the
spontaneous-symmetry-breaking vacuum at \|M\|² = v².

## Implementation steps

1. **Bump field count.** In `Sim` struct, change `phi[NFIELDS]`,
   `theta[NFIELDS]` (and their _vel, _acc) to `M[8]`, `M_vel[8]`,
   `M_acc[8]`. `NFIELDS=3` becomes `NCOMP=8`.
2. **Update operators.** `compute_grads_and_lap_all` already takes a
   `[6]` slot; bump to `[8]`. Same Laplacian, just one more component.
3. **Replace force computation.** `compute_forces` body becomes:
   ```c
   for (uint32_t i = 0; i < N; i++) {
       double Mb = 0;
       for (int k = 0; k < 4; k++) Mb -= s->M[k][i] * s->M[k][i];
       for (int k = 4; k < 8; k++) Mb += s->M[k][i] * s->M[k][i];
       double drive = lambda * (Mb - vsq);
       static const double g[8] = {-1,-1,-1,-1, +1,+1,+1,+1};
       for (int k = 0; k < 8; k++)
           s->M_acc[k][i] = lap[k][i] - drive * g[k] * s->M[k][i];
   }
   ```
4. **Update SFA columns.** `n_columns = 8` (was 6). Names: `M0..M7`
   with semantic codes (POSITION for spacelike, CUSTOM for the
   timelike rotor parts). Cell-native FCEL/FCEP machinery already
   handles arbitrary `n_columns`.
5. **New seeds.** Drop `init_braid`. Add:
   - `init_vacuum`: M[3]=v, all others 0. Should be stationary.
   - `init_qball`: M[3] = v + A·exp(−r²/2R²)·cos(ωt) at t=0,
     others 0. Localised perturbation, candidate Q-ball.
   - `init_skyrme_hedgehog`: see DERIVATION.md §7. Topological
     winding number 1.
6. **Diagnostics**. Replace P/phi_max with bulk-norm stats:
   `bulk_min`, `bulk_max`, `bulk_mean`, plus per-component RMS and
   max. Energy decomposition: kinetic, gradient, mass, quartic.
7. **No curl coupling**. Drop the `eta * curl(θ)` term from forces.
8. **Volview**. With 8 columns instead of 6, the existing renderer
   needs a new column→RGB mapping. Defer for the first run; we can
   render via cellsfa_to_voxel into a 6-column SFA by computing
   derived scalars (\|rotor\|, \|spatial\|, etc.).

## Run config (proposed)

`v56/foam/run_qball_L40.cfg`:

```
# V56 stage A — Q-ball candidate on existing L=40 mesh
# Field: 8-component multivector in PGA ℝ(3,1,1)
# Dynamics: Klein-Gordon + Higgs quartic
m = 1.5                  # bare mass scale (used in initial transient + diag)
lambda = 5.0             # quartic coupling strength
v = 0.5                  # vacuum expectation value sqrt(|M|²_bulk = v²)

# Q-ball seed parameters
A = 0.3                  # perturbation amplitude on top of vacuum
R_seed = 4.0             # localisation radius
omega_seed = 1.0         # oscillation frequency of the Q-ball test mode

T = 50
dt_factor = 0.025
init = qball

# Cell-native SFA output (FMSH + FCEL + FCEP + FMTL)
cell_native = 1
sfa_output = 0
cell_iframe_interval = 20
cell_model_interval = 5     # write FMTL every 5 I-frames (5x save vs embedded)
cell_delta_tol = 0.01
cell_omega = 1.0            # Q-ball mode frequency

output = qball_L40.sfa
diag_file = qball_L40_diag.tsv
snap_dt = 1.0               # 50 frames over T=50
diag_dt = 1.0
```

Mesh used: `/space/sfa/v55/meshes/foam_L40.bin` (3.26M cells, L=40).

## Expected wall time

Same per-step cost as v55 (8 components vs 6 ≈ 33% more work, but the
combined gradient pass already dominates). Estimate **~25-30 min wall**
for T=50 at 9.4× speedup baseline.

## What we look for

1. **Vacuum stability**: a baseline run with `init = vacuum` should
   show \|M\|² = v² to machine precision, no dynamics.
2. **Q-ball persistence**: the seeded perturbation either persists as
   a localised structure (Q-ball) or radiates away. Track the
   localisation radius vs time.
3. **Energy conservation**: the new potential is Lorentz-invariant by
   construction; drift should be small. Compare to v55's bounded drift.
4. **Visual character**: the rotor part `(M[0..3])` is what carries
   the spinor structure; the spatial part `(M[4..7])` is what carries
   the Higgs amplitude. Expect to see spatial cores with a rotor halo,
   roughly inverted from v55's φ-cores + θ-halo pattern.

## Stage-A vs stage-B decision point

Run the Q-ball seed for T=50. If:

- **The Q-ball persists with bounded drift**: stage A produces a
  particle-like state; proceed to richer experiments (two-Q-ball
  interaction, Lorentz boost test, etc.).
- **The Q-ball dissolves quickly**: stage-A Klein-Gordon dynamics
  don't admit non-topological solitons in this regime; jump to
  stage B (Dirac equation) or stage C (Skyrme topological term).
- **The Q-ball survives but drift is large**: discretization issue;
  re-derive the Hamiltonian explicitly and check operator consistency.

## Files this plan creates

```
v56/foam/
  foam_sim.c               (modified: 8-component, Higgs potential, drop curl)
  init.h                   (new: separate seed implementations)
  run_qball_L40.cfg        (run config)
  run_vacuum_L40.cfg       (sanity baseline)
v56/
  RUN_PLAN.md              (this doc)
  RESULTS.md               (filled in after the run)
```
