# V52 Immediate Next Steps

## CRITICAL: Re-run Diffuse Annealing with TRUE Periodic BC

The kernel audit (v52/kernel_audit.md) found that ALL our "periodic BC"
runs had the absorbing boundary layer active (default damp_width=3.0,
damp_rate=0.01). This caused 99.99% energy loss even with bc_type=2.

**Fix**: Add `damp_width=0` and `damp_rate=0` to the config.

The diffuse annealing run showed genuine condensation at t=50 (sphere
forming in volview) before the damping killed it. Re-running with true
periodic BC should let this condensation proceed to completion.

### Run Config
```
N = 384
L = 50
T = 2000
dt_factor = 0.025
m = 1.5
m_theta = 0
eta = 0.5
mu = -41.345
kappa = 50.0
alpha_cs = 0.1
beta_h = 0.5
bc_type = 2
damp_width = 0
damp_rate = 0
init = sfa
init_sfa = diffuse_seed.sfa
A_bg = 0.1
output = diffuse_true_periodic.sfa
diag_file = diffuse_true_periodic_diag.tsv
precision = f16
snap_dt = 50.0
diag_dt = 5.0
```

### Expected Behavior
- Energy should be conserved to O(dt²) — no secular drift
- The condensation at t=50 should continue rather than dissipate
- If V(P) binding captures the condensing structure, a stable
  particle should emerge
- phi_max should grow and stabilize (not decay to 0.005)

### Also Re-run
- Spherical seed (gen_spherical_seed) with damp_width=0
- Hedgehog seed with damp_width=0
- All previous "periodic BC" tests were actually absorbing — results invalid

## Pending: PDF Conversion Agent
Agent converting ref/chapter*.pdf to markdown — check status.

## Pending: Vectorizer Volview Integration
The polynomial vectorizer works (sfa_vectorize.c). Next:
- 3D tensor-product patches instead of 1D slices
- Write vectors as SFA overlay columns for volview
- Iterative error minimization loop
