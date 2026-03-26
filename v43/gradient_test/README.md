# V43 Gradient Force Test

Test whether a compact proton (UUD composite) responds to a density gradient.
Control: single z-aligned braid in the same gradient.

## Experiment Design

### Test Object A: UUD Proton (from averaged template)
- Pre-converged proton from v43/proton_formation
- Near-spherical, isotropic coupling to background
- The physical particle

### Test Object B: Single Braid (control)
- z-aligned braid from gen_braid
- Anisotropic coupling (z-preferred)
- The V33 test object

### Grid and Gradient
- N=512, L=100 (large domain, good separation from boundaries)
- bc_type=1 (gradient pinned)
- gradient_A_high=0.15, gradient_A_low=0.05
- Gradient margin=5 (5 pinned slabs each side)
- Object placed at center (0,0,0)

### Runs (8 total)
For each test object (proton, braid):
- eta=0.0 (3-field gravity only): T=200
- eta=0.5 (6-field gravity+EM): T=200

For each eta, 2 gradient strengths:
- Gentle: A_high=0.12, A_low=0.08 (50% variation)
- Steep: A_high=0.15, A_low=0.05 (200% variation)

### Measurements
- Track |P|-weighted centroid every snap_dt=1.0
- Fit x_cm(t) for drift velocity and acceleration
- Compare proton vs braid response
- Compare eta=0 vs eta=0.5

### Expected GPU time
- N=512 on V100: ~1.9 ms/step, dt≈0.01, T=200 → 20000 steps → ~38 sec per run
- 8 runs → ~5 minutes total compute
- Add seed generation and download time: ~30 minutes total instance time
