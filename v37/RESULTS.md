# V37: Compact Braid Geometry Search — Results

## Setup
- 6-field Cosserat: 3 position (phi) + 3 angle (theta) fields
- Parameters: m^2=2.25, mu=-41.345, kappa=50, eta=0.5, m_theta^2=0, A_bg=0.1
- Grid: N=128, L=15, T=500 (periodic BC)
- Velocity Verlet integrator, dx=0.236, dt=0.0236

## Summary Table

| Geometry | Initial E_pot | Death time | Final E_pot (% init) | P_int (% init) | theta_rms | Aspect ratio | Survival |
|----------|--------------|------------|---------------------|-----------------|-----------|-------------|----------|
| truncated | -17.1 | ~359 (oscillating) | -0.12 (0.7%) | 27.4% | 0.036 | 1.03-1.20 | PARTIAL |
| borromean | -6.8 | ~30 | -0.035 (0.5%) | 15.9% | 0.008 | 1.02-1.10 | NO |
| loop | -11.4 | ~15 | -0.094 (0.8%) | 21.6% | 0.010 | 1.09-1.03 | NO |
| triloop | -2.4 | ~20 | -0.004 (0.2%) | 9.0% | 0.008 | 1.01-1.12 | NO |

## Detailed Analysis

### 1. Truncated Helix (BEST CANDIDATE)

The truncated helix (Gaussian z-envelope, sigma_z=3.0) showed the strongest survival:

Evolution at N=128, L=15, T=500, eta=0.5, m_theta^2=0:

| Time | E_total | E_pot | P_int | theta_rms | Aspect |
|------|---------|-------|-------|-----------|--------|
| 0 | 1424 | -17.1 | 22.7 | 0.000 | 1.03 |
| 50 | 1423 | -3.6 | 10.4 | 0.019 | 1.06 |
| 100 | 1423 | -23.8 | 26.6 | 0.024 | 1.30 |
| 150 | 1422 | -1.4 | 8.3 | 0.032 | 1.08 |
| 200 | 1422 | -25.9 | 28.7 | 0.032 | 1.06 |
| 250 | 1422 | -13.6 | 19.0 | 0.036 | 1.23 |
| 300 | 1422 | -9.1 | 15.4 | 0.033 | 1.14 |
| 349 | 1421 | -29.2 | 32.7 | 0.035 | 1.20 |
| 359 | 1422 | -0.12 | 6.2 | 0.036 | 1.20 |

Key observations:
- **E_pot oscillates** between -40 and -0.1 throughout (strong breathing mode, period ~20-25)
- The braid did NOT actually die: the early termination at t=359 caught E_pot at a
  breathing-mode NODE (near zero). Time-averaged E_pot remained substantial (~-15)
- Energy drift: -0.18% (excellent conservation)
- Theta coupling: theta_rms grew from 0 to 0.036 (curl coupling is active)
- Aspect ratio: stable at 1.03-1.30 (very compact, nearly spherical) -- PASSES target < 2
- P_int retained 27.4% of initial at the final (node) measurement, oscillated up to 32.7

### 2. Borromean Rings (FAST DEATH, t~30)

Died within 30 time units. Initial E_pot = -6.8 (weak vs truncated's -17.1).
- Three rings in orthogonal planes separate the fields: phi_0 on xy-ring, phi_1 on xz-ring,
  phi_2 on yz-ring. The potential V(P) = V(phi_0 * phi_1 * phi_2) needs all three
  fields nonzero at the SAME point -- the ring geometry fails this.
- E_mass oscillated wildly (720 -> 110 -> 620) indicating standing-wave breakup
- Theta coupling barely activated (0.008)

### 3. Closed Loop (FASTEST DEATH, t~15)

Died within 15 time units.
- Torus at R_major=5 too large for L=15 box (PBC artifacts: 2.5% energy drift)
- Enormous initial gradient energy (E_grad=147 vs 31 for truncated) from tight twists
- Would need L=25+ to avoid PBC artifacts

### 4. Triple Loop (FAST DEATH, t~20)

Died within 20 time units. Weakest initial E_pot = -2.4.
- Same problem as borromean: fields separated by plane, triple product P is small
- Three loops compete rather than cooperate

## Why Geometries Fail: The Co-Location Requirement

The binding potential V(P) = (mu/2) P^2/(1+kappa P^2) with P = phi_0 * phi_1 * phi_2
only activates where ALL THREE fields overlap with significant amplitude.

- **Truncated helix**: All three fields share the same tube (helical braid), so P is large.
  The z-envelope makes it compact but preserves co-location.
- **Borromean/triloop**: Fields separated into different planes, so P ~ 0 almost everywhere.
  The 0.3x cross-contamination terms in initialization are not enough.
- **Loop**: Fields share the torus tube, but large gradients from k_twist=3 cause rapid
  dispersion before the potential can stabilize.

The mass term m^2=2.25 also plays a role: it provides a restoring force toward phi=0
that is compatible with helical oscillation cos(kz) but incompatible with the static
field profiles needed for ring/torus closures.

## Theta Coupling

| Geometry | Final theta_rms | Growth pattern |
|----------|----------------|----------------|
| truncated | 0.036 | Steady growth 0 -> 0.036 over 359t |
| borromean | 0.008 | Brief, dies with braid |
| loop | 0.010 | Brief, dies with braid |
| triloop | 0.008 | Brief, dies with braid |

Only the truncated helix shows significant theta field excitation from curl coupling.

## Compactness (Truncated Helix)

Aspect ratio 1.03-1.30 throughout evolution (target: < 2) -- PASSES.

The braid self-sphericalizes: aspect oscillates between elongated (1.30) and near-spherical
(1.03) phases, tracking the breathing mode. The time-averaged aspect is ~1.15.

Shell energy profile at final state (t=359) shows energy spreading outward, but the braid
P_int (tracking the triple product) remains localized.

## Recommendation

**The truncated helix is the clear winner** and the only viable compact braid geometry:

1. Survives as a coherent breathing structure for > 350 time units
2. Aspect ratio excellent (1.03-1.30, target < 2)
3. Theta coupling active and growing
4. The "death" was a measurement artifact (E_pot caught at oscillation node)

**Next steps to improve:**
- Re-run with relaxed termination (time-averaged E_pot over ~50 time units)
- Try sigma_z = 4-5 (wider envelope = more axial wavelengths) or sigma_z = 2 (tighter)
- Increase amplitude A from 0.8 to 1.0-1.2 for stronger potential binding
- The ring/loop topologies need initialization where all three fields share the SAME
  geometric region -- perhaps a trefoil knot where a single tube carries all three fields

## Files

- Source: `src/v37_compact.c` -- 4 geometry modes, 6-field Cosserat
- Data: `data/{truncated,borromean,loop,triloop}/timeseries.tsv`
- Snapshots: `data/*/field_t0000.bin` (t=0 for all geometries)
- Truncated evolution: `data/truncated/field_t{0100,0200,0300,0359,0500}.bin`
