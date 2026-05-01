# V53 Grid Resolution Tests — N=192 Validation

## Summary

Three tests at N=192 (2x the original 96^3 grid) confirm that the template init
physical coordinate mapping works correctly, and produce identical behavior at both
resolutions. However, the results demonstrate that these objects are **breathing
oscillons, not stable particles**. The localized field structure periodically disperses
to near-zero and reforms, H_cross is not conserved, and no opposite-chirality
configuration has been found to be stable.

## Test Configuration

| Parameter | Test A | Test B | Test C |
|-----------|--------|--------|--------|
| N | 192 | 192 | 192 |
| L | 30 | 30.16 | 30 |
| dx | 0.314 | 0.316 (exact seed match) | 0.314 |
| Init | Template from 96^3 seed | Template from 96^3 seed | Fresh Lissajous at 192^3 |
| A_bg | 0 | 0 | 0 (after cooling) |
| BC | Periodic | Periodic | Absorbing T=300, then Periodic T=600 |
| T | 600 | 600 | 300 (cool) + 600 (run) |
| m^2 | 9.0 | 9.0 | 9.0 |

## Energy Conservation

| Test | E_initial | E_final | Drift |
|------|-----------|---------|-------|
| A | 5189 | 5182 | -0.14% |
| B | 5271 | 5263 | -0.15% |
| C cool | 68520 | 3390 | -95% (absorbing BC) |
| C run | 3380 | 3336 | -1.30% |

Energy is well-conserved under periodic BC. But energy conservation does not imply
stability — periodic BC recycles radiated energy back into the simulation volume,
masking dispersal.

## Critical Assessment: These Are Oscillons, Not Stable Particles

### Evidence against true stability

**1. Mass nearly vanishes during breathing cycle.**
Test A particle tracking shows total localized mass swinging from 183 down to 9
(a 95% drop) before reforming:

| t | total mass | E_pot | H_cross |
|---|-----------|-------|---------|
| 275 | 152.0 | -181.8 | -25.4 |
| 300 | 183.1 | -220.1 | -24.9 |
| 325 | 188.0 | -218.7 | -13.1 |
| 375 | 111.8 | -131.1 | -14.2 |
| 425 | 19.5 | -12.4 | -1.4 |
| 450 | **9.1** | **-4.6** | **-2.8** |
| 475 | 12.2 | -6.3 | -2.3 |
| 500 | 49.8 | -56.9 | -12.1 |
| 525 | 102.5 | -113.1 | -16.0 |

At t=450, E_pot=-4.6 means essentially no binding. The "particle" has dispersed
into low-amplitude waves filling the periodic box, then reconcentrates. A truly
stable soliton would maintain a persistent core.

**2. H_cross is not conserved — it tracks mass.**
H_cross drops from -25 to -2.8 when the mass disperses, and recovers when mass
reconcentrates. This is not a topological charge; it's simply proportional to the
amount of localized field present. A real chiral quantum number would remain constant
(or quantized) regardless of the field's spatial distribution.

**3. No stable opposite chirality exists.**
Multiple attempts to construct positive-H_cross particles (flip_theta, rev_phase seeds)
all failed under absorbing BC. The dynamics spontaneously drives any configuration
toward negative H_cross. If chirality were a genuine topological property, both signs
would be equally stable. The single-sign preference suggests it's a dynamical bias
from the equation structure, not a conserved charge.

**4. Periodic BC creates an illusion of stability.**
Under periodic BC, energy that radiates away wraps around and re-enters the
simulation volume. The oscillon "survives" because it's bathed in its own radiation,
which periodically reconcentrates. Under absorbing BC, these objects eventually die
(as seen in earlier V52/V53 experiments where particles under absorbing BC had
finite lifetimes).

**5. The theta growth is generic, not chiral.**
The eta*curl(phi) coupling means ANY localized phi oscillation will pump theta.
This is linear response from the coupling term, not evidence of a chiral mechanism.
The theta field grows monotonically regardless of the phi configuration's "chirality."

### What the |mu|/m^2 < 7 condition actually means

The derived stability condition is better understood as a **slow-dispersal condition**.
When alpha = |mu|/m^2 is small enough, the nonlinear potential V(P) produces
sub-threshold breathing: the oscillon's amplitude stays below the level where
radiation losses dominate. The oscillon disperses slowly rather than quickly. Under
periodic BC this appears as indefinite survival; under absorbing BC it has a long
but finite lifetime.

## Technical Results (Still Valid)

Despite the negative physics conclusion, the technical infrastructure works:

1. **Template init physical coordinate mapping**: Correctly maps seeds across
   different grid resolutions. Tests A (dx=0.314) and B (dx=0.316) produce
   identical behavior.

2. **Resolution independence**: N=192 reproduces N=96 dynamics exactly. The
   previous N=192 failure was caused by the index-based mapping bug, now fixed.

3. **Fresh seed generation**: Lissajous seeds generated and cooled at N=192
   produce the same breathing oscillons as template-transplanted seeds.

4. **CUDA kernel correctness**: Float atomicAdd, Chebyshev basis, physical
   coordinate mapping all verified at N=192.

## Particle Tracking Details

### Test A (template init)
- t=0-75: particle fragments into 2-3 clusters
- t=100+: settles to 2 clusters oscillating in anti-phase
- Total mass swings 9-183 with period ~50 time units
- All clusters have negative H_cross (proportional to mass, not conserved)

### Test C (fresh Lissajous seed)
- 3 persistent sub-clusters form from cooling
- Dominant cluster grows to mass 227 then oscillates
- Same non-conserved H_cross behavior

## Data Files

All data at `/space/scp/v53/grid_tests/`:
- `test_a.sfa` (6.5 GB, 26 frames) — template init, L=30
- `test_b.sfa` (6.5 GB, 26 frames) — template init, L=30.16
- `test_c_cool.sfa` (1.5 GB, 7 frames) — Lissajous seed cooling
- `test_c_run.sfa` (6.7 GB, 26 frames) — cooled seed periodic run
- `test_*_diag.tsv` — full diagnostic time series
- `test_*_particles.tsv` — per-particle tracking analysis
