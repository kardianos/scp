# Theta Field Characterization: Results

Three experiments testing whether the theta field around a Cosserat braid
behaves like electromagnetism. Data from 6-field Cosserat simulations
(N=80 grids, Verlet integrator, curl coupling with eta=0.5).

## Experiment 1: Radial Decay of Circular theta

**Setup**: Analyzed frame 200 (t=37.97) of `sfa_hires.sfa` (N=80, L=25, T=50,
eta=0.5, m_theta=0, 264 frames). Decomposed theta into cylindrical components
(theta_r, theta_phi, theta_z) around the braid axis. Averaged in shells dr=0.5.

**Key findings**:

1. **Azimuthal dominance is weak**: theta_phi^2 / theta_total^2 ~ 0.45-0.58.
   The azimuthal component exceeds radial by only 1.2-3.2x. This is not the
   overwhelming dominance expected for a magnetic field (which would be >95%
   azimuthal far from the wire).

2. **Radial decay is much slower than 1/r^2**: Power-law fit gives
   theta_phi^2 ~ r^(-0.50) (r=2-12 fit), far from the Biot-Savart prediction
   of r^(-2). Even a tighter fit (r=3-8) gives exponent 1.74, still not 2.0.
   The log-log residual is 0.33 -- the fit is poor because the profile is
   not a clean power law.

3. **Radial standing wave**: The signed theta_phi alternates sign with radius
   at any given time. At frame 200: negative at r=2, positive at r=4, mixed
   at r=6, positive at r=8, negative at r=10. Angular analysis confirms each
   shell is highly coherent in angle (|mean|/rms ~ 0.97), so the sign is real.

4. **Temporal oscillation**: Scanning frames 50-260, the signed theta_phi at
   every radius oscillates with period ~4 time units, flipping sign repeatedly.
   The time-averaged DC component is only 0.2% of the RMS oscillation amplitude.

**Verdict**: NOT CONFIRMED as Biot-Savart. The theta field response to the
braid's curl source is a propagating/standing cylindrical wave, not a static
1/r magnetic-like circulation. The oscillation dominates the DC component by
500:1.

## Experiment 2: Winding Reversal (Charge Conjugation)

**Setup**: Ran two T=30 simulations (N=80, L=25, eta=0.5, m_theta=0):
- W=+1: delta = (0, +3.0005, +4.4325)
- W=-1: delta = (0, -3.0005, -4.4325)

Compared signed theta_phi averaged over 22 frames each (t=10-30).

**Key findings**:

1. **Single-snapshot comparison is misleading**: At a single frame (t=23.73),
   the signed theta_phi has correlation +0.78 between windings (SAME sign
   pattern). This is because the oscillatory wave has the same phase at this
   snapshot.

2. **Time-averaged DC residual**: After averaging over 22 frames, the DC
   component (mean theta_phi) is tiny: |mean|/rms ~ 0.002 for both windings.
   However, the DC residuals between the two windings have correlation -0.68.
   The mean |sum| = 3.93e-5 vs mean |diff| = 1.46e-4 (cancellation ratio 0.48).

3. **Partial sign flip in DC**: For r > 3, the positive-winding DC is negative
   for r=3-9 and positive for r>10, while the negative-winding DC has the
   opposite pattern. But this signal is 0.2% of the oscillation.

4. **Energy is the same**: theta_phi^2 (energy density) is identical for both
   windings (ratio 0.978), confirming the same dynamics with different phase.

**Verdict**: WEAKLY CONFIRMED. The tiny DC component (~0.2% of oscillation)
does flip sign with winding (correlation -0.68). But the dominant theta
response is an oscillating wave whose direction is independent of winding.
Winding sets a small DC offset, not a strong static circulation.

## Experiment 3: Two-Braid Force Comparison

**Setup**: Three T=150 simulations (N=80, L=35, braids=2, D=15, diag=5):
1. **3-field baseline** (eta=0): pure phi gravity, no theta coupling
2. **6-field same winding** (eta=0.5): both braids W=+1
3. **6-field opposite winding** (eta=0.5): braid 1 W=+1, braid 2 W=-1

Tracked braid separation D(t) via energy-weighted centroid in left/right halves.

**Key findings**:

1. **theta ADDS attraction for same-winding braids**:
   - Mean D (late time): 3-field=9.22, 6f-same=7.76 (delta=-1.46)
   - Minimum D: 3-field=6.77, 6f-same=5.52 (delta=-1.25)
   - Early infall rate: 6f-same is 23% faster than 3-field
   - The braids fall together faster with theta coupling enabled.

2. **Opposite winding attracts less**:
   - Mean D (late time): 6f-opp=9.38 vs 6f-same=7.76 (delta=+1.62)
   - Early infall rate: 6f-opp is 36% SLOWER than 3-field
   - Opposite-winding braids barely fall inward at all at early times.

3. **Quantitative force hierarchy**: D_late: same < 3-field < opposite
   - Same-winding: theta adds ~1.5 units of attraction
   - Opposite-winding: theta adds ~0.2 units of REPULSION (or reduces attraction)

4. **Late-time dynamics are chaotic**: After the braids approach closely
   (D<8), the motion becomes oscillatory and hard to interpret. The early
   infall rate (t<50) is the cleanest measure.

**Verdict**: CONFIRMED. Same-winding braids attract 23% faster through theta
(parallel "currents" attract). Opposite-winding braids attract 36% slower
(anti-parallel "currents" repel, partially canceling the phi gravity).
This is the correct qualitative pattern for EM-like behavior.

## Overall Assessment

| Test | Expected for EM | Observed | Status |
|------|-----------------|----------|--------|
| theta_phi ~ 1/r | Clean 1/r power law | Standing wave, n=0.5 | NOT confirmed |
| Sign flip with winding | Full flip of DC field | 0.2% DC flip, 99.8% wave | WEAKLY confirmed |
| Same-winding force | Stronger attraction | 23% stronger infall | CONFIRMED |
| Opposite-winding force | Weaker attraction/repulsion | 36% slower infall | CONFIRMED |

**Summary**: The theta field is NOT a static magnetic field. It is an
oscillating wave field sourced by curl(phi). The dominant response is a
cylindrical standing wave, not a 1/r DC circulation.

However, the FORCE measurements clearly show that:
- Theta coupling adds attraction between same-winding braids
- Theta coupling reduces attraction between opposite-winding braids

This is qualitatively consistent with the EM analogy (parallel currents
attract, antiparallel repel), but the mechanism is wave-mediated rather
than static-field-mediated. The theta field acts more like a radiation
coupling (exchanging wave energy) than like a magnetostatic field.

The key distinction: in real EM, B ~ 1/r around a wire (static, long-range).
Here, theta_phi oscillates at every radius with a 0.2% DC bias. The force
effect is real but the field structure is fundamentally different from
classical electromagnetism.
