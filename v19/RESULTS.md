# V19 Results: Three-Mode Resonance on S³ + Skyrme

## Summary

**Phase 0 (1D radial)**: Marginal oscillon signature. Non-topological bump disperses
rapidly; dominant frequency omega* = 14.2 with quasi-periodic ringing but total energy
drops by 98% over 500 time units.

**Phase 1 (3D triad)**: **Definitive null result.** No B=0 oscillon, no triad phase
locking, no persistent resonance. All non-topological configurations disperse at the
same rate as the free wave equation, regardless of:
- Skyrme coupling (c4=0 vs c4=2)
- Number of active modes (1 vs 3)
- Initialization method (bump, kick, displaced)
- Amplitude (A=0.1 to A=1.0)

**Phase 2 (density coupling)**: Not attempted — precondition (persistent resonance) not met.

---

## Phase 0: 1D Radial Frequency Calibration

**Code**: `src/radial1d.c`
**Parameters**: A=0.8, sigma=1.5, c4=2 (pure sigma model dynamics), Nr=2000, rmax=30

| Quantity | Value |
|----------|-------|
| Dominant spectral peak | omega* = 14.21 |
| Secondary peaks | omega = 4.16, 0.88 |
| E_kin contrast (oscillon test) | 0.997 (marginal) |
| Energy at t=500 | ~2% of initial (98% radiated) |

The 1D bump rings quasi-periodically but loses nearly all energy to radiation. The
"oscillon detection" is based on E_kin contrast but is misleading — the configuration
is not trapped, just reverberating in a dispersing shell.

---

## Phase 1: 3D Triad Experiments

**Code**: `src/triad3d.c`
**Grid**: N=100, L=16, dx=0.162, dt=0.023, c4=2.0

### Dispersal Timescale Comparison

Dispersal time T_disp defined as time when f_core drops below 0.05:

| Test | Description | c4 | Modes | A | T_disp | Relative |
|------|-------------|-----|-------|---|--------|----------|
| 1A   | Hedgehog bump | 2 | 3 (sym) | 0.3 | 4.2 | 1.0x (baseline) |
| 2    | Kicked triad | 2 | 3 (triad) | 0.3 | 5.8 | 1.4x |
| 2'   | Kicked triad (high-A) | 2 | 3 (triad) | 1.0 | 5.6 | 1.3x |
| 3    | Displaced triad | 2 | 3 (triad) | 0.2/0.2/0.28 | 5.7 | 1.4x |
| 4    | Single mode | 2 | 1 | 0.28 | 5.3 | 1.3x |
| 5    | Bump (no Skyrme) | 0 | 3 (sym) | 0.3 | 4.2 | 1.0x |
| 6    | Triad (no Skyrme) | 0 | 3 (triad) | 0.3 | 5.7 | 1.4x |

**Key finding**: The triad is ~1.4x slower to disperse than the symmetric bump, but this
is entirely explained by the different spatial structure (sinc envelope vs Gaussian, larger
sigma=2.0 vs 1.5). The single-mode control (Test 4) disperses at 1.3x — nearly identical
to the triad. **Three modes provide no advantage over one mode.**

The Skyrme term (c4=0 vs c4=2) makes **no measurable difference** to dispersal time for
B=0 configurations.

### Phase-Lock Diagnostic

For all triad tests (2, 2', 3, 6), the phase diagnostic Phi(t) = phi_3 - phi_1 - phi_2
wanders randomly across [-pi, +pi] with no evidence of locking:

| Test | Phi behavior | Verdict |
|------|-------------|---------|
| 2 (A=0.3) | Random walk, range [-2.7, +2.2] | No lock |
| 2' (A=1.0) | Random walk, range [-2.7, +2.0] | No lock |
| 3 | Random walk, range [-2.5, +2.8] | No lock |
| 6 (c4=0) | Random walk, range [-2.7, +2.6] | No lock |

### Per-Mode Energy Exchange

In all triad tests, mode energies E1, E2, E3 decrease **monotonically together**.
No periodic energy sloshing was observed at any amplitude or coupling strength.

| Test | E1/E2/E3 behavior | Energy exchange? |
|------|-------------------|-----------------|
| 2 (A=0.3) | Equal, monotonic decrease | No |
| 2' (A=1.0) | Slight asymmetry (E3 > E1,E2 by ~10%), monotonic decrease | No |
| 3 | Persistent asymmetry (E3 > E1 > E2), monotonic decrease | No |

The slight asymmetry in Test 2' and 3 reflects the different initial energies of the
three modes (different radial profiles g_1, g_2, g_3) — not nonlinear coupling.

### Frequency Spectrum

The DFT of q_1(origin, t) for Test 2 shows a **flat noise spectrum** with no discernible
peaks. No coherent oscillation frequency is established.

### Far-Field Analysis

Not meaningful at N=100 — the far-field fitting range (r > 2*2*sigma = 8.0) exceeds
the usable box radius (L/2 - boundary ≈ 6.1). Higher resolution runs were not attempted
because no persistent resonance was found.

### Numerical Stability Note

The Test 1A bump with A=0.8 on N=200 grid was **numerically unstable** — energy diverged
to 10^8 within one diagnostic step. This is caused by extreme Skyrme forces (E4/E2 ~ 5)
on the steep Gaussian gradient. The A=0.3 bump (E4/E2 ~ 0.34) and all kick initializations
were stable. The instability is a CFL violation for the specific gradient structure, not
a physics result.

---

## Interpretation

### Why B=0 Oscillons Don't Form on S³ + Skyrme

The null result has a clean theoretical explanation:

1. **Derrick's theorem (static)**: For the sigma model + Skyrme Lagrangian, a static
   B=0 configuration has E = E_2 + E_4 with both terms positive. Under spatial rescaling
   x -> lambda*x: E_2 -> lambda*E_2, E_4 -> lambda^{-1}*E_4. The minimum is at
   E_2 = E_4 (finite size) — but ONLY if the configuration is topologically prevented
   from continuously deforming to vacuum. For B=0, the field can shrink E_2 and E_4
   to zero by smoothly unwinding to vacuum. No static bound state.

2. **Dynamic extension**: A time-dependent B=0 configuration has E_kin available as a
   third energy channel. Energy in the core can always escape as radiation because there
   is no topological barrier. The S³ constraint couples the three quaternion components,
   but this coupling is too weak (cubic in epsilon) to overcome the dispersive tendency
   of the wave equation. The three-wave resonance condition omega_1 + omega_2 = omega_3
   requires sustained amplitude to maintain the parametric coupling, but the amplitude
   itself disperses.

3. **Skyrme term is topology-specific**: The Skyrme quartic term prevents topological
   solitons (B != 0) from collapsing to zero size, but for B=0 it acts as an additional
   dispersive pressure — it ACCELERATES spreading, not prevents it.

### Implications for the v19 Program

- **Phase 2 is blocked**: Without a persistent resonance, there is no source to drive
  density depletion. The v6-style conserved density coupling requires a long-lived
  localized structure that continuously sources the density equation. A transient
  that disperses in ~5 time units cannot produce a steady 1/r depletion profile.

- **The "hydrocarbon chain" analogy fails**: In the NEXT_v19.md proposal, three-mode
  shared resonance was supposed to create higher energy density than dispersed vacuum.
  The numerical result shows the opposite — the S³ nonlinearity provides no binding
  mechanism for B=0 configurations. Energy concentration requires topology (B != 0).

- **Far-field 1/r from constraint**: The hypothesis that S³ incompressibility would
  force long-range field adjustments is not supported. The constraint is satisfied
  locally at each point; it does not create global pressure. The time-averaged DC
  rectification (Section 2.9 of the method) does exist but falls as r^{-2n} where
  n >= 2 (the mode amplitude falloff), giving r^{-4} or steeper — not 1/r.

---

## What Would Be Needed for a Positive Result

For a B=0 resonance to persist on S³, one would need:

1. **A confining mechanism beyond the Skyrme term**. Options:
   - Pion mass term (m_pi > 0): creates a mass gap, but also makes the far field Yukawa
     (exponential), killing any hope of 1/r.
   - Higher-derivative stabilizer (L_6 sextic): already explored in v2, provides more
     confinement but still topology-dependent.
   - External potential V(q): breaks the sigma-model structure.

2. **A conserved charge for B=0 configurations**. The Skyrmion is stable because B is
   a topological invariant. For B=0, there is no conserved quantity preventing dispersal.
   A Noether charge (like Q-ball charge) could work but requires a U(1) symmetry not
   present in the S³ sigma model.

3. **A fundamentally different field theory**. The S³ sigma model + Skyrme is optimized
   for topological solitons. A three-body resonance might require a field theory with
   intrinsic three-wave coupling (e.g., a non-abelian gauge theory with three-gluon
   vertex, or a conformal field theory with specific OPE structure).

---

## Files

| File | Description |
|------|-------------|
| `src/radial1d.c` | Phase 0: 1D radial Skyrme PDE solver |
| `src/triad3d.c` | Phase 1: 3D triad resonance simulator |
| `data/phase0_timeseries.tsv` | 1D time series: f(0,t), energies |
| `data/phase0_spectrum.tsv` | 1D frequency spectrum |
| `data/test{N}_timeseries.tsv` | 3D time series for test N |
| `data/test{N}_modes.tsv` | Mode projections and phase-lock diagnostic |
| `data/test{N}_tavg_profile.tsv` | Time-averaged radial profiles |
| `data/test{N}_spectrum.tsv` | DFT power spectrum of q_1(origin) |
| `data/test{N}_farfield.txt` | Far-field power-law fits |
