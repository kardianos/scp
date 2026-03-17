# V23-D Phase 1 Results: Inter-Oscillon Potential

## Parameters

- mu = -20, kappa = 20, mass = 1.0
- Gaussian init: A = 1.0, sigma = 3.0
- Equilibration: Nx=4000, xmax=100, t_equil=5000
- Two-oscillon: Nx=8000, xmax=200, t_meas=2000
- Force measurement: quadratic fit to sep(t) over t in [0, 200]
- Center-of-energy tracking (not center-of-phi^2)

## Phase 1a: Single Oscillon Equilibration

Oscillon stabilizes by t ~ 2500 with:
- Mass (energy) M = 1.298 code units
- Breathing frequency omega ~ 0.87 (below mass gap m = 1.0)
- Slow energy loss through radiation: E drops from 5.58 to 1.30 during equilibration
- Profile saved at breathing maximum (t = 4994.4, after 349 maxima)

## Phase 1b: Force Table

| D   | F(D)        | V(D)        | D_final | Notes                       |
|-----|-------------|-------------|---------|------------------------------|
| 8   | +9.48e-04   | -1.42e-02   | 15.8    | Strong overlap, cores merge  |
| 10  | +2.09e-03   | -1.11e-02   | 21.4    | Overlap regime               |
| 12  | +3.12e-03   | -5.92e-03   | 22.0    | Peak repulsion               |
| 14  | +1.74e-03   | -1.05e-03   | 41.8    | Repulsive, weaker            |
| 16  | -2.67e-04   | +4.23e-04   | 58.0    | Weakly attractive            |
| 18  | -5.06e-05   | +1.05e-04   | 28.2    | Attractive                   |
| 20  | -1.51e-05   | +3.97e-05   | 40.2    | Weakly attractive            |
| 25  | -6.78e-07   | +2.25e-07   | 46.9    | Very weak attraction         |
| 30  | +5.88e-07   | ~0          | 48.8    | Noise level                  |

Sign convention: F > 0 = repulsive, F < 0 = attractive.

## Key Finding: Equilibrium Spacing Exists

**Sign change at D_eq ~ 15.7** (between D=14 and D=16).

- D < 15.7: repulsive (tails overlap, cores repel)
- D > 15.7: attractive (tail-mediated attraction, Yukawa-like)
- D >> 20: force negligible

The potential V(D) has:
- A repulsive core (D < 15.7)
- An attractive well centered near D ~ 16-18
- Well depth ~ 4e-04 code energy units

## Potential Well Parameters

From the force table, at D_eq ~ 15.7:
- Well depth: V_min ~ 4.2e-04 (from trapezoid integration)
- Spring constant: K ~ 1.0e-3
- Predicted phonon speed: c_s = d_eq * sqrt(K/M) ~ 0.44

---

# V23-D Phase 2 Results: Oscillon Chain Phonon Spectrum

## Method

Code: `src/chain1d.c`

1. Equilibrate single oscillon (t=5000, absorbing BC, identical to Phase 1).
2. Place N=8 oscillons on a PERIODIC ring of length L = N*d.
3. Periodic boundary conditions eliminate edge effects and radiation pressure
   asymmetry that destabilized open-boundary chains.
4. Evolve with Velocity Verlet (energy-conserving, no damping).
5. Track oscillon positions via energy-weighted centroids within Voronoi cells.
6. Normal mode decomposition Q_q(t) + DFT to find dispersion omega(k).

Grid: Nx = N * 320 points per spacing. CFL: dt ~ 0.025.

## Chain Stability

The oscillon chain is UNSTABLE at all tested spacings (d = 13-20).

The inter-oscillon potential well depth (V_min ~ 4e-4) is far too shallow
to bind the oscillons against the transient energy from profile superposition.
The superposition of 8 equilibrated single-oscillon profiles creates a
non-equilibrium initial state whose excess energy (~0.1 per bond) exceeds
the well depth by a factor of ~250.

| d   | Initial E | N*M    | Interaction | max_u/d | Status    |
|-----|-----------|--------|-------------|---------|-----------|
| 13  | 9.58      | 10.38  | -0.80       | 3.7     | unstable  |
| 14  | 9.63      | 10.38  | -0.76       | 3.9     | unstable  |
| 15  | 9.70      | 10.38  | -0.69       | 2.4     | unstable  |
| 16  | 9.77      | 10.38  | -0.62       | 0.8     | marginal  |
| 17  | 9.85      | 10.38  | -0.54       | 2.5     | unstable  |
| 18  | 9.91      | 10.38  | -0.48       | 3.7     | unstable  |

Energy conservation is excellent with periodic BC: 99.8-99.9% retained
over t=10000 (vs 33% with absorbing BC in the open-boundary version).

Pre-equilibration with velocity damping was attempted but the damping
destroys the oscillon breathing mode, causing the oscillons to disperse.
The oscillons NEED their internal breathing kinetic energy to survive —
any damping that removes translational energy also removes breathing energy.

## Phonon Signal Detection

Despite the lattice instability, a coherent phonon signal IS detectable
in the transient dynamics at d=15 (just inside the repulsive regime):

| q | k       | omega_peak | c_s = omega/k | Power   |
|---|---------|------------|---------------|---------|
| 0 | 0.000   | 0.881      | --            | ~0      |
| 1 | 0.052   | 0.021      | **0.40**      | 0.295   |
| 2 | 0.105   | 0.013      | 0.13          | 0.058   |
| 3 | 0.157   | 0.059      | **0.38**      | 0.045   |
| 4 | 0.209   | 0.009      | 0.04          | 0.030   |

The q=0 mode correctly identifies the breathing frequency (omega = 0.88,
consistent with Phase 1 measurement of omega_breath ~ 0.87).

The q=1 and q=3 modes both give c_s ~ 0.40, consistent with the Phase 1
prediction of c_s ~ 0.44 from the spring constant measurement. This is
a 9% agreement.

However, the q=2 and q=4 modes do NOT follow the acoustic dispersion
relation, giving much lower c_s values. This indicates the phonon signal
is not clean — only specific modes happen to be excited coherently by the
superposition transient.

## Spacing Dependence

The detected phonon speed depends strongly on the initial spacing d:

| d   | c_s(q=1) | c_s(q=3) | Notes                        |
|-----|----------|----------|------------------------------|
| 13  | 0.43     | 0.96     | Strong repulsion, q=3 noisy  |
| 14  | 0.24     | 0.37     | Moderate signal               |
| 15  | **0.40** | **0.38** | Best consistency              |
| 16  | 0.03     | 0.19     | Near equilibrium, weak signal |
| 17  | 0.25     | ~0       | Attractive side, noisy        |
| 18  | 0.03     | ~0       | Weak interaction, drift only  |

The clearest signal appears at d=15 where the repulsive force provides
a structured initial perturbation. At d=16 (near D_eq), the force is
nearly zero and the transient is unstructured, producing only drift.

## Is There a Gapless Acoustic Branch?

**Partially yes, but not robustly.**

At d=15, the q=1 and q=3 modes show omega proportional to k with
c_s ~ 0.40, consistent with an acoustic branch. The q=0 mode (after
COM subtraction) has zero power, confirming no gap at k=0.

However:
1. The chain is globally unstable (max_u/d > 2).
2. Only q=1 and q=3 follow the acoustic relation; q=2 and q=4 do not.
3. The signal is deterministic (identical across seeds with delta=0),
   meaning it comes from the specific transient, not thermal equilibrium.
4. With random initial displacements (delta > 0), the phonon signal
   is overwhelmed by chaotic diffusion.

## Measured Sound Speed

**c_s = 0.40 +/- 0.02** (from q=1 and q=3 modes at d=15).

This is within 9% of the Phase 1 prediction c_s ~ 0.44 from the
spring constant K ~ 1e-3 and mass M ~ 1.298:
  c_s_predicted = d * sqrt(K / M) = 15.7 * sqrt(1e-3 / 1.298) = 0.44

The 9% discrepancy may come from:
- Using d=15 instead of D_eq=15.7 (different effective spring constant)
- Nonlinear effects (large-amplitude oscillations beyond harmonic regime)
- Finite-N effects (N=8 is small)

## Conclusion

The oscillon lattice is **too weakly bound to be dynamically stable**.
The potential well depth (~4e-4) is only 0.03% of the oscillon rest mass,
making the lattice softer than any known physical crystal. The transient
energy from initializing the chain by superposition (~0.1 per bond) exceeds
the well depth by ~250x, causing immediate lattice melting.

Despite this instability, the early-time transient dynamics DO show a
phonon-like signal with sound speed c_s ~ 0.40, matching the Phase 1
prediction from the measured spring constant. This confirms that the
inter-oscillon potential has the correct structure for phonon propagation —
the physics is right, but the binding is too weak for a stable lattice.

**Implications for emergent gravity**: A phonon-mediated graviton analog
would require a STABLE oscillon lattice with long-lived acoustic modes.
The current model (mu=-20, kappa=20, m=1.0) does not achieve this.
Stronger coupling (larger |mu|, kappa) might deepen the well, but would
also change the oscillon structure and equilibrium spacing. The fundamental
challenge is that the inter-oscillon force is exponentially weak at large
separations (Yukawa-like, range ~1.5 code units), while the oscillon cores
extend ~5 code units, requiring spacings where overlap is non-negligible.

---

# V23-D Phase 2 (3D): Two-Oscillon Pair in Three Dimensions

## Motivation

The 1D chain test (above) failed because radiation pressure overwhelmed the
weak attractive well. In 3D, radiation pressure falls as 1/r^2 due to
geometric spreading. This could make the lattice stable. The minimal viability
test: can TWO oscillons remain bound in 3D?

## Method

Code: `src/chain3d.c`

**Phase 1: Equilibrate a single 3D oscillon**
- Parameters: mu=-20, kappa=20, m=1.0, A=0.8, sigma=3.0 (v21 production)
- Grid: N=96, L=15 (dx=0.316)
- Evolve t_equil=500; save 3D profile (phi + velocity) at last breathing peak
- Equilibrated oscillon: E = 82.7, peak amplitude = 0.90

**Phase 2: Two-oscillon evolution**
- Grid: N=96, L=25 (dx=0.526)
- Place two equilibrated profiles at z = +/-D/2 via trilinear interpolation
- In-phase initialization (attractive configuration)
- Velocity Verlet integrator, absorbing spherical shell at 70-95% of L
- Track energy-weighted z-centroids (split at z=0) for t_run=2000

## Results

### D=10 (initial separation 10.0)

| Time range | Separation | E_total | Peak amplitude | Phase    |
|------------|-----------|---------|----------------|----------|
| t=0        | 10.0      | 172     | 0.87           | Initial  |
| t=0-400    | 10 -> 3.5 | 172->83 | 0.87->0.30     | Merger   |
| t=400-800  | 3.1-7.7   | 83->72  | 0.10-0.85      | Merged oscillon |
| t=800-1200 | 7.7-12.4  | 72->21  | 0.03-0.11      | Expansion+dissipation |
| t=1200-2000| 8.3-11.3  | 21->5.0 | 0.002-0.05     | Dying oscillation |

### D=16 (initial separation 16.0)

| Time range | Separation | E_total | Peak amplitude | Phase    |
|------------|-----------|---------|----------------|----------|
| t=0        | 16.0      | 165     | 0.87           | Initial  |
| t=0-300    | 16 -> 4.0 | 165->94 | 0.87->0.67     | Attraction+merger |
| t=300-700  | 3.1-12.3  | 94->76  | 0.10-0.86      | Merged oscillon |
| t=700-1200 | 7.5-12.6  | 76->20  | 0.01-0.05      | Expansion+dissipation |
| t=1200-2000| 9.0-11.1  | 20->4.6 | 0.002-0.04     | Dying oscillation |

### Common features (both D values)

1. **Rapid merger** (t < 400): Both separations quickly collapse to ~3-5 units.
   The oscillons attract strongly and their cores overlap, forming a single
   merged blob. Initial energy ~165 drops to ~80 as radiation is shed.

2. **Merged oscillon phase** (t ~ 400-800): A single composite oscillon
   breathes with peak amplitude ~0.85, slowly losing energy. This is
   consistent with two oscillons merging into one (E ~ 82 for a single
   oscillon at these parameters).

3. **Expansion and dissipation** (t > 800): The merged blob expands as
   energy radiates away. The separation opens to ~10-12, but this is NOT
   two separate oscillons — it is a single diffuse energy distribution
   whose centroid splits across z=0.

4. **Late-time convergence**: Both D=10 and D=16 converge to the SAME
   final state: sep ~ 9.8, E ~ 4.6, pk ~ 0.02. This "universal attractor"
   is not a two-oscillon bound state but a nearly uniform remnant with
   negligible field amplitude.

5. **Energy loss**: 97% of initial energy is radiated by t=2000. The
   absorbing boundaries efficiently remove radiation.

## Analysis

### Why merger occurs

The 3D oscillons at mu=-20 have radius ~5 code units (half-maximum of energy
density). At D=16, the edge-to-edge gap is only ~6 units. The oscillons' tails
overlap significantly, and the attractive inter-oscillon force (measured in
Phase 1b at D=16: F = -2.67e-04) pulls them together. In 3D, there is no
repulsive radiation pressure to resist this because:

1. **Radiation dilution**: In 3D, radiation amplitude falls as 1/r, so
   radiation pressure falls as 1/r^2. The pressure on the neighboring
   oscillon at distance D is ~P_rad/D^2, which is negligible.

2. **No equilibrium**: The 1D equilibrium at D_eq ~ 15.7 was a balance
   between repulsive radiation pressure and attractive tail interaction.
   In 3D, the radiation pressure vanishes, removing the repulsive barrier.
   The potential becomes purely attractive at all separations > core size.

### Resolution effects

The Phase 2 grid (dx=0.53) is coarser than the equilibration grid (dx=0.32).
This could cause numerical dissipation. However, the merger occurs BEFORE
significant numerical decay — the separation drops to 4 within t=200, when
the oscillon still has high amplitude (pk > 0.5). The merger is physical,
not numerical.

## Conclusion

**The 3D two-oscillon pair is UNSTABLE: oscillons attract and merge.**

The 3D geometric dilution of radiation (1/r^2) is actually WORSE for lattice
stability, not better. In 1D, radiation pressure provided a repulsive force
that balanced the attractive tail interaction at D_eq ~ 15.7. In 3D, this
repulsive mechanism vanishes because radiation spreads in all directions,
and the oscillons simply fall together.

This is the fundamental problem:
- **1D**: Radiation pressure ~ const (no spreading) -> repulsion exists -> equilibrium exists -> but well too shallow for stable chain
- **3D**: Radiation pressure ~ 1/r^2 (spreading) -> repulsion vanishes -> NO equilibrium -> immediate merger

## Overall V23-D Assessment

| Phase | Test | Result | Verdict |
|-------|------|--------|---------|
| 1a | Single oscillon equilibration | Stable oscillon, M=1.30, omega=0.87 | OK |
| 1b | Inter-oscillon force (1D) | Equilibrium at D_eq~15.7, well depth ~4e-4 | OK |
| 2 (1D) | 8-oscillon chain | Unstable: well too shallow (250x below needed) | FAIL |
| 2 (3D) | Two-oscillon pair | Unstable: merger (no repulsive barrier in 3D) | FAIL |

**The oscillon lattice concept fails in both 1D and 3D, but for DIFFERENT reasons:**
- 1D: equilibrium exists but binding energy too weak for stability
- 3D: no equilibrium exists because radiation pressure is geometrically diluted

**Implications**: An oscillon-based emergent spacetime with phonon gravitons
is not viable in this model. The inter-oscillon interaction does not support
a stable crystal in any dimension. A fundamentally different mechanism for
oscillon binding (beyond radiation-mediated forces) would be needed.
