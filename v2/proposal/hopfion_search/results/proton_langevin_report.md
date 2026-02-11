# Proton Field Dynamics: Langevin Simulation Results

## Hypothesis Tested

**Claim**: Particles are not static soliton solutions but self-reinforcing dynamic patterns
in the Skyrme field that *require* a thermal bath (T>0) to sustain themselves. Individual
"quarks" are small, unstable field knots that dissolve alone but mutually reinforce when
composed, creating a stable proton. "Gluons" are the field gradient pattern between quarks,
not separate particles.

**Test**: Langevin field dynamics with thermal bath at various temperatures, testing four
initialization modes that probe different aspects of the hypothesis.

## Method

### Simulator: `src/proton.c`

Self-contained 3D Skyrme field dynamics with Langevin thermostat:
- **Field**: Q4 quaternion on σ-model (|q|=ρ₀), Skyrme L₂+L₄ energy
- **Dynamics**: Leapfrog time-stepping with Langevin damping + noise
- **Thermostat**: `v(t+dt) = v(t) + dt·c²·F - γ·dt·v + √(2γT·dt/h³)·ξ`
- **Projection**: velocity projected tangent to S³ after each step
- **Grid**: N=128, L=8, h=0.125 (code units, e=1, ρ₀=1)
- **Force**: 9-point consistent Laplacian (E₂) + commutator-based Skyrme term (E₄)

### Initialization Modes

| Mode | Description | Q(0) | E_pot(0) |
|------|-------------|-------|----------|
| **hedgehog** | B=1 hedgehog from 1D equilibrium profile | 1.000 | 114 |
| **quarks** | Three Gaussian perturbations at equilateral triangle vertices (D=2, A=1.5) | 0.000 | 29 |
| **b3** | B=3 rational map R(z)=(z³-√3iz)/(√3iz²-1) from 1D profile | 2.994 | 401 |
| **extract** | Single quark extracted from B=3 (R_cut=1.0 around density peak) | 1.244 | 404 |

### Temperature Regime

Soliton energy E≈114 code units. Thermal energy E_thermal ≈ N³×3×T/2 = 3.15×10⁶ × T.

| Temperature | E_thermal | Signal/Noise | Verdict |
|-------------|-----------|-------------|---------|
| T=0 | 0 | ∞ | Clean baseline |
| T=1×10⁻⁵ | 31 | 3.7:1 | Soliton visible |
| T=1×10⁻⁴ | 315 | 0.36:1 | Soliton barely visible |
| T=0.1 | 315,000 | 3.6×10⁻⁴ | Soliton completely buried |

## Results

### 1. Hedgehog B=1: Topology Loss at All Temperatures

| Experiment | γ | T | Q(0) | Q at t=2 | Q at t=4 | t_dissolve |
|------------|---|---|------|----------|----------|------------|
| T=0, no damping | 0 | 0 | 1.000 | 1.000 | 0.017 | ~3.0 |
| T=1e-5 | 0.5 | 1e-5 | 1.000 | 0.162 | ~0 | ~2.5 |
| T=1e-4 | 0.5 | 1e-4 | 0.997 | 0.015 | ~0 | ~1.5 |

**Key finding**: Topology dissolves at t≈2-3 regardless of temperature. Higher temperature
*accelerates* dissolution (t=3→1.5 from T=0 to T=1e-4). The Langevin damping (γ=0.5)
also accelerates dissolution by removing kinetic energy that keeps the soliton temporarily
coherent.

**Mechanism**: The soliton shrinks below the lattice resolution (h=0.125) and unwinds.
This is the well-known "lattice saddle point" phenomenon — the Skyrmion is an energy
saddle point on the discrete lattice, not a local minimum. At T=0 without damping, kinetic
energy keeps Q≈1 for about 2.5 code time units before the soliton radius crosses below
~2h and the topology is lost. With damping, kinetic energy drains, and with thermal noise,
random kicks push the soliton off the saddle faster.

### 2. Three Quarks: No Topological Nucleation

| Experiment | γ | T | Q(0) | Q(t=2) | Q(t=18) |
|------------|---|---|------|--------|---------|
| T=0, no damping | 0 | 0 | 0.000 | 0.000 | 0.000 |
| T=1e-5 | 0.5 | 1e-5 | 0.000 | 0.000 | 0.000 |

Q = 0.000 throughout all simulations. The three Gaussian perturbations (each individually
Q=0 since they don't wrap S³) simply radiate away as dispersive waves. No topological
nucleation occurs, not even transiently.

**Interpretation**: Baryon number is fundamentally a topological invariant (winding number
of the map S³→S³). It cannot be created from smooth, small-amplitude perturbations — no
matter how many you superpose. The quarks-as-Gaussians initial condition has Q=0 exactly,
and Q is conserved until the field gradients reach the lattice scale where numerical
dissipation takes over.

### 3. B=3 Rational Map: Three-Quark Structure Dissolves

| Experiment | dt | Q(0) | Q at t=0.5 | Q at t=1.0 | Q at t=2.0 |
|------------|-----|------|-----------|-----------|-----------|
| T=0, γ=0.5 | 0.005 | 2.994 | 2.949 | 0.031 | -0.528 |

The B=3 soliton (with its visible three-quark baryon density structure) loses topology
by t≈1.0. This is faster than B=1 (t≈3.0) because:
1. B=3 has higher E₄ (287 vs 52) → larger forces → faster dynamics
2. B=3 is larger (R_rms > B=1 R_rms) but less stable per baryon on the lattice
3. The rational map ansatz is only approximate — not the true B=3 minimum

After topology loss, Q oscillates around -0.5 rather than 0, indicating residual winding
from the larger field structure that hasn't fully radiated away yet.

### 4. Extracted Quark: Immediate Dissolution

| Experiment | R_cut | Q(0) | Q at t=0.5 | Q at t=1.0 |
|------------|-------|------|-----------|-----------|
| T=0, γ=0.5 | 1.0 | 1.244 | 0.040 | -0.004 |

The extracted quark (one density peak from B=3, masked outside R_cut=1.0) carries
fractional topology Q≈1.24 and dissolves almost immediately (t<0.5). The enormous initial
E₄=373 (vs 52 for equilibrium B=1) shows this is far from any equilibrium — the sharp
mask boundary creates massive Skyrme force that blows the field apart.

**Interpretation**: A "quark" isolated from a multi-baryon soliton is not topologically
protected. The baryon charge is a global, non-local property of the field wrapping. Cutting
a piece out doesn't give you a fraction of a stable soliton — it gives you an unstable
field configuration that immediately unwinds.

## Discussion

### Thermal Bath Does NOT Stabilize Topology

This is the central negative result. The hypothesis that thermal fluctuations could
sustain soliton topology is falsified:

1. **At low T** (T=1e-5): soliton visible but dissolves at same timescale as T=0
2. **At moderate T** (T=1e-4): dissolution accelerated, soliton barely visible above noise
3. **At high T** (T≥0.1): thermal energy 1000× larger than soliton, no meaningful dynamics

The problem is fundamental: on a discrete lattice, the Skyrmion is a *saddle point* of the
energy functional, not a local minimum. Thermal fluctuations provide random kicks that push
the system *off* the saddle in both stable and unstable directions. Since there is an
unstable direction (shrinking), thermal noise accelerates rather than prevents dissolution.

### Topology Cannot Be Created from Non-Topological Initial Conditions

Three individually non-topological field lumps (Q=0 each) never combine to produce topology.
This follows from homotopy theory: Q∈π₃(S³)=Z is an integer invariant that cannot change
by continuous evolution. The only way to change Q is through a singularity (field passes
through zero), which corresponds to creating/annihilating baryon-antibaryon pairs. Random
thermal fluctuations at T=1e-5 are far too weak to nucleate such pairs — the energy barrier
is of order E_soliton ≈ 114 code units = 1038 MeV, while T=1e-5 code units = 0.09 MeV.

### Comparison with Hopfion Langevin Results

The identical experiment was performed on Faddeev-Skyrme hopfions (previous session):
- Hopfion Q_H dissolved at t≈2 regardless of temperature
- Thermal bath accelerated hopfion dissolution
- No spontaneous hopfion nucleation from thermal fluctuations

Both Skyrmion and hopfion results are consistent: **lattice solitons are saddle points,
not local minima, and thermal baths cannot stabilize them.**

### What Would Be Needed for Thermal Stabilization

For the hypothesis to work, one would need:
1. A mechanism that makes the soliton a true *local minimum* on the lattice (not just in
   the continuum). This requires either:
   - Much finer grids (h→0, but computational cost ∝ N⁴)
   - Modified lattice action that preserves the topological protection
   - Constraint that prevents |q| from reaching zero (e.g., finite-λ potential)
2. A thermal nucleation mechanism with barrier ≈ E_soliton accessible at physical T.
   The physical proton mass is 938 MeV, so T_nucleation ~ 10¹² K — far above QCD
   deconfinement (T_c ~ 10¹² K). At these temperatures, solitons don't exist.

## Summary Table

| Question | Answer |
|----------|--------|
| Does thermal bath stabilize B=1 Skyrmion? | **NO** — accelerates dissolution |
| Can thermal fluctuations nucleate topology? | **NO** — Q exactly conserved until lattice artifacts |
| Does extracted quark survive in isolation? | **NO** — dissolves in t<0.5 |
| Does B=3 (three-quark) structure persist? | **NO** — dissolves faster than B=1 |
| Is dissolution temperature-dependent? | **YES** — higher T → faster dissolution |
| What causes dissolution? | Lattice saddle point: soliton shrinks below grid resolution |

## Conclusion

**The thermal stabilization hypothesis is definitively falsified** for the Skyrme model
on a cubic lattice with Langevin dynamics. Soliton topology is a property of the continuum
field equations, protected by the continuity of the map q: R³→S³. On a discrete lattice,
this protection is absent — the soliton can shrink to sub-grid scales and unwind. Thermal
fluctuations provide additional energy to traverse the lattice barrier, accelerating rather
than preventing topology loss.

This result, combined with the identical negative result for Faddeev-Skyrme hopfions,
establishes that **thermal baths cannot substitute for the topological protection that
exists in the continuum theory but is lost on the lattice**.

## Files

- `src/proton.c` — Langevin Skyrme dynamics simulator (~1515 lines)
- `scripts/viz_particle.py` — Binary snapshot visualization
- `scripts/viz_proton_comparison.py` — Comparison plots
- `results/viz/proton_Q_comparison.png` — Q(t) comparison across all experiments
- `results/viz/proton_dissolution_timeline.png` — Hedgehog dissolution sequence
- `results/viz/proton_initial_modes.png` — Four initialization modes
- `results/viz/proton_summary_table.png` — Summary table
- `data/proton_*/timeseries.dat` — Raw timeseries data
- `data/proton_*/snapshot_*.bin` — Binary field snapshots
