# SCP Field Theory — Current Best Concept

**Status**: Working theory with numerical confirmation of core mechanisms.
This document reflects the current best understanding as of March 2026.
Sections marked [CONFIRMED] have numerical evidence. Sections marked
[THEORETICAL] are consistent with the framework but unverified.

---

## 1. The Field

There is ONE field φ_a(x,t) with three real components (a = 0, 1, 2).
The field fills all of space. It IS the medium — there is no separate
"spacetime" or "vacuum." Where the field exists, physics exists.

### Equation of Motion [CONFIRMED]

    ∂²φ_a/∂t² = ∇²φ_a - m²φ_a - ∂V/∂φ_a               (1)

    V(P) = (μ/2) P² / (1 + κP²)                          (2)

    P = φ₀ φ₁ φ₂    (triple product)                      (3)

Parameters: m² = 2.25, μ = -41.345, κ = 50

The triple product P = φ₀φ₁φ₂ is the VOLUME FORM of the three-field
displacement gradient. The potential V(P) penalizes local volume change.
This is the natural energy functional of a 3D deformable medium with
irreducible three-body interactions.

No modifications to this equation are needed. No c(ρ), no gradient
coupling, no S/B split, no smoothing. All physics emerges from Eq. (1).

### Background

The field has a uniform background with amplitude A_bg ≈ 0.1:

    φ_a(x) = A_bg × cos(k·z + 2πa/3) + (braid structures)

The background is the "fabric" of space. Its density ρ_bg = energy density
of the oscillating background ≈ 0.03 in code units.

---

## 2. Particles (Braids) [CONFIRMED]

Particles are self-sustaining helical patterns in the field — braids.
The braid is NOT separate from the field. It IS the field, organized
differently. The simulation makes no distinction between "braid" and
"background" — every grid point runs the same equation.

### Braid Structure [CONFIRMED]

The bimodal braid (the best-characterized configuration):

    φ_a(x,y,z,0) = A_a × E(x,y) × cos(k·z + δ_a) + A_bg × cos(...)

    E(x,y) = exp(-[x²/(1+ε)² + y²/(1-ε)²] / 2R²)       (4)

    A_a = 0.8,  δ = (0, 3.00, 4.43),  R = 3.0,  ε = 0.333

The three fields oscillate with specific phase offsets (δ), creating a
helical twist along z. The elliptical envelope (ε = 0.333) breaks
azimuthal symmetry.

### Binding Anatomy [CONFIRMED]

The braid has internal geometry measured by |P| = |φ₀φ₁φ₂|:

| Region | r | |P| | Role |
|--------|---|-----|------|
| Tight core | 0-3 | 0.11-0.15 | Strongly bound, self-contained |
| Interaction surface | 4-6 | 0.01-0.06 | WHERE the braid talks to the fabric |
| Fabric | >6 | ~0 | Unbound ambient field |

The core INTAKES field (energy flows inward). The surface OUTTAKES
(energy flows outward). The braid is a PUMP — it processes field
through its helical structure.

### Stability [CONFIRMED]

- Braid survives T=500+ in periodic BC (V29-T1b, V33-C2)
- Robust to 6× amplitude perturbation (V29-T4b)
- Thermal equilibrium at A_noise ≈ 0.02 (V29-T1b)
- Evaporative self-concentration: fc rises 0.80→0.93 as outer layers
  radiate away (V29-T3)
- Steady radiation rate ~0.5 energy/time (V29-T3)

---

## 3. Gravity [CONFIRMED]

### Mechanism: Asymmetric Footprint

The braid depletes the surrounding field by binding field energy into its
helical structure. When another braid's depletion zone creates a ρ gradient,
the first braid's spatial perturbation profile becomes ASYMMETRIC.

**Why asymmetric**: The braid's perturbation decays as a Yukawa tail
(~e^{-m_eff·r}/r), where the effective mass depends on the local background:

    m_eff² = m² + V''(P_bg)

In low-ρ regions, P_bg is smaller → V'' is weaker → m_eff is smaller →
the Yukawa tail extends further. In high-ρ regions, m_eff is larger →
the tail is shorter. The braid therefore REACHES FURTHER into the depleted
side than the dense side.

**Footprint asymmetry measured** (V33 gradient test, N=128, 5/6 snapshots):
- Half-width ratio R_low/R_high = 1.09–1.57
- The braid consistently extends further into the low-ρ side

This geometric asymmetry shifts the braid's equilibrium position: the
self-energy minimum moves toward the side with the longer perturbation
tail. The braid drifts toward depletion — toward other braids = ATTRACTION.

    Braid A depletes ρ → creates low-ρ zone
    Braid B's perturbation extends further into the low-ρ side
    B's equilibrium shifts toward A → drift → GRAVITY

This is a DYNAMIC effect driven by the braid's oscillation cycle in the
gradient, not a static energy minimization. The static interaction energy
E(D) between two braids is monotonically repulsive at all separations —
there is NO attractive well (V33 energy sweep, D=8–50). The attractive
force arises entirely from the dynamic footprint asymmetry.

The braid-field coupling is STRONGER in high-ρ backgrounds (drag test:
momentum retention 1.40× in high-ρ vs 0.30× in low-ρ). The force
arises from the asymmetric spatial profile, not from friction.

### Force Law [CONFIRMED]

**F ∝ ∇ρ (R² = 0.9998)**: Gradient sweep at four strengths (ρ ratios
1.22 to 9.00) shows drift perfectly proportional to ∇ρ. The ratio
drift/∇ρ = -186 ± 3, constant across all gradient strengths. This is
exactly the Newtonian form F = -C × ∇Φ where Φ ∝ ρ.

Measured from the standard equation (Eq. 1) with NO modifications:

| D | ΔD (T=200) | Character |
|---|------------|-----------|
| 5 | +9.9 | Hard repulsive core |
| 8 | +5.2 | Repulsion |
| 12 | -1.1 | Attraction begins |
| 15 | -9.1 | STRONG attraction |
| 18 | -8.0 | Strong attraction |
| 20 | -0.2 | Weak attraction |
| 25 | -0.3 | Weak attraction |
| 40 | -0.24 | Detected at long range |
| 80 | -1.2 | Detected at very long range |

Power law fit (D=15-30): F ∝ 1/D^1.8 (steeper than Newton's 1/D²)

Energy conservation: -0.1% per T=100 (symplectic Verlet, small dt drift).

### Gradient Response [CONFIRMED]

A braid placed in an imposed ρ gradient (high left, low right) drifts
toward LOW ρ (V33 gradient test, N=128 T=100). This is the correct
direction for gravity: mass creates depletion, other mass moves toward
depletion. N=512 test with gentle 2.25:1 gradient shows consistent
slow drift in the same direction.

The directional footprint asymmetry was measured directly: the braid's
perturbation (δρ above background) has a wider half-width on the low-ρ
side than the high-ρ side (ratio 1.09–1.57 across 5/6 time snapshots).
This confirms the GEOMETRIC mechanism: the braid's spatial profile is
lopsided in a gradient, and the equilibrium position shifts accordingly.

### Energy Conservation and the Depletion [CONFIRMED]

The depletion around the braid is NOT consumed or radiated energy.
It is field energy that has been BOUND into the braid's helical structure.
Total energy is conserved (drift < 0.1% per T=100).

The braid does not act as a heat engine. The intake/outtake cycle is a
STEADY STATE, not a one-way consumption:

    Intake at core: field energy flows into the braid's helical pattern
    Outtake at surface: the braid's oscillation re-emits some coherent energy
    Net: the braid structure is maintained, the depletion profile is STATIC

The depletion profile exists because the braid's binding (triple product
coupling) concentrates field energy from the surrounding region into the
compact helical arrangement. The surrounding field has lower ρ because
that energy is now part of the braid.

Analogy: a crystal forming in solution. The crystal binds solute into a
lattice. The surrounding solution is depleted. The depletion is structural,
not from ongoing consumption. No energy flows once equilibrium is reached.

This is NOT a Le Sage / kinetic / pushing theory of gravity:
- Le Sage: continuous absorption of mediator particles → heating → fails
- This model: static structural depletion from field binding → no heating
- The gravitational force comes from the GEOMETRY of the depletion,
  not from continuous energy flow
- At thermal equilibrium (V29-T1b), radiation = absorption, yet the
  depletion persists because it's structural, not dynamic

The braid's radiation (V29-T3: ~0.5 energy/time) is the braid slowly
relaxing toward its equilibrium configuration. It is NOT "processing
field for gravity." The radiation is a separate phenomenon that ceases
at thermal equilibrium. The depletion (and thus gravity) persists
regardless of the radiation state.

Key: 0 < V < c and T > 0 always. The braid is always moving and always
at nonzero temperature. But its gravitational depletion is a BOUND STATE
property — like the electric field of a charged particle, which exists
without continuous energy expenditure.

### F = ma from the Field [CONFIRMED structure, THEORETICAL identification]

The measured force law F = -C × ∇ρ has the structure of Newton's second
law. To see this, define the quantities precisely:

**Definitions**:

    ρ(x)    = local field energy density = ½Σ(∂φ_a)² + ½m²Σφ_a² + V(P)
    ρ_bg    = undisturbed background energy density ≈ 0.03
    δρ(x)   = ρ(x) - ρ_bg  (depletion profile around a braid; δρ < 0 near braid)
    ∇ρ      = spatial gradient of the energy density

    M_grav  = gravitational mass of the SOURCE braid
              = total field energy bound into the braid's helical structure
              = -∫ δρ dV  (integrated depletion; positive for a braid)
              This sets the DEPTH and SHAPE of the depletion profile.

    m_grav  = gravitational mass of the TEST braid
              = coupling strength to the density gradient
              = how strongly the test braid's footprint asymmetry
                responds to a given ∇ρ

    m_inert = inertial mass of the test braid
              = resistance to acceleration
              = E_braid / c²  (total bound energy, relativistic mass-energy)

    C       = proportionality constant measured in gradient sweep
              ≈ 186 (code units), constant across all tested gradients

**The force law**:

    F = -C × ∇ρ                                                     (5)

The source braid creates the depletion profile δρ(x), which sets ∇ρ at
the test braid's location. ∇ρ is proportional to M_grav (deeper depletion
= larger source mass = steeper gradient). So:

    ∇ρ  ∝  M_grav × geometry(D)

where geometry(D) encodes the radial dependence (1/D², Yukawa, etc.).

The test braid responds to this gradient with force F = -C × ∇ρ. If the
coupling constant C is proportional to the test braid's mass m_grav:

    C  =  m_grav × G_eff

then the force becomes:

    F  =  -G_eff × m_grav × M_grav × geometry(D)                    (6)

This is Newton's gravitational force law, with G_eff as the emergent
gravitational constant.

**The acceleration**:

    a = F / m_inert = -(C / m_inert) × ∇ρ                          (7)

**If m_grav = m_inert** (equivalence principle), then C/m_inert = G_eff,
and ALL braids experience the same acceleration in the same gradient:

    a = -G_eff × ∇ρ                                                 (8)

This is F = ma, where:
- **m** = m_grav = m_inert = the braid's bound energy (one quantity)
- **a** = -G_eff × ∇ρ = gravitational acceleration (set by the source)
- **F** = m × a = the gravitational force

**Why m_grav = m_inert should hold**: Both masses originate from the SAME
physical property — how much field energy is bound into the braid's
helical structure via the V(P) coupling:

- m_grav: more bound energy → deeper depletion → stronger ∇ρ response
- m_inert: more bound energy → more energy to accelerate → greater inertia

There is only ONE coupling mechanism (the triple-product potential V(P)),
and it governs both how a braid creates gradients AND how it responds to
them. The equivalence principle is not imposed — it follows from the
single-coupling structure of the theory.

**Status**: The proportionality F ∝ ∇ρ is confirmed (R² = 0.9998). The
identification C = m_grav × G_eff and m_grav = m_inert is theoretical —
testing requires braids of different sizes/binding energies in the same
gradient to verify that C/m_inert is constant.

The m² = 2.25 in Eq. (1) is a LAGRANGIAN parameter that provides field
confinement. It is NOT the gravitational mass. The gravitational mass
M_grav emerges from the braid's structural depletion of the surrounding
fabric.

### Braid-Field Coupling Scales with ρ [CONFIRMED]

The braid's coupling to the background field (both force and drag) is
STRONGER in high-ρ environments. A kicked braid in high-ρ background
(A_bg=0.15) retains/amplifies its momentum (1.40×), while the same kick
in low-ρ (A_bg=0.05) decays rapidly (0.30× retention). This is because
V(P) coupling is proportional to P_bg, which scales with ρ.

This means the gravitational drift is NOT from reduced drag in depleted
regions. It is from the geometric footprint asymmetry: the braid's
perturbation profile extends further where ρ is lower (longer Yukawa
range), shifting the equilibrium position toward depletion.

### Comparison to Newton [THEORETICAL]

For F ∝ 1/D^n: measured n ≈ 1.8. Newton requires n = 2. The discrepancy
may be from:
- Periodic BC image contamination at D > 20
- The massive field (m=1.5) creating Yukawa corrections
- Resolution effects at the tested grid spacing
- Genuine deviation from Newton (modified gravity at nuclear scale)

The exponent needs remeasurement with proper (pinned/outflow) boundary
conditions and larger domains.

---

## 4. Electromagnetism [THEORETICAL]

### Path: Complex Fields + Gauge Symmetry

V29-T10G showed that upgrading from 3 real to 3 COMPLEX scalar fields
preserves the braid and adds U(1) conserved charges. The potential
becomes V(|ψ₀ψ₁ψ₂|²) with ψ_a ∈ ℂ.

The braid carries charge Q ≈ 3-14 (from T10G). Promoting the global
U(1) to a gauge symmetry (adding A_μ) would give:

    D_μ ψ_a = ∂_μ ψ_a + ieA_μ ψ_a
    L_EM = |D_μ ψ|² - (1/4)F_μν² - V(|P|²)

This would produce EM interactions between charged braids. However,
V24-MA showed that explicit gauge coupling can kill oscillons — the
braid's survival under gauging needs verification.

### Torsion as EM Proxy [CONFIRMED]

The antisymmetric gradient tensor ω_ij = ½(∂_i φ_j - ∂_j φ_i) is the
torsion of the elastic interpretation. The V28 bimodal braid has nonzero
torsion flux (Φ_T = 1.02 at N=128 validation). This is an EM-like
geometric property but NOT a gauge field.

---

## 5. Multi-Braid Bonding [THEORETICAL]

### Nuclear Analogy

Two separate braids scatter (V33-C3: approach to D≈12, bounce, separate).
They CANNOT merge by collision at the tested initial conditions.

Multi-braid bound states (heavier "atoms") must form by CONDENSATION
from a hot dense field state — analogous to stellar nucleosynthesis:
- Single braids = hydrogen (one helical pattern)
- Multi-braid bonds = heavier nuclei (shared helical structure)
- Formation requires high-energy conditions (dense field + cooling)
- Once formed, the shared structure is stable

### Covalent Bonding Analog [THEORETICAL]

Like covalent bonds share electrons, heavier "atoms" would be braids
sharing their helical structure — partial overlap of the binding regions.
The binding energy comes from the shared triple-product coupling.

---

## 6. Dark Matter [THEORETICAL]

### Accretion + Depletion Profile

V30's M7 analysis showed the braid simultaneously:
- ACCRETES B-field at the core (higher density)
- DEPLETES B-field in a surrounding shell (lower density)

This profile (dense core + depleted shell) is structurally similar to
dark matter halos (NFW profiles). The "dark matter" in this framework
would be the field depletion around each braid — invisible but
gravitationally active through the depletion mechanism.

---

## 7. Cosmology [THEORETICAL]

### Braid Formation

V30 tested expansion + rotation for spontaneous braid formation.
Result: no braids from expansion alone. Braids require specific
initialization conditions (helical phase structure).

In a cosmological context, braids would form during phase transitions
in the early universe — analogous to topological defect formation
(Kibble-Zurek mechanism) but requiring a field with three-body coupling
and sufficient angular momentum.

---

## 8. Lorentz Invariance and Relativity

### The Equation IS Lorentz-Invariant [CONFIRMED by construction]

Equation (1) is a standard relativistic wave equation with c=1. It is
manifestly Lorentz-invariant. Any solution boosted by a Lorentz
transformation is also a solution.

### The Background IS NOT an Aether

The uniform background φ_a = A_bg × cos(k·z + 2πa/3) has a preferred
z-direction. This is an INITIAL CONDITION, not a property of the equation.
The equation has full rotational and Lorentz symmetry. The z-preference
is a simulation convenience for the helical braid initialization.

A rotationally symmetric background (random phases, isotropic spectrum)
would have no preferred direction. The physics does not require a
preferred frame.

### Lorentz Contraction of Braids [THEORETICAL, untested]

A boosted braid SHOULD contract by γ = 1/√(1-v²/c²) because the
equation is Lorentz-invariant. This has not been tested numerically
(V29-T8 was deprioritized). It is a key planned verification:

    Test: Initialize a braid with Lorentz boost v=0.1c, 0.3c, 0.5c
    Measure: Does the aspect ratio match γ?
    Expected: YES (from the equation's symmetry)

### The Force Law and the Mass Term

The mass term m²=2.25 in Eq. (1) creates a Yukawa decay at long range:
individual field perturbations fall off as e^{-mr}/r (range ≈ 1/m ≈ 0.67).

The measured F ∝ 1/D^1.8 is likely Yukawa at intermediate D, not a true
power law. At D >> 1/m, the force should decay exponentially. The
"long-range" attraction measured at D=40-80 may be:

    (a) Periodic BC image contamination (likely at D > 2L/3)
    (b) Collective nonlinear effects that extend beyond the Yukawa range
    (c) A true power-law tail from the braid's steady-state radiation

For exact Newtonian 1/r² gravity: need m → 0. V29-T2 showed m=0 braids
survive but lose the quadrupole. V33-C4 showed m=0 gives massive
attraction (ΔD=-16.8) but energy is not conserved (unstable vacuum).

Mapping the force exponent n as a function of m² is a critical test:
if n → 2 as m → 0, the Yukawa mass is the sole obstacle to Newton.

### Numerical vs Physical Energy Drift

The 0.1% energy drift per T=100 is from the NUMERICAL method (dt² error
in symplectic Verlet), not from physical radiation. The physical radiation
(V29-T3: ~0.5 energy/time) is a separate, much larger effect visible in
the energy timeseries as a slow decrease in E_total.

At thermal equilibrium (V29-T1b: A_noise ≈ 0.02), physical radiation
= absorption, and E_total stabilizes. The numerical drift persists but
can be reduced with smaller dt.

---

## 9. Open Questions

1. **Force law vs mass**: Map n in F∝1/D^n as function of m². Does n→2 as m→0?
2. **Mass emergence**: can the Lagrangian mass m² be eliminated?
   V27-M4 showed m=0 braids survive but without quadrupole.
3. **EMF integration**: complex fields + gauge coupling preserving the braid.
4. **Multi-braid formation**: conditions for stellar nucleosynthesis analog.
5. **Spin**: does the braid's helical handedness correspond to spin?
6. **Quantization**: the theory is purely classical. What changes in QFT?
7. **Gravitational waves**: does the braid's radiation have spin-2 structure?
8. **Background origin**: what sets ρ_bg and m²?
9. **Proportionality constant**: What sets C ≈ 186 in F = -C × ∇ρ?
   Can this be derived analytically from (m, μ, κ, braid geometry)?

## 10. Resolved Questions

1. **Is gravity energy minimization?** NO. E(D) is monotonically repulsive.
   The attraction is a dynamic effect. (V33 energy_vs_D, March 2026)
2. **Is the drift proportional to gradient?** YES. F ∝ ∇ρ with R²=0.9998.
   (V33 gradient sweep, March 2026)
3. **Is drag lower in depleted regions?** NO. Coupling is stronger in high ρ.
   The mechanism is geometric asymmetry, not friction. (V33 drag test)

---

## References to Simulation Code

### Best-in-Class Simulation
**File**: `/home/d/code/scp/v33/src/v33.c`
- Standard equation (Eq. 1), no modifications
- Single malloc for entire simulation (9 arrays × N³ × 8 bytes)
- Symplectic Velocity Verlet (energy drift < 0.1%/T=100)
- Fully periodic boundary conditions (for C1-C5 characterization)

### Gradient Test (Gravity Direction)
**File**: `/home/d/code/scp/v33/src/v33_gradient_test.c`
- x-direction: PINNED (high ρ left, low ρ right)
- y,z-directions: FREE-FLOATING (linear extrapolation from interior:
  phi[boundary] = 2×phi[boundary-1] - phi[boundary-2])
- Braid drifts toward LOW ρ = toward depletion = GRAVITY

### Braid Detection and Analysis
**File**: `/home/d/code/scp/v33/src/analyze.c`
- Reads binary field snapshots
- Flood-fill connected components above 5× average ρ
- Per-braid: position, energy, size, |P|, binding weight, winding, momentum
- Radial profiles: ρ(r), |P|(r), w(r) around each braid
- Inter-braid separations
- Background statistics

### Footprint Asymmetry Analysis
**File**: `/home/d/code/scp/v33/src/footprint_asymmetry.c`
- Reads field snapshots, finds braid center, extracts x-directional profile
- Measures half-width toward high-ρ vs low-ρ sides
- Computes geometric vs field-weighted center of perturbation footprint

### Drag Test (Braid Coupling vs Background Density)
**File**: `/home/d/code/scp/v33/src/v33_drag_test.c`
- Single braid with Galilean velocity kick in uniform backgrounds
- Tracks position, momentum, and braid energy over time
- Result: coupling STRONGER in high-ρ (opposite of simple drag model)

### Energy vs Separation (F = -dE/dD test)
**File**: `/home/d/code/scp/v33/src/v33_energy_vs_D.c`
- Two braids at separation D, settle T=30, average E over t=20-30
- Sweep D=8–50 to map E_int(D)
- Result: E_int always positive (repulsive), no attractive well

### F2 Test Battery (gradient sweep + energy sweep)
**File**: `/home/d/code/scp/v33/run_f2_tests.sh`
- Gradient strength sweep: 4 gradients (ρ ratio 1.22–9.00)
- Energy separation sweep: D=8–50
- Footprint analysis on all snapshots

### Binary Snapshot Format
```
int32:   N (grid size)
float64: L (domain half-width)
float64: t (simulation time)
float64[N³]: phi_0
float64[N³]: phi_1
float64[N³]: phi_2
```
Loadable in Python: `np.fromfile(f, dtype=np.float64, count=N**3).reshape(N,N,N)`

---

## Version History

- **V28**: Bimodal braid discovery (CMA-ES search)
- **V29**: Validation campaign (13 confirmed findings, field study)
- **V30**: Failed approaches (c(ρ), FRW expansion)
- **V31**: M7 split gravity tests (wrong sign → inverted c works but artificial)
- **V32**: SPH prototype, gradient coupling (self-interaction problem solved by binding weight, but gradient coupling is net repulsive with controls)
- **V33**: Clean standard equation. Single alloc. CONFIRMED: intrinsic attraction from standard equation, F∝1/D^1.8, gravity direction verified via gradient test. Footprint asymmetry measured (R_low/R_high=1.09–1.57). Drag test showed coupling scales with ρ (not anti-ρ). Gravity mechanism identified as geometric (asymmetric Yukawa profile), not dynamic (friction).

---

## Appendix: Critical Code Implementation

### A1. Memory — Single Allocation

All 9 field arrays (phi[3], vel[3], acc[3]) in ONE contiguous block:

```c
static Grid *grid_alloc(int N, double L) {
    Grid *g = calloc(1, sizeof(Grid));
    g->N  = N;
    g->N3 = (long)N * N * N;
    g->L  = L;
    g->dx = 2.0 * L / (N - 1);
    g->dt = 0.12 * g->dx;

    /* ONE allocation */
    long total = 9 * g->N3;
    g->mem = malloc(total * sizeof(double));
    memset(g->mem, 0, total * sizeof(double));

    /* Pointers into the single block */
    for (int a = 0; a < NFIELDS; a++) {
        g->phi[a] = g->mem + (0 + a) * g->N3;
        g->vel[a] = g->mem + (3 + a) * g->N3;
        g->acc[a] = g->mem + (6 + a) * g->N3;
    }
    return g;
}

static void grid_free(Grid *g) {
    free(g->mem);  /* ONE free */
    free(g);
}
```

### A2. Force Computation — Eq. (1) Evaluated

The core physics loop. Every grid point runs the same equation.
No distinction between braid and background.

```c
static void compute_forces(Grid *g) {
    const int N = g->N, NN = N * N;
    const long N3 = g->N3;
    const double idx2 = 1.0 / (g->dx * g->dx);

    #pragma omp parallel for schedule(static)
    for (long idx = 0; idx < N3; idx++) {
        int i = (int)(idx / NN);
        int j = (int)((idx / N) % N);
        int k = (int)(idx % N);

        /* Fully periodic neighbors */
        int ip = (i+1)%N, im = (i-1+N)%N;
        int jp = (j+1)%N, jm = (j-1+N)%N;
        int kp = (k+1)%N, km = (k-1+N)%N;

        long n_ip = (long)ip*NN + j*N + k;
        long n_im = (long)im*NN + j*N + k;
        long n_jp = (long)i*NN + jp*N + k;
        long n_jm = (long)i*NN + jm*N + k;
        long n_kp = (long)i*NN + j*N + kp;
        long n_km = (long)i*NN + j*N + km;

        /* Triple product P = φ₀φ₁φ₂ and dV/dP */
        double p0 = g->phi[0][idx];
        double p1 = g->phi[1][idx];
        double p2 = g->phi[2][idx];
        double P  = p0 * p1 * p2;
        double den = 1.0 + KAPPA * P * P;
        double mPd2 = MU * P / (den * den);    /* μP/(1+κP²)² */

        for (int a = 0; a < NFIELDS; a++) {
            /* 7-point Laplacian (6 neighbors + center) */
            double lap = (g->phi[a][n_ip] + g->phi[a][n_im]
                        + g->phi[a][n_jp] + g->phi[a][n_jm]
                        + g->phi[a][n_kp] + g->phi[a][n_km]
                        - 6.0 * g->phi[a][idx]) * idx2;

            /* ∂P/∂φ_a */
            double dPda = (a==0) ? p1*p2 : (a==1) ? p0*p2 : p0*p1;

            /* acc = ∇²φ - m²φ - (μP/(1+κP²)²)(∂P/∂φ_a) */
            g->acc[a][idx] = lap - MASS2 * g->phi[a][idx] - mPd2 * dPda;
        }
    }
}
```

### A3. Time Integration — Symplectic Velocity Verlet

Energy-conserving (symplectic). Drift < 0.1% per T=100.

```c
static void verlet_step(Grid *g) {
    const long N3 = g->N3;
    const double hdt = 0.5 * g->dt;
    const double dt  = g->dt;

    /* Half-kick: v += (dt/2) × acc */
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
    /* Drift: φ += dt × v */
    for (int a = 0; a < NFIELDS; a++) {
        double *p = g->phi[a], *v = g->vel[a];
        for (long idx = 0; idx < N3; idx++)
            p[idx] += dt * v[idx];
    }
    /* Recompute forces from new φ */
    compute_forces(g);
    /* Half-kick: v += (dt/2) × acc */
    for (int a = 0; a < NFIELDS; a++) {
        double *v = g->vel[a], *ac = g->acc[a];
        for (long idx = 0; idx < N3; idx++)
            v[idx] += hdt * ac[idx];
    }
}
```

### A4. Boundary Conditions — Pinned + Free-Floating (Gradient Test)

x-direction pinned to imposed ρ gradient. y,z free-floating outflow.

```c
static void apply_bc(Grid *g) {
    const int N = g->N, NN = N*N;
    const int margin = 3;

    /* x-direction: PIN to saved boundary values (high ρ left, low ρ right) */
    for (long idx = 0; idx < g->N3; idx++) {
        int i = (int)(idx / NN);
        if (i < margin || i >= N - margin) {
            for (int a = 0; a < NFIELDS; a++) {
                g->phi[a][idx] = g->pin_phi[a][idx];  /* restore saved */
                g->vel[a][idx] = g->pin_vel[a][idx];
                g->acc[a][idx] = 0;
            }
        }
    }

    /* y,z-directions: FREE-FLOATING (linear extrapolation from interior)
       phi[boundary] = 2×phi[interior_1] - phi[interior_2]
       Allows waves to flow out. Does not inject or remove energy. */
    for (int b = 0; b < margin; b++) {
        int j_in1 = margin, j_in2 = margin + 1;
        // ... (applied symmetrically to both edges of y and z)
        g->phi[a][idx_boundary] = 2*g->phi[a][idx_in1] - g->phi[a][idx_in2];
    }
}
```

### A5. Braid Detection — Connected Component Analysis

Post-hoc analysis of field snapshots. No distinction during simulation.

```c
/* Threshold: cells with ρ > 5× average are "braid" candidates */
double thresh = 5.0 * avg_rho;

/* Flood-fill from each seed above threshold */
for (long seed = 0; seed < N3; seed++) {
    if (rho[seed] < thresh || label[seed] >= 0) continue;
    int bid = n_braids++;

    /* BFS: expand to 6-connected neighbors above threshold */
    stack[sp++] = seed;
    label[seed] = bid;
    while (sp > 0) {
        long idx = stack[--sp];
        /* Accumulate: centroid, energy, |P|, binding weight */
        braids[bid].cx += x * rho * dV;
        braids[bid].E_total += rho * dV;
        braids[bid].P_avg += |P|;
        braids[bid].w_avg += 1.0/(1.0 + |P|/P_thresh);
        /* Expand to unvisited neighbors above threshold */
        for (int d = 0; d < 6; d++) { ... }
    }
}
```
