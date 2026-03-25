# V41 Phase Confinement: Detailed Analysis Results

## Frequency Analysis

### UUD (Proton): Clear breathing mode

| Observable | Period | Frequency | Notes |
|------------|--------|-----------|-------|
| E_pot | **90 t** | 0.011 | Strong coherent oscillation |
| P_int | 90 t | 0.011 | In phase with E_pot |
| phi_max | 90 t | 0.011 | In phase |
| theta_rms | — | — | Monotonically declining, no oscillation |

The proton has a **coherent breathing mode at period 90 time units**. All
binding observables (E_pot, P_int, phi_max) oscillate in phase. The E_pot
swings between -167 (binding peak) and -52 (breathing minimum). This is the
same breathing pattern seen in single braids but at a longer period (single
braid ~50 t vs composite ~90 t) — the composite is larger and breathes slower.

Theta_rms monotonically declines from 0.011 to 0.007 — the theta field is
slowly draining as expected, but not oscillating.

### UDD (Neutron): No detectable breathing period

| Observable | Period | Frequency | Notes |
|------------|--------|-----------|-------|
| E_pot | **0** (none) | — | Irregular oscillation |
| P_int | 0 | — | No periodicity |
| phi_max | 0 | — | Irregular |
| theta_rms | — | — | Declining |

The neutron has **NO coherent breathing mode**. The autocorrelation finds no
periodic signal. The E_pot fluctuates between -38 and -90 irregularly. This is
consistent with the zero-net-theta prediction: without a coherent electromagnetic
coupling (θ cancels), the three braids oscillate independently rather than
coherently. Each braid breathes on its own timescale, creating an irregular
composite signal.

**This is a fundamental UUD/UDD difference**: The proton breathes coherently
(all three braids synchronized), the neutron breathes incoherently (braids
desynchronized). The net θ coupling in UUD acts as a synchronization mechanism.

## Phase Relationships

### UUD (Proton): Two-group phase locking

At T=200, the 6 clusters show phases clustering into **two groups**:
- Group A: φ ≈ +2.35 (clusters 0, 1)
- Group B: φ ≈ -0.75 (clusters 2, 3)

Phase difference: 2.35 - (-0.75) = **3.10 ≈ π**. The two groups are in
ANTI-PHASE — when one group peaks, the other is at minimum. This is the
breathing mode: half the structure expands while the other half contracts.

The initial phases at t=0 were {-0.69, 1.49, 2.28, 2.31} — four clusters with
mixed phases. By T=200, they've self-organized into two anti-phase groups.
This is spontaneous phase locking driven by the V(P) coupling.

### UDD (Neutron): Phase convergence (loss of diversity)

At T=200, all 4 clusters have **nearly identical phases**: {2.39, 2.35, 2.35, 2.36}.
The phase diversity has COLLAPSED — all clusters converged to the same phase.

The initial phases at t=0 were {-0.69, 1.49, 2.31, 2.28} — similar to UUD.
But while UUD developed anti-phase structure, UDD lost all phase differences.

**Interpretation**: Without net θ coupling to maintain phase offsets, the braids
synchronize to the same phase and effectively become indistinguishable. This is
the LOSS OF COLOR CHARGE — the carrier phase offsets (which were supposed to
prevent merging) have been erased by the dynamics.

**Critical finding**: The phase confinement mechanism (P → 0 at triple overlap)
requires the carrier phases to remain DISTINCT. In UDD, the phases converged
→ the confinement mechanism is weakening. Given more time, the UDD structure
will merge.

In UUD, the anti-phase locking PRESERVES phase diversity → confinement persists.

## Multi-Structure Patterns

### UUD Cluster Evolution

```
t=0: 4 clusters → t=30: 9 (fragmentation) → t=120: 6 → t=200: 6 (stable)
```

Final structure at T=200 (6 clusters):

| Cluster | Voxels | P_peak | E_pot | θ_rms | Interpretation |
|---------|--------|--------|-------|-------|----------------|
| 1 | 42,114 | 0.331 | -53.1 | 0.057 | **Main core** (largest, strongest) |
| 3 | 21,385 | 0.365 | -27.6 | 0.041 | **Secondary core** |
| 2 | 17,684 | 0.266 | -15.9 | 0.045 | **Tertiary structure** |
| 6 | 6,726 | 0.060 | -1.1 | 0.044 | Peripheral (weak) |
| 5 | 7,502 | 0.050 | -0.6 | 0.028 | Peripheral (weak) |
| 4 | 2,821 | 0.019 | -0.1 | 0.021 | Debris (dying) |

**Three distinct bound clusters** (1, 2, 3) with P_peak > 0.2 and significant
E_pot. These are the three confined braids. They're at separations of:
- Cl1 to Cl2: ~13 code units
- Cl1 to Cl3: ~22 code units
- Cl2 to Cl3: ~15 code units

The three main clusters maintain separation (> 10 code units) — **confinement
is working**. They haven't merged despite 200 time units of evolution.

### UDD Cluster Evolution

```
t=0: 4 → t=30: 12 (heavy fragmentation) → t=120: 4 → t=200: 4
```

Final structure at T=200 (4 clusters):

| Cluster | Voxels | P_peak | E_pot | θ_rms | Interpretation |
|---------|--------|--------|-------|-------|----------------|
| 2 | 35,962 | 0.392 | **-54.8** | 0.046 | **Dominant blob** |
| 3 | 14,864 | 0.222 | -11.6 | 0.041 | Secondary |
| 4 | 9,957 | 0.060 | -1.7 | 0.038 | Peripheral |
| 1 | 4,810 | 0.031 | -0.3 | 0.039 | Debris |

**One dominant cluster** (Cl2 at 36K voxels) contains most of the binding energy.
The structure has partially merged — instead of three distinct braids, there's
one large blob + one secondary + two peripherals. The confinement is FAILING
for UDD as the phase offsets converge.

### P_int Distribution Over Time

**UUD per-cluster P_int shows oscillating redistribution:**
```
t=0:   cl0=57.7  cl1=0.7   cl2=0.7    (one dominant cluster)
t=30:  cl0=9.5   cl1=0.1   cl2=108.0  (energy shifted to cl2)
t=120: cl0=102.8 cl1=10.8  cl2=8.2    (back to cl0)
t=200: cl0=44.5  cl1=14.2  cl2=24.0   (more evenly distributed)
```

The binding energy CIRCULATES between clusters — this is the breathing mode.
Energy flows from cluster to cluster and back, with period ~90 t. By T=200,
the distribution is more even (44:14:24 vs initial 58:1:1), suggesting the
structure is equilibrating.

**UDD per-cluster P_int shows consolidation:**
```
t=0:   cl0=57.7  cl1=0.7   cl2=0.0    (one dominant)
t=120: cl0=57.7  cl1=14.7  cl2=2.3    (still dominant, some redistribution)
t=200: cl0=1.0   cl1=47.5  cl2=10.6   (dominance shifted to cl1)
```

The binding energy moved from cl0 to cl1 — the structure is reorganizing as
the dominant blob absorbs material. This is the merging process.

## Summary of Findings

### What the phase confinement produces:

| Property | UUD (Proton) | UDD (Neutron) |
|----------|-------------|---------------|
| Breathing mode | **90 t period, coherent** | None (irregular) |
| Phase structure | **Anti-phase locking (π separation)** | Phase convergence (all same) |
| Cluster count at T=200 | 6 (3 main + 3 peripheral) | 4 (1 dominant + 3 weak) |
| Confinement | **3 distinct cores maintained** | Partially merged |
| P_int distribution | Circulating between clusters | Consolidating into one |
| θ coupling | Synchronizes breathing | Fails to synchronize |

### Key insight: Net θ drives synchronization and confinement

The UUD proton's net θ ≠ 0 creates a coupling that:
1. **Synchronizes** the three braids' breathing (period 90 t)
2. **Locks** them in anti-phase (π phase offset between groups)
3. **Maintains** cluster separation (confinement persists)

The UDD neutron's net θ = 0 fails to:
1. Synchronize → irregular breathing, no coherent period
2. Lock phases → phases converge to same value
3. Maintain separation → dominant blob absorbs material

### Implications for V42

The phase confinement works for UUD but is weakening for UDD. Before attempting
nuclear binding (V42), we need to:

1. **Strengthen UDD confinement**: The carrier phase offsets need a stronger
   enforcement mechanism. The current linear superposition allows the dynamics
   to erase phase differences. A potential approach: use different amplitude
   or tube radius for the D-chirality braid to maintain distinctness.

2. **Verify UUD at longer T**: Run UUD to T=500+ to confirm the anti-phase
   locking is truly stable (not just a transient).

3. **Measure inter-braid forces**: Use the acceleration mode (viewer key 6)
   to map the forces between the three confined braids. Is the force profile
   consistent with confinement (short-range repulsion + medium-range attraction)?

4. **Test charge response**: Place a UUD proton in a θ gradient (from another
   structure) and verify it responds to the "electromagnetic" field.
