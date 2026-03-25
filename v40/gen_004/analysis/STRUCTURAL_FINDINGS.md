# Gen 4 Structural Analysis: Inner-Outer Field Characterization

## Data Source

Per-cluster radial profiles from `cluster_profile.c`, run on all three T=200
survivors (S20, CB15, UDD_R4). Each cluster's volume is divided into inner third,
middle third, and outer third by radius from the cluster centroid.

## Universal Radial Structure of Bound Clusters

All surviving clusters (>1000 voxels, E_pot < -0.1) share this profile:

```
        inner        mid         outer
ρ:      HIGH   →    medium   →   LOW      (outer/inner = 0.09-0.22)
|P|:    HIGH   →    medium   →   VERY LOW (outer/inner = 0.02-0.11)
θ_rms:  HIGH   →    medium   →   LOW      (outer/inner = 0.38-0.68)
|v|:    LOW    →    MEDIUM   →   HIGH     (outer/inner = 1.5-2.1)
|a|:    HIGH   →    medium   →   LOW      (outer/inner = 0.30-0.42)
E_pot:  HIGH   →    PEAKED   →   ~ZERO    (outer/inner = 0.01-0.12)
E_kin:  LOW    →    PEAKED   →   MODERATE (mid > inner AND mid > outer)
```

### Key finding: The "breathing shell" structure

Stable clusters have a THREE-LAYER architecture:

1. **Core** (inner third): High ρ, high |P|, high θ, LOW velocity, HIGH forces.
   The core is the binding engine — V(P) coupling is strongest here. The core
   is nearly stationary (|v| = 0.17-0.69) despite experiencing the strongest
   forces (|a| = 1.7-2.8). This means the forces are nearly BALANCED — the
   core is in dynamic equilibrium.

2. **Breathing shell** (middle third): Peak E_kin and E_pot, moderate everything.
   This is where energy concentrates during the breathing cycle. The shell
   oscillates between expansion and contraction, carrying energy between the
   core and the outer region.

3. **Outer halo** (outer third): Low ρ, very low |P|, LOW θ, HIGH velocity.
   The halo is kinetically dominated — field energy is being radiated outward.
   The velocity is high but the density is low, so the kinetic energy density
   is moderate.

## Discriminating Stable vs Unstable Clusters

### θ gradient is the clearest discriminant

| Cluster type | θ inner→outer | Δθ sign | Interpretation |
|-------------|---------------|---------|----------------|
| **Stable** (grew, >10K vox) | 0.06→0.03 | **NEGATIVE** | θ confined to core |
| **Unstable** (shrank) | 0.02→0.05 | **POSITIVE** | θ escaping outward |

**When θ_rms increases from inner to outer (outer/inner > 1.0), the cluster
is losing its electromagnetic coupling.** The curl(φ) source is in the core,
but the theta waves propagate outward faster than the binding can contain them.
This is the V39 result at the per-cluster level: theta radiation extracts energy.

### Velocity structure discriminates breathing vs dispersal

| Cluster type | |v| inner→outer | Δv sign | Interpretation |
|-------------|-----------------|---------|----------------|
| **Stable** | 0.17→0.39 | **POSITIVE** | Coherent breathing (calm core, active shell) |
| **Unstable** | 1.2→1.2 | **~ZERO** | Random motion (no organized dynamics) |

Stable clusters have a velocity GRADIENT: the core is calm, the shell moves.
This is coherent breathing — the structure that sustains the braid. Unstable
clusters have uniform velocity — they're just turbulent debris with no
organized motion.

### |P| concentration ratio

| Cluster type | |P| inner | |P| outer | Ratio |
|-------------|-----------|-----------|-------|
| **Stable** | 0.14-0.43 | 0.005-0.018 | **10-100×** |
| **Unstable** | 0.004-0.006 | 0.003-0.003 | **~1-2×** |

Stable clusters concentrate triple product in the core (ratio > 10). Unstable
clusters have nearly uniform (low) |P| — no concentrated binding region.

### Energy distribution

| Cluster type | E_pot location | E_kin location |
|-------------|---------------|---------------|
| **Stable** | **Core-dominated** (>80% in inner+mid) | **Mid-shell peaked** |
| **Unstable** | Uniformly weak | Uniformly weak |

The E_kin peaking in the mid-shell is a signature of organized breathing
dynamics. The potential energy is in the core (binding), the kinetic energy
is in the shell (oscillation). This separation of roles is what makes the
structure self-sustaining.

## Quantitative Thresholds

From the data across all three structures and all clusters:

| Metric | Stable threshold | Unstable indicator |
|--------|-----------------|-------------------|
| ρ inner | > 0.5 | < 0.1 |
| |P| inner | > 0.1 | < 0.01 |
| θ outer/inner | < 0.7 | > 1.0 |
| |v| outer/inner | > 1.3 | ≈ 1.0 |
| |a| inner | > 1.5 | < 1.0 |
| E_pot inner/total | > 0.5 | < 0.2 |
| Cluster volume | > 10,000 voxels | < 1,000 voxels |

## Per-Structure Summary

### S20 (braid 2.0×) — 3 clusters at T=200

Two large clusters (79K, 85K voxels) with nearly identical profiles:
- ρ: 1.6 inner → 0.18 outer (9× ratio)
- |P|: 0.41 inner → 0.017 outer (24× ratio)
- |v|: 0.18 inner → 0.39 outer (2.1× — strong breathing)
- θ: 0.065 inner → 0.030 outer (0.46 ratio — theta confined)
- **Stable architecture confirmed**

One small cluster (3.7K voxels):
- ρ: 0.23 → 0.11 (2× — weak gradient)
- |P|: 0.020 → 0.007 (3× — moderate)
- |v|: 1.41 → 1.28 (0.91 — no breathing, turbulent)
- θ: 0.032 → 0.053 (1.64× — **theta escaping**)
- **Unstable — will likely disperse**

### CB15 (counter-braid 1.5×) — 6 clusters at T=200

One dominant cluster (110K voxels, aspect=2.93 — elongated):
- ρ: 0.58 → 0.13 (4.5×)
- |v|: 0.83 → 1.24 (1.5× — breathing)
- θ: 0.041 → 0.034 (0.83 — theta mostly confined)
- E_kin strongly mid-peaked (62 vs 16 vs 47) — active breathing shell
- **Moderately stable — the elongation (aspect 2.93) is a concern**

Two medium clusters (48K, 27K) — both show stable architecture.

Three small clusters (<6K) — two show θ increasing outward (unstable signature).

### UDD_R4 (3-braid UDD) — 8 clusters at T=100

Main cluster (73K, aspect=2.57):
- ρ: 0.84 → 0.10 (8.4×)
- |P|: 0.16 → 0.007 (23× — very concentrated)
- |v|: 0.69 → 0.26 (0.38 — velocity DECREASING outward — NOT typical breathing)
- θ: 0.051 → 0.020 (0.38 — theta very well confined)
- **Strong binding but unusual velocity profile — the 3D structure creates
  a different dynamic than the 1D braid breathing**

Three secondary clusters (25K-30K each) — all show:
- Strong ρ and |P| gradients (stable)
- θ well-confined (outer/inner < 0.62)
- |v| decreasing outward (0.32-0.43 ratio) — **contracting, not breathing**

**Key UDD_R4 insight**: The 3-braid composite has CONTRACTING outer shells
(|v| decreases outward) unlike the single-braid structures (|v| increases).
The three perpendicular braids create an inward-directed flow pattern.
This is a fundamentally different stability mechanism — **gravitational
contraction rather than breathing oscillation**.

## Implications for Gen 4 Step 2 (Modifications)

1. **Target the θ confinement**: Modifications that reduce θ outer/inner ratio
   will stabilize clusters. Adding structure that absorbs outgoing theta waves
   (a "theta mirror") could prevent the radiation loss.

2. **Preserve the breathing shell**: Don't add high-amplitude perturbations to
   the mid-shell region — this is where E_kin peaks and the breathing dynamics
   live. Perturbations should target the core or the far outer halo.

3. **UDD_R4's contraction mechanism is novel**: The inward flow in the 3-braid
   structure is different from the braid's breathing. This could be the basis
   for a self-binding composite particle that doesn't rely on oscillation.

4. **Merge the small unstable clusters**: The small clusters in CB15 and UDD_R4
   are turbulent debris. If they could be absorbed into the main cluster
   (e.g., by initializing them closer to the core), the structure might
   consolidate into fewer, more stable units.

5. **Critical |P| threshold**: Clusters with inner |P| > 0.1 survive.
   Below 0.01, they disperse. The modification should ensure the added
   structure contributes to |P| in the core region.
