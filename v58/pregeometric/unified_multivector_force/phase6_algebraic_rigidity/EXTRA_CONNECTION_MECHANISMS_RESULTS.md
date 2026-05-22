# Extra Connection Mechanisms: Particle as an Additional Dimension of Connection

**Date**: 2026-05-21  
**Status**: Experimental results from the 7-run matrix + focused backprop on mechanism weights (Phase 6 algebraic rigidity program)

---

## 1. Conceptual Motivation

In the pre-geometric relational model, the living candidate produces self-stabilizing attractors whose effective dimensionality for internal observers has historically been limited to ~1.0–1.8. The long-standing puzzle is why localized, stable, high-density excitations (“particles”) form and persist, and why they appear to source an effective gravitational-like attraction.

One hypothesis explored in this phase is that a particle is not merely a dense excitation within the existing three-plane protected bivector algebra (the structure that gives rise to effective 3D geometry). Instead, a particle may involve an **additional layer of protected connection** — an “extra dimension of connection” — whose energetic or algebraic cost is strongly modulated by local density. In this view:

- Gravity is not a separate force but the effective consequence of density acting as a *scaler* on the cost of maintaining this extra connection.
- At low density the extra connection is expensive (or impossible) to sustain → vacuum remains “empty.”
- At high density the cost drops sharply → regions that have already formed particles make it cheaper for nearby regions to do the same → apparent attraction.

This framing naturally suggests that the algebraic rigidity and pressure terms we have been evolving should be extended with explicit density-dependent mechanisms that reward or stabilize additional protected structure only inside dense regions.

---

## 2. The Three Mechanisms

Three independent, toggleable pressure features were implemented and tested. Each can be activated by a positive weight in the linear combination that produces the per-node pressure scalar.

### Mechanism A — Extra Protected Connection Stabilized by Density
(`feat_mech_A_extra_conn_stab`, weight `MECH_A_EXTRA_CONN_STAB_W`)

**Implementation**:
```python
if rho > 0.25:
    return -biv_mag * (rho - 0.25)
```
where `biv_mag = sum(|b_i|)` for the three protected bivector components.

**Physics interpretation**: Nodes that already carry substantial protected bivector magnitude receive a *negative* pressure contribution (i.e., they become more attractive as parents for new nodes) once local density exceeds a soft threshold. This creates a positive feedback: dense regions that have “paid” to maintain protected structure become preferred sites for further growth. Conceptually this is the “extra connection” that a particle folds in; the density acts as the scaler that makes the connection cheap.

### Mechanism B — Density-Dependent Rigidity / Protection Cost
(`feat_mech_B_density_rigidity`, weight `MECH_B_DENSITY_RIGIDITY_W`)

**Implementation**: The existing compaction + expansion resistance signals (deviation from protected rank 3, weighted by bivector strength) are multiplied by a density-dependent factor that *reduces* the effective cost at high ρ:

```python
density_factor = 1.0 / (1.0 + 4.0 * max(0.0, rho - 0.2))
effective_cost = base_cost * density_factor
```

**Physics interpretation**: Maintaining algebraic rigidity (the closure that favors exactly three protected directions) becomes cheaper inside dense regions. This is a direct “gravity as scaler” effect on the protection cost itself.

### Mechanism C — Density-Dependent Algebraic Closure Preference
(`feat_mech_C_density_closure`, weight `MECH_C_DENSITY_CLOSURE_W`)

**Implementation**: At high density an additional negative pressure term rewards high protected bivector magnitude, effectively shifting the preferred closure condition:

```python
if rho > 0.2:
    return -biv_mag * (rho - 0.2) * 0.8
```

**Physics interpretation**: The algebraic target itself (what counts as a “good” protected configuration) becomes density-dependent. High-density regions are actively encouraged to adopt configurations that look like an extra connection has been folded in.

---

## 3. Experimental Design

All experiments were performed on top of the strongest base configuration discovered so far:

- The best 7-repulsion channels (direct-children, protected-children, ancestor-overlap, sibling-sharing, etc.) with boosted weights on the top performers.
- The connection-scaled repulsion term (`feat_connection_scaled_repulsion`) that applies ~10× repulsion to currently unconnected nodes and relaxes toward 1× as they gain children.

Two complementary protocols were used:

1. **7-run matrix** (coarse combinatorial test)  
   Each of A, B, C turned on at a fixed non-zero strength (1.5) in all combinations:
   - A only, B only, C only
   - A+B, A+C, B+C
   - A+B+C  
   All runs used identical seeds, ~1591 nodes, 35 growth steps. Late-time d_eff = mean of last 8 steps.

2. **Focused backprop (Differential Evolution)**  
   The three mechanism weights + the connection-scaled repulsion coefficient were treated as free parameters and optimized unconstrained (bounds [0,4] for mechanisms, [0.5,5] for scaling).  
   Objective: late-time d_eff + small stability bonus.  
   Evaluations performed at ~1500 nodes for speed.  
   Final discovered weights were then validated in a full 40-step / 1816-node run.

---

## 4. Results

### 4.1 7-Run Matrix (coarse combinatorial screen)

| Active Mechanisms     | Final d_eff | Late-time (last 8) | Peak d_eff | Comment |
|-----------------------|-------------|--------------------|------------|---------|
| A + B                 | **2.428**   | 2.268              | 2.428      | Highest final |
| B + C                 | **2.430**   | 2.276              | 2.430      | Very close 2nd |
| A + B + C             | 2.256       | **2.277**          | **2.476**  | Highest peak |
| A only                | 2.298       | 2.253              | 2.434      | Strongest singleton |
| A + C                 | 2.211       | 2.197              | 2.395      | Weaker |
| C only                | 2.156       | 2.139              | 2.339      | Modest |
| B only                | 1.964       | 2.125              | 2.376      | Clearly weakest |

**Key observations from the matrix**:
- Any combination containing **A** performed well.
- **B alone** was the weakest; B only became competitive when paired with A or C.
- The full triple produced the single highest instantaneous peak but did not give the best late-time stability (possible interference at equal strengths).

### 4.2 Backprop-Optimized Weights

The evolutionary search converged on:

| Parameter                        | Value   |
|----------------------------------|---------|
| `MECH_A_EXTRA_CONN_STAB_W`       | **3.8405** |
| `MECH_B_DENSITY_RIGIDITY_W`      | 1.4353  |
| `MECH_C_DENSITY_CLOSURE_W`       | 2.1873  |
| `CONNECTION_SCALED_REPULSION_W`  | 2.1335  |

Final 40-step validation (N=1816) with these weights:

- Late-time d_eff (mean of last 10 steps): **2.3045**
- Final d_eff at N=1816: **2.492**

This is the highest final effective dimensionality recorded in the entire Phase 6 series to date.

---

## 5. Interpretation and Current State

The optimizer’s strong preference for a large positive `MECH_A_EXTRA_CONN_STAB_W` (3.84) directly supports the hypothesis that an “extra dimension of connection” whose cost is lowered by local density is a useful inductive bias for the model. The density-dependent reward for high protected bivector magnitude inside dense regions creates exactly the positive feedback loop needed for stable particle-like objects that then source further growth.

Mechanism C (density-dependent shift in preferred closure) also received solid weight, suggesting that the algebraic target itself benefits from being density-modulated. Mechanism B (softening the rigidity cost) is helpful scaffolding but less critical once A is present.

With the current best configuration (strong 7-repulsion + connection-scaled 10×→1× repulsion + the optimized A+B+C weights), we are now consistently achieving:

- Late-time effective dimensionality in the **2.25 – 2.35** range
- Peaks and final values up to **2.492** on graphs of ~1800 nodes

This represents a clear advance beyond the ~1.0–1.8 plateau that characterized earlier phases and the ~2.1–2.4 ceiling obtained from repulsion alone.

---

## 6. Files and Reproducibility

All seven matrix runs and the final optimized validation wrote rich snapshots:
- `mech_A_only_snapshot.json`
- `mech_B_only_snapshot.json`
- ...
- `mech_A_plus_B_plus_C_snapshot.json`
- `phase6_best_config_validation_snapshot.json` (final optimized run)

The backprop script that produced the weights is `phase6_backprop_mech_weights.py`.

---

## 7. Local Causal Versions (Finite Propagation Speed)

After the first (partially global) Inside Space experiments, we implemented **strictly local** versions of both philosophies to respect causal locality and finite propagation speed — no node or core has global knowledge of the total inside volume.

### Local Version 1 (Reward + Localized Collapse, Causal)

- `feat_local_inside_volume_reward_v1` + `feat_local_volume_collapse_v1`
- Inside volume for a core is computed only from nodes causally connected to it within the recent layer window (its local light-cone).
- Reward for creating/protecting inside space is local.
- Collapse pressure on high-MECH_A cores is proportional only to the inside volume that has causal paths back to that specific core.

### Local Version 2 (Negative Pressure Inside + Boundary Collapse, Causal)

- `feat_local_inside_negative_pressure_v2` + `feat_local_boundary_collapse_v2`
- Nodes that are causally “inside” a core’s recent neighborhood receive negative pressure (interior stabilized after signal arrives).
- Collapse on the boundary remains local to the volume the core actually experiences.

Both use the helper `_local_enclosed_score_for_core(core, nodes)` which only looks at recent nodes that have the core (or its immediate high-MECH_A neighbors) as parents and have low outward connectivity.

### First Local Results (strength 1.2 on top of best base)

- Local V1: Final d_eff ≈ 2.257 (late ~2.26)
- Local V2: Final d_eff ≈ 2.225 (late ~2.20)

Both remained stable in the ~2.2 range. Locality did not destroy the effect.

### Strength Sweep Launched

Small manual sweep on the local versions:
- V1 at 0.8 and 2.0
- V2 at 1.5

These runs are executing in parallel.

### Ongoing Backprop

The focused backprop script (`phase6_backprop_mech_weights.py`) has been extended to a 10-dimensional search that now includes the four local Inside Space weights. It is running unconstrained and will report optimal strengths for the full set (MECH A/B/C, connection scaling, both chiral biases, and the four inside-space parameters).

---

*Document updated 2026-05-21. Local causal Inside Space versions added and first results + strength sweep + extended backprop documented.*