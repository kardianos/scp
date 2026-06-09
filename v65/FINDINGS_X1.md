# v65 FINDINGS — X1: homeostatic self-tuning (mechanism B) WORKS

**Date**: 2026-06-05
**Goal**: test mechanism (B) of [`SELF_TUNING.md`](SELF_TUNING.md) — a homeostatic control
law that drives the coupling κ to the stability edge with **no external hand-tuning and no
kernel edit**. The law: run a short relaxation at the current κ, sense collapse from the
core density P_max, then **push κ up on collapse / drift κ down when stable** (the action's
pull, ∂E/∂κ>0). The fixed point should be the self-organized-critical edge.

---

## Result — CONVERGED to a common edge κ\* ≈ 1.46 (SOC attractor is real)

Autonomous controller run *on the GPU box* (`remote_x1.py` via `sim_exec`), 128³, settled
T=8, from two opposite starts:

| start | trajectory | parked band | κ\* |
|---|---|---|---|
| κ=4.0 (high, stable) | drifts down (P_max 0.12→0.0014) to κ≈1.48 → collapse → bounce → oscillate | 1.19–1.57 | **1.36** |
| κ=0.5 (low, collapsed) | climbs out (P_max 14.7→4.3→2.9) to κ≈1.86 stable → oscillate | 1.18–1.94 | **1.56** |

**Both basins converge to a common band centered κ\*≈1.46.** The homeostatic law drives κ to
the marginal-stability edge *by itself, from both directions* — the defining signature of a
self-organized-critical attractor. The transition is razor-sharp (P_max jumps from ~0.001
stable to ~3–4 at collapse across the edge); the residual ±0.3 oscillation is the expected
critical fluctuation (SOC parks *near* the edge, not at a fixed real number — the caveat
predicted in `SELF_TUNING.md` §6). κ\*≈1.46 is consistent with the E3-measured stability
threshold (collapse between κ=0.5 and κ=1.5).

**Interpretation**: the coupling we use by hand (κ=50) is *not* what an intrinsic
self-tuning dynamics selects — the dynamics select the **edge of stability**, κ_crit≈1.46,
for this seed/charge. The self-tuned value is an *output* of the mechanism, not an input.

---

## Two bugs found and fixed en route (recorded so they aren't repeated)

1. **Coarse grid suppresses collapse (CPU 48³ null).** The first attempt ran the controller
   on CPU at 48³/T=1.5. P_max was **flat (0.0147) across all κ** — even κ=0.03, which
   collapses violently at 128³ — so the sensor carried *zero* stability information and the
   controller drifted κ monotonically to 0 with no bounce. Diagnosis: the density runaway is
   a **small-scale (UV) phenomenon**; a coarse grid (dx=0.43) cannot resolve the
   concentrating core and artificially stabilizes everything. **The collapse sensor requires
   a fine grid (≥128³).** CPU at 128³ is too slow (~12 min/step), hence the on-GPU loop.
2. **Transient compression mislabeled as collapse (GPU T=2).** At T=2 the unequilibrated
   oscillon seed *transiently compresses* (P_max~1.3–2.4) before settling, so a P_max>1.5
   sensor fired on moderate κ. The high-start run then parked at a spurious κ≈6. **Fix**:
   settle the transient (T=8, where stable κ give P_max~0.001 and true collapse gives
   P_max>10), threshold P_max>2. With this, both starts converge (above).

---

## What X1 establishes, and what it does not

- **Establishes**: mechanism (B) is real — an intrinsic-style feedback law self-organizes κ
  to the stability edge, basin-independently. "Self-tune to a point of stability" is
  demonstrated. The control law is automated (no hand-tuning), though it is still an *outer
  loop* around the kernel (discrete-time controller), not yet a term in the action.
- **Does not**: (i) it is *feedback*, not *action-derived* — the κ-update is imposed, not a
  δS/δκ EOM (that is mechanism A / X2); (ii) it lands on a *band*, not a sharp constant
  (SOC criticality, not a prediction to many digits); (iii) κ_crit here depends on the seed
  (its conserved content) — the falsifiable κ\*(Q) relation needs the charge-controlled
  version. Next: mechanism (A) — minimize the action **at fixed conserved charge** — which
  derives the feedback from a conservation law and is the bridge to the intrinsic X2 kernel
  field. See `SELF_TUNING.md` §3A.

**Artifacts**: `remote_x1.py` (the controller), run log in this file's tables.
