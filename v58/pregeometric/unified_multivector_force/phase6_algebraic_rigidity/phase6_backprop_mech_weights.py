#!/usr/bin/env python3
"""
Focused backprop (Differential Evolution) on the three new mechanism weights
(A, B, C) + the connection-scaled repulsion strength.

Goal: discover good relative strengths for the "extra dimension of connection"
ideas after the 7-run matrix showed that combinations involving A (and to some
extent B) are promising.

Base configuration kept fixed:
- Best 7-repulsion channels (with their current good weights)
- Connection-scaled 10x→1x repulsion is always active (its weight is also optimized)

Search space (10 dimensions):
- MECH_A_EXTRA_CONN_STAB_W
- MECH_B_DENSITY_RIGIDITY_W
- MECH_C_DENSITY_CLOSURE_W
- CONNECTION_SCALED_REPULSION_W
- CHIRAL_HELICAL_BIAS_W          (direct, ∝ connections)
- CHIRAL_HELICAL_BIAS_INVERSE_W  (inverse, 1/∝ connections)
- INSIDE_VOLUME_REWARD_V1_W
- VOLUME_COLLAPSE_V1_W
- INSIDE_NEGATIVE_PRESSURE_V2_W
- BOUNDARY_COLLAPSE_V2_W

Each evaluation runs a ~1400-1600 node validation and scores on late-time d_eff
(mean of last 8 steps) + a small stability bonus.

Run:
    python phase6_backprop_mech_weights.py
"""

import sys
import random
import statistics
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import phase6_best_config_validation as val

# --- Fixed base (the best repulsion setup we have) ---
# These are left at the values that performed well in prior ablations.
def lock_best_repulsion():
    val.REPULSION_DIRECT_W   = 1.8
    val.REPULSION_PROT_W     = 1.8
    val.REPULSION_ANCESTOR_W = 1.4
    val.REPULSION_SIBLING_W  = 1.1
    val.REPULSION_CONE_W     = 0.7
    val.REPULSION_LAYER_W    = 0.6
    val.REPULSION_BIVEC_W    = 0.5
    # Keep outer pressure modest
    val.OUTER_UNCONNECTED_PRESSURE_W = 1.0

# Reasonable search bounds for the six parameters we care about
BOUNDS = [
    (0.0, 5.0),   # MECH_A
    (0.0, 4.0),   # MECH_B
    (0.0, 4.0),   # MECH_C
    (0.5, 5.0),   # CONNECTION_SCALED_REPULSION_W
    (0.0, 4.0),   # CHIRAL_HELICAL_BIAS_W (direct)
    (0.0, 3.0),   # CHIRAL_HELICAL_BIAS_INVERSE_W
    (0.0, 3.0),   # INSIDE_VOLUME_REWARD_V1_W
    (0.0, 3.0),   # VOLUME_COLLAPSE_V1_W
    (0.0, 3.0),   # INSIDE_NEGATIVE_PRESSURE_V2_W
    (0.0, 3.0),   # BOUNDARY_COLLAPSE_V2_W
]

def evaluate(weights, seed_base=20260530):
    """
    weights = [mech_A, mech_B, mech_C, conn_scaled, chiral_direct, chiral_inverse,
               inside_reward_v1, volume_collapse_v1, inside_neg_v2, boundary_collapse_v2]
    Returns negative score (because DE minimizes).
    """
    (mech_a, mech_b, mech_c, conn_w, chiral_direct, chiral_inverse,
     inside_reward_v1, volume_collapse_v1, inside_neg_v2, boundary_collapse_v2) = weights

    lock_best_repulsion()
    val.CONNECTION_SCALED_REPULSION_W = conn_w
    val.MECH_A_EXTRA_CONN_STAB_W = mech_a
    val.MECH_B_DENSITY_RIGIDITY_W = mech_b
    val.MECH_C_DENSITY_CLOSURE_W = mech_c
    val.CHIRAL_HELICAL_BIAS_W = chiral_direct
    val.CHIRAL_HELICAL_BIAS_INVERSE_W = chiral_inverse
    val.INSIDE_VOLUME_REWARD_V1_W = inside_reward_v1
    val.VOLUME_COLLAPSE_V1_W = volume_collapse_v1
    val.INSIDE_NEGATIVE_PRESSURE_V2_W = inside_neg_v2
    val.BOUNDARY_COLLAPSE_V2_W = boundary_collapse_v2

    # Use a different seed per evaluation for robustness
    seed = seed_base + random.randint(0, 10000)

    try:
        result = val.run_best_config_validation(
            target_nodes=1500,
            growth_steps=32,
            seed=seed,
        )
        hist = result.get("history", [])
        if len(hist) < 5:
            return 0.0   # bad run

        # Primary objective: late-time d_eff (mean of last 8)
        late_vals = [h["d_eff_proxy"] for h in hist[-8:]]
        late_mean = statistics.mean(late_vals)

        # Small bonus for stability (low variance in last 8)
        late_std = statistics.stdev(late_vals) if len(late_vals) > 1 else 0.0
        stability = max(0.0, 0.05 - late_std) * 2.0

        score = late_mean + stability
        return -score   # DE minimizes

    except Exception as e:
        print(f"  Evaluation failed: {e}")
        return 0.0

def main():
    print("Starting focused backprop on mechanism weights (A, B, C + connection scaling)...")
    print("Bounds:", BOUNDS)
    print("Population will be evaluated at ~1500 nodes for speed.\n")

    # Use a modest population size for speed
    result = differential_evolution(
        evaluate,
        BOUNDS,
        workers=4,            # parallel evaluations
        updating='deferred',
        popsize=12,
        mutation=0.7,
        recombination=0.8,
        tol=0.015,
        seed=20260530,
        disp=True,
        polish=False,         # we can polish later if needed
    )

    print("\n" + "="*70)
    print("BACKPROP FINISHED")
    print("="*70)
    print("Best weights found:")
    print(f"  MECH_A_EXTRA_CONN_STAB_W       = {result.x[0]:.4f}")
    print(f"  MECH_B_DENSITY_RIGIDITY_W      = {result.x[1]:.4f}")
    print(f"  MECH_C_DENSITY_CLOSURE_W       = {result.x[2]:.4f}")
    print(f"  CONNECTION_SCALED_REPULSION_W  = {result.x[3]:.4f}")
    print(f"  CHIRAL_HELICAL_BIAS_W          = {result.x[4]:.4f}")
    print(f"  CHIRAL_HELICAL_BIAS_INVERSE_W  = {result.x[5]:.4f}")
    print(f"  INSIDE_VOLUME_REWARD_V1_W      = {result.x[6]:.4f}")
    print(f"  VOLUME_COLLAPSE_V1_W           = {result.x[7]:.4f}")
    print(f"  INSIDE_NEGATIVE_PRESSURE_V2_W  = {result.x[8]:.4f}")
    print(f"  BOUNDARY_COLLAPSE_V2_W         = {result.x[9]:.4f}")
    print(f"\nBest (negative) score = {result.fun:.6f}")
    print(f"Equivalent late-time d_eff ≈ {-result.fun:.4f}")

    # One final full-size validation with the discovered weights
    print("\nRunning one final ~1800-node validation with the discovered weights...")
    lock_best_repulsion()
    val.CONNECTION_SCALED_REPULSION_W = result.x[3]
    val.MECH_A_EXTRA_CONN_STAB_W = result.x[0]
    val.MECH_B_DENSITY_RIGIDITY_W = result.x[1]
    val.MECH_C_DENSITY_CLOSURE_W = result.x[2]

    final = val.run_best_config_validation(
        target_nodes=1800,
        growth_steps=40,
        seed=20260599,
    )
    hist = final.get("history", [])
    if hist:
        late = statistics.mean(h["d_eff_proxy"] for h in hist[-10:])
        print(f"\nFinal validation late-time d_eff (last 10 steps) = {late:.4f}")
        print(f"Final d_eff at N={hist[-1]['N']} = {hist[-1]['d_eff_proxy']:.4f}")

if __name__ == "__main__":
    main()