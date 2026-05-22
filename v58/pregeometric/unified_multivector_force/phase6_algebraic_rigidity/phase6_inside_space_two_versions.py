#!/usr/bin/env python3
"""
Launch the two "Inside Space" mechanisms independently in parallel.

Version 1: Explicit inside-volume reward + localized collapse on cores
Version 2: Negative pressure inside enclosed regions + collapse on boundary

Both use the current best base (7-repulsion + connection-scaled + optimized MECH A/B/C + direct chiral).

Each writes its own snapshot:
- inside_space_v1_snapshot.json
- inside_space_v2_snapshot.json
"""

import multiprocessing as mp
import phase6_best_config_validation as val

BEST_WEIGHTS = {
    "MECH_A": 3.8405,
    "MECH_B": 1.4353,
    "MECH_C": 2.1873,
    "CONNECTION_SCALED": 2.1335,
    "CHIRAL_DIRECT": 1.5,
}

INSIDE_STRENGTH = 1.2   # reasonable starting strength for the new mechanisms

def _run_one(func):
    """Helper to avoid pickling a lambda."""
    return func()


def run_v1():
    print("\n=== Launching Inside Space Version 1 (reward + collapse) ===")
    val.MECH_A_EXTRA_CONN_STAB_W = BEST_WEIGHTS["MECH_A"]
    val.MECH_B_DENSITY_RIGIDITY_W = BEST_WEIGHTS["MECH_B"]
    val.MECH_C_DENSITY_CLOSURE_W = BEST_WEIGHTS["MECH_C"]
    val.CONNECTION_SCALED_REPULSION_W = BEST_WEIGHTS["CONNECTION_SCALED"]
    val.CHIRAL_HELICAL_BIAS_W = BEST_WEIGHTS["CHIRAL_DIRECT"]
    val.CHIRAL_HELICAL_BIAS_INVERSE_W = 0.0

    # Version 1 active
    val.INSIDE_VOLUME_REWARD_V1_W = INSIDE_STRENGTH
    val.VOLUME_COLLAPSE_V1_W = INSIDE_STRENGTH
    val.INSIDE_NEGATIVE_PRESSURE_V2_W = 0.0
    val.BOUNDARY_COLLAPSE_V2_W = 0.0

    result = val.run_best_config_validation(
        target_nodes=1700,
        growth_steps=36,
        seed=20260620,
    )

    import os
    old = "phase6_best_config_validation_snapshot.json"
    new = "inside_space_v1_snapshot.json"
    if os.path.exists(old):
        os.rename(old, new)
        print(f"  -> wrote {new}")
    print("=== Finished V1 ===\n")
    return "V1", result


def run_v2():
    print("\n=== Launching Inside Space Version 2 (negative inside + boundary collapse) ===")
    val.MECH_A_EXTRA_CONN_STAB_W = BEST_WEIGHTS["MECH_A"]
    val.MECH_B_DENSITY_RIGIDITY_W = BEST_WEIGHTS["MECH_B"]
    val.MECH_C_DENSITY_CLOSURE_W = BEST_WEIGHTS["MECH_C"]
    val.CONNECTION_SCALED_REPULSION_W = BEST_WEIGHTS["CONNECTION_SCALED"]
    val.CHIRAL_HELICAL_BIAS_W = BEST_WEIGHTS["CHIRAL_DIRECT"]
    val.CHIRAL_HELICAL_BIAS_INVERSE_W = 0.0

    # Version 2 active
    val.INSIDE_VOLUME_REWARD_V1_W = 0.0
    val.VOLUME_COLLAPSE_V1_W = 0.0
    val.INSIDE_NEGATIVE_PRESSURE_V2_W = INSIDE_STRENGTH
    val.BOUNDARY_COLLAPSE_V2_W = INSIDE_STRENGTH

    result = val.run_best_config_validation(
        target_nodes=1700,
        growth_steps=36,
        seed=20260621,
    )

    import os
    old = "phase6_best_config_validation_snapshot.json"
    new = "inside_space_v2_snapshot.json"
    if os.path.exists(old):
        os.rename(old, new)
        print(f"  -> wrote {new}")
    print("=== Finished V2 ===\n")
    return "V2", result


if __name__ == "__main__":
    print("Spinning up both Inside Space versions in parallel...")
    with mp.Pool(2) as pool:
        results = pool.map(_run_one, [run_v1, run_v2])

    print("\nBoth Inside Space experiments completed.")
    for name, _ in results:
        print(f"  - inside_space_{name.lower()}_snapshot.json")