#!/usr/bin/env python3
"""
Launch the exact 7-run matrix requested:

- 3 single mechanisms (A, B, C)
- 3 pairwise combinations (A+B, A+C, B+C)
- 1 all-three (A+B+C)

All runs keep the base best 7-repulsion + connection-scaled repulsion.
Each writes a distinct rich snapshot named after the active mechanisms.

All seven are run in parallel using multiprocessing.
"""

import sys
import multiprocessing as mp
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import phase6_best_config_validation as val

# Reasonable non-zero strength for "on" mechanisms in this matrix
MECH_STRENGTH = 1.5

CONFIGS = [
    # (name, mech_A, mech_B, mech_C, extra_seed_offset)
    ("A_only",            MECH_STRENGTH, 0.0,          0.0,          0),
    ("B_only",            0.0,           MECH_STRENGTH, 0.0,          1),
    ("C_only",            0.0,           0.0,           MECH_STRENGTH, 2),
    ("A_plus_B",          MECH_STRENGTH, MECH_STRENGTH, 0.0,          3),
    ("A_plus_C",          MECH_STRENGTH, 0.0,           MECH_STRENGTH, 4),
    ("B_plus_C",          0.0,           MECH_STRENGTH, MECH_STRENGTH, 5),
    ("A_plus_B_plus_C",   MECH_STRENGTH, MECH_STRENGTH, MECH_STRENGTH, 6),
]

def run_one(cfg):
    name, a, b, c, seed_off = cfg
    print(f"\n=== Starting {name} ===")
    # Each run gets its own seed for reproducibility across the matrix
    result = val.run_best_config_validation(
        target_nodes=1600,          # slightly smaller for faster parallel runs
        growth_steps=35,
        seed=20260521 + seed_off,
        mech_A_w=a,
        mech_B_w=b,
        mech_C_w=c,
    )
    # Rename the snapshot so we don't overwrite each other
    old_path = SCRIPT_DIR / "phase6_best_config_validation_snapshot.json"
    new_path = SCRIPT_DIR / f"mech_{name}_snapshot.json"
    if old_path.exists():
        old_path.rename(new_path)
        print(f"  -> wrote {new_path}")
    print(f"=== Finished {name} ===")
    return name, result

if __name__ == "__main__":
    print("Launching 7 parallel mechanism experiments...")
    with mp.Pool(processes=7) as pool:
        results = pool.map(run_one, CONFIGS)

    print("\nAll seven runs completed.")
    for name, _ in results:
        print(f"  - mech_{name}_snapshot.json")