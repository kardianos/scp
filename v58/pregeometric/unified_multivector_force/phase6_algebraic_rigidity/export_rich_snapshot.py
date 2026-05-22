#!/usr/bin/env python3
"""
export_rich_snapshot.py

Helper to take a graph produced by phase6_rigidity_model.py (or similar)
and export a Lean-friendly JSON that contains per-node protected bivector
data + computed rigidity features.

This is what the Lean side (Phase6Rigidity.lean) will consume to
actually measure "what our current system dims resolves to".
"""

from pathlib import Path
import json
import sys

# Adjust path to import from the phase6 model if needed
SCRIPT_DIR = Path(__file__).resolve().parent
PHASE6_DIR = SCRIPT_DIR
if str(PHASE6_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE6_DIR))

# We assume the graph is already built by running phase6_rigidity_model.py
# and we can import the functions we need, or we re-implement the minimal
# rigidity computation here for export purposes.

# For now, this script is a template that shows the expected JSON schema
# that Phase6Rigidity.lean will be written to consume.

# Expected output format (per node):
# {
#   "id": int,
#   "rho": float,
#   "protected": bool,
#   "layer": int,
#   "biv_protected": [float, float, float],   # the three protected bivector coeffs (e12, e13, e23)
#   "local_rank": int,                        # number of "active" protected directions
#   "rigidity_cost": float
# }

# The full snapshot will also contain the living candidate params and graph summary.

def compute_simple_rigidity_from_biv(biv: list[float]) -> tuple[int, float]:
    """
    Replicates the Python-side rigidity heuristic.
    Returns (local_rank, cost).
    """
    b = [abs(x) for x in biv]
    active = [x for x in b if x > 1e-4]
    r = len(active)

    if r == 3:
        base = 0.0
    elif r == 2:
        base = 1.5
    elif r == 1:
        base = 4.0
    else:
        base = 6.0

    if r >= 2:
        s = sum(active)
        balance = 1.0 - (max(active) - min(active)) / (s + 1e-9) if s > 0 else 0.0
        closure_bonus = 1.0 - balance
    else:
        closure_bonus = 2.0

    cost = base + 1.2 * closure_bonus
    return r, cost


def main():
    # Example: load a previous snapshot that has node data
    # For the first real use, run the phase6_rigidity_model with an extended exporter,
    # or patch it to also dump per-node biv + rank + cost.

    print("This is a helper template.")
    print("In practice you will call this (or integrate the logic) after running a Phase 6 graph evolution")
    print("to produce a rich JSON with per-node protected bivector data for Lean analysis.")


if __name__ == "__main__":
    main()