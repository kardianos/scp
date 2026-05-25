#!/usr/bin/env python3
"""
targeted_density_protection_scan.py

Targeted experiments for v59/density_algebra work.

Deliberately varies:
  - "Algebraic content / strength": the quadratic coefficient lambda_nl
    (and optionally mu_nl) in the living-candidate multivector equation.
  - "Protection strength": the protected_biv flag (bivector support restriction
    on the test excitation M_t, corresponding to the "protected chirality"
    mechanism validated in the v58 runs).

For each combination we measure:
  - Density achieved: proxy via average scalar part (s_part) or effective rho
    of the high-density lumps during the retarded evolution.
  - Force cleanliness: cross-term percentage (neutral vs charged), commutation
    / leakage error, and whether the run stays inside the known safe band.

Results are printed as a compact table and written to a JSON report that can
be directly referenced from the hypotheses and CONSTRAINTS_AND_TROUGHS.md.

This is the first concrete quantitative feedback loop from the Python side
of the unified multivector force experiments into the physical-motivation work.

Run (from anywhere):
    python -m v59.density_algebra.python.targeted_density_protection_scan
or simply:
    cd v59/density_algebra/python && python targeted_density_protection_scan.py
"""

import sys
import json
import os
from datetime import datetime
from typing import Dict, List, Any

# Make the existing v58 experiment code importable
V58_PYTHON = "/home/d/code/scp/v58/pregeometric/unified_multivector_force/python"
if V58_PYTHON not in sys.path:
    sys.path.insert(0, V58_PYTHON)

# Import the reusable simulation driver and helpers
try:
    from retarded_dynamic_scan import (
        run_retarded_simulation,
        f_winning,
    )
except ImportError as e:
    print("ERROR: Could not import from retarded_dynamic_scan.py")
    print(f"  Make sure the path {V58_PYTHON} is correct and the file exists.")
    print(f"  Original error: {e}")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Configuration of the targeted sweep
# ---------------------------------------------------------------------------

# Algebraic strength (quadratic coefficient λ in the living candidate)
LAMBDA_VALUES = [0.0, 0.003, 0.005, 0.008]   # 0.0 = pure linear baseline

# Protection strength
PROTECTION_LEVELS = [
    {"name": "unprotected", "protected_biv": False},
    {"name": "protected_biv", "protected_biv": True},
]

# Other fixed parameters for this first scan (can be varied later)
N_STEPS = 35
DT = 0.06

OUTPUT_DIR = "/home/d/code/scp/v59/density_algebra/reports"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Helper to extract a simple "density achieved" proxy from a run
# ---------------------------------------------------------------------------

def extract_density_proxy(result: Dict[str, float]) -> float:
    """
    Use the average |force| on the charged test particle as a rough proxy
    for the strength of the density-sourced interaction it feels.
    In a fuller version we would compute the actual scalar part of M~M
    (ρ_M) at multiple probe locations and return max/avg.
    For now the force magnitude is directly tied to the ambient density
    gradients in the model, so it serves as a monotonic proxy.
    """
    # run_retarded_simulation returns 'avg_force_c' among other keys
    return float(result.get("avg_force_c", result.get("avg_force", 0.0)))


def extract_cleanliness_metrics(result: Dict[str, float]) -> Dict[str, float]:
    return {
        "cross_term_pct": float(result.get("cross_pct", 0.0)),
        "comm_err": float(result.get("comm_err", result.get("comm_errs", [0.0])[-1] if "comm_errs" in result else 0.0)),
    }


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("TARGETED DENSITY vs PROTECTION STRENGTH SCAN")
    print("For v59/density_algebra hypotheses (forcing relationships)")
    print(f"Date: {datetime.now().isoformat(timespec='seconds')}")
    print("=" * 78)

    results: List[Dict[str, Any]] = []

    for lam in LAMBDA_VALUES:
        for prot in PROTECTION_LEVELS:
            print(f"\nRunning: lambda={lam:.4f}  protection={prot['name']}")
            res = run_retarded_simulation(
                lambda_nl=lam,
                use_retarded=True,
                protected_biv=prot["protected_biv"],
                n_steps=N_STEPS,
                dt=DT,
            )

            dens = extract_density_proxy(res)
            clean = extract_cleanliness_metrics(res)

            record = {
                "lambda_nl": lam,
                "protection": prot["name"],
                "protected_biv": prot["protected_biv"],
                "density_proxy": dens,
                "cross_term_pct": clean["cross_term_pct"],
                "comm_err": clean["comm_err"],
                "raw_stats": {k: float(v) for k, v in res.items() if isinstance(v, (int, float))},
            }
            results.append(record)

            print(f"  density_proxy = {dens:.5f}   "
                  f"cross_term = {clean['cross_term_pct']:.2f}%   "
                  f"comm_err = {clean['comm_err']:.4f}")

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("SUMMARY TABLE (density achieved vs force cleanliness)")
    print("-" * 78)
    print(f"{'lambda':>8}  {'protection':<14}  {'dens_proxy':>10}  {'cross_%':>8}  {'comm_err':>8}")
    print("-" * 78)
    for r in results:
        print(f"{r['lambda_nl']:8.4f}  {r['protection']:<14}  "
              f"{r['density_proxy']:10.5f}  {r['cross_term_pct']:8.2f}  {r['comm_err']:8.4f}")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # Write report for feeding back into the hypotheses
    # -----------------------------------------------------------------------
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join(OUTPUT_DIR, f"protection_density_scan_{timestamp}.json")

    report = {
        "timestamp": datetime.now().isoformat(),
        "description": "Targeted sweep of quadratic strength (lambda) and protected-chirality (bivector support restriction) measuring density proxy vs force-channel cleanliness. First quantitative data point for v59/density_algebra hypotheses.",
        "parameters": {
            "lambda_values": LAMBDA_VALUES,
            "protection_levels": [p["name"] for p in PROTECTION_LEVELS],
            "n_steps": N_STEPS,
            "dt": DT,
        },
        "results": results,
        "links": {
            "hypotheses": "../HYPOTHESES.md",
            "constraints": "../CONSTRAINTS_AND_TROUGHS.md",
            "source_script": "v58/pregeometric/unified_multivector_force/python/retarded_dynamic_scan.py",
        },
    }

    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nReport written to: {report_path}")
    print("These numbers can now be referenced from the hypotheses when discussing")
    print("the protection-budget vs density-well trade-off and the forcing power")
    print("of the protected-chirality mechanism.")

    return report_path


if __name__ == "__main__":
    main()