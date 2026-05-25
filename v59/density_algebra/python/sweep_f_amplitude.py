#!/usr/bin/env python3
"""
Sweep over f_amplitude (F-content strength in mixed sources)
to see when/ how the different protection strategies differentiate.

This directly tests the "stacking vs isolation" hypothesis for L vs F protection technologies.
"""

import sys
import json
import os
from datetime import datetime
import numpy as np

sys.path.insert(0, '/home/d/code/scp/v59/density_algebra/python')
from improved_density_protection_scan import run_improved_scan

F_AMPLITUDES = [0.1, 0.25, 0.4, 0.55, 0.7, 0.9]
PROTECTIONS = ['unprotected_full', 'protected_L', 'protected_F', 'protected_LF']
LAMBDA = 0.005
N_STEPS = 22   # balance speed vs stability of metrics
RHO_CRIT = 2.5

def main():
    print("=" * 85)
    print("f_amplitude SWEEP — Mixed L+F Sources under Different Protection Strategies")
    print(f"lambda = {LAMBDA}, steps = {N_STEPS}")
    print("Goal: Find when F-content in the source makes pure-L protection worse than L+F stacking,")
    print("      and how F-only protection behaves as F-content increases.")
    print("=" * 85)

    all_results = []

    for f_amp in F_AMPLITUDES:
        print(f"\n=== f_amplitude = {f_amp:.2f} ===")
        row = {"f_amplitude": f_amp}
        for prot in PROTECTIONS:
            res = run_improved_scan(
                lambda_nl=LAMBDA,
                protection_name=prot,
                n_steps=N_STEPS,
                f_amplitude=f_amp,
                rho_crit=RHO_CRIT
            )
            row[prot] = {
                "peak_density": round(res["peak_density"], 5),
                "cross_term_pct": round(res["cross_term_pct"], 2),
                "comm_err_avg": round(res["comm_err_avg"], 4),
                "num_active": res["num_active"]
            }
            print(f"  {prot:18s}  peak_ρ={res['peak_density']:.5f}  cross={res['cross_term_pct']:.2f}%  comm={res['comm_err_avg']:.4f}")

        all_results.append(row)

    # Save structured report
    os.makedirs("/home/d/code/scp/v59/density_algebra/reports", exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = f"/home/d/code/scp/v59/density_algebra/reports/f_amplitude_sweep_{ts}.json"

    report = {
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "lambda": LAMBDA,
            "n_steps": N_STEPS,
            "rho_crit": RHO_CRIT,
            "f_amplitudes": F_AMPLITUDES,
            "protections": PROTECTIONS
        },
        "results": all_results,
        "interpretation_note": "Look for crossover points where protected_L peak_density drops relative to protected_LF as f_amplitude rises (evidence for stacking benefit). Also note when protected_F cross_term stays low while its density remains competitive."
    }

    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nReport saved to: {out_path}")

    # Quick summary table for the user
    print("\n" + "=" * 85)
    print("PEAK DENSITY vs f_amplitude (key observable for stacking benefit)")
    print("-" * 85)
    print(f"{'f_amp':>6} | {'unprot_full':>12} | {'prot_L':>10} | {'prot_F':>10} | {'prot_LF':>10}")
    print("-" * 85)
    for r in all_results:
        fa = r["f_amplitude"]
        print(f"{fa:6.2f} | {r['unprotected_full']['peak_density']:12.5f} | "
              f"{r['protected_L']['peak_density']:10.5f} | "
              f"{r['protected_F']['peak_density']:10.5f} | "
              f"{r['protected_LF']['peak_density']:10.5f}")
    print("=" * 85)

if __name__ == "__main__":
    main()
