#!/usr/bin/env python3
"""
Phase 6 Ablation Study: Multi-orientation Torsion + Internal Repulsion

Systematically turns components on/off and varies their strengths to determine
which parameters are operational, neutral, or hindering the high-d_eff regime
(~2.1–2.35) discovered in the long run with repulsion + torsion.

Each condition is run with 3 different seeds for basic statistical robustness.
Late-time d_eff is defined as the average of the last 8 steps (or all steps if fewer).

Run:
    python phase6_ablation_repulsion_torsion.py
"""

import sys
import statistics
from pathlib import Path

# Import the validation module so we can override its weights
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import phase6_best_config_validation as val


def run_condition(name, repulsion_w, t0, t1, t2, t_ortho, seeds, target=1400, steps=30):
    """Run one ablation condition across multiple seeds and return stats."""
    results = []
    for s in seeds:
        # Override the relevant weights (globals in the imported module)
        val.INTERNAL_REPULSION_WEIGHT = repulsion_w
        val.TORSION_ORIENT_0_WEIGHT = t0
        val.TORSION_ORIENT_1_WEIGHT = t1
        val.TORSION_ORIENT_2_WEIGHT = t2
        val.TORSION_ORTHO_WEIGHT = t_ortho
        val.EXTERNAL_PRESSURE_BOOST = 1.0   # keep neutral

        print(f"\n=== {name} | seed={s} ===")
        out = val.run_best_config_validation(
            target_nodes=target,
            growth_steps=steps,
            seed=s
        )
        hist = out["history"]
        final_deff = hist[-1]["d_eff_proxy"]
        # Late-time: average of last min(8, len) steps
        k = min(8, len(hist))
        late_deff = statistics.mean(h["d_eff_proxy"] for h in hist[-k:])
        peak_deff = max(h["d_eff_proxy"] for h in hist)

        results.append({
            "seed": s,
            "final": final_deff,
            "late": round(late_deff, 3),
            "peak": peak_deff,
            "N": out["final_N"]
        })
        print(f"  -> final={final_deff:.3f}  late_avg={late_deff:.3f}  peak={peak_deff:.3f}")

    # Aggregate
    late_vals = [r["late"] for r in results]
    mean_late = statistics.mean(late_vals)
    std_late = statistics.stdev(late_vals) if len(late_vals) > 1 else 0.0

    return {
        "name": name,
        "results": results,
        "mean_late": round(mean_late, 3),
        "std_late": round(std_late, 3),
        "mean_final": round(statistics.mean(r["final"] for r in results), 3),
    }


def main():
    seeds = [20260521, 20260522, 20260523]

    # Define the ablation matrix (most diagnostic contrasts first)
    conditions = [
        # (name, repulsion, t0, t1, t2, t_ortho)
        ("Full recipe (rep 1.9 + 4 torsion)", 1.9, 1.05, 1.15, 0.95, 1.40),
        ("Repulsion only (no torsion)", 1.9, 0.0, 0.0, 0.0, 0.0),
        ("Torsion only (no repulsion)", 0.0, 1.05, 1.15, 0.95, 1.40),
        ("Base pressure only (neither)", 0.0, 0.0, 0.0, 0.0, 0.0),
        ("Rep + ortho torsion only", 1.9, 0.0, 0.0, 0.0, 1.40),
        ("High repulsion (3.0 + torsion)", 3.0, 1.05, 1.15, 0.95, 1.40),
        ("Low repulsion (0.8 + torsion)", 0.8, 1.05, 1.15, 0.95, 1.40),
        ("No orthogonal term (rep + 3 torsion)", 1.9, 1.05, 1.15, 0.95, 0.0),
    ]

    print("=" * 90)
    print("PHASE 6 ABLATION — Internal Repulsion + Multi-orientation Torsion")
    print("Target ~1400 nodes | 3 seeds per condition | Late-time = avg of last 8 steps")
    print("=" * 90)

    all_results = []
    for name, r, t0, t1, t2, to in conditions:
        res = run_condition(name, r, t0, t1, t2, to, seeds)
        all_results.append(res)

    # Final summary table
    print("\n" + "=" * 90)
    print("ABLATION SUMMARY (sorted by mean late-time d_eff)")
    print("=" * 90)
    print(f"{'Condition':<45} {'Mean Late d_eff':>16} {'±std':>8} {'Mean Final':>12}")
    print("-" * 90)

    sorted_results = sorted(all_results, key=lambda x: -x["mean_late"])
    for r in sorted_results:
        print(f"{r['name']:<45} {r['mean_late']:>16.3f} {r['std_late']:>8.3f} {r['mean_final']:>12.3f}")

    print("-" * 90)
    print("\nInterpretation guidance:")
    print("  • Large drop when turning repulsion off → repulsion is the main driver.")
    print("  • Large drop when turning torsion off (while keeping repulsion) → torsion adds value.")
    print("  • Base pressure only should reproduce the old ~0.95–1.0 ceiling.")
    print("  • Comparing 'Rep + ortho only' vs full torsion tells if the three plane terms matter.")
    print("  • High vs low repulsion shows whether there is an optimum strength.")
    print("\nRaw per-seed data saved in memory above. Re-run with more seeds or larger N if needed.")


if __name__ == "__main__":
    main()