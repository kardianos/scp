#!/usr/bin/env python3
"""
Phase 6 Ablation: 7 distinct Torsion (spin) types

This script keeps the strong 7-repulsion configuration fixed and systematically
ablates the new 7 torsion pressure features (expanded from the original 4).

It compares:
- Old 4-torsion baseline (the previous best torsion setup)
- Various subsets and strength variations of the new 7 torsion types

Each condition uses 3 seeds. Late-time d_eff = avg of last 8 steps.

Run:
    python phase6_ablation_7_torsion.py
"""

import sys
import statistics
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import phase6_best_config_validation as val


def set_7_repulson_weights():
    """Lock the 7 repulsion weights to the values that performed well."""
    val.REPULSION_DIRECT_W   = 1.10
    val.REPULSION_SIBLING_W  = 0.85
    val.REPULSION_CONE_W     = 0.70
    val.REPULSION_LAYER_W    = 0.65
    val.REPULSION_PROT_W     = 1.25
    val.REPULSION_BIVEC_W    = 0.55
    val.REPULSION_ANCESTOR_W = 0.80


def run_torsion_ablation(name, t0, t1, t2, t_ortho, t_signed, t_div, t_mag, seeds, target=1400, steps=30):
    """Run one torsion ablation condition."""
    results = []
    for s in seeds:
        set_7_repulson_weights()

        # Set the 7 torsion weights for this condition
        val.TORSION_PLANE_0_W       = t0
        val.TORSION_PLANE_1_W       = t1
        val.TORSION_PLANE_2_W       = t2
        val.TORSION_ORTHOGONAL_W    = t_ortho
        val.TORSION_SIGNED_W        = t_signed
        val.TORSION_DIVERSITY_W     = t_div
        val.TORSION_MAGNITUDE_GROW_W = t_mag

        print(f"\n=== {name} | seed={s} ===")
        out = val.run_best_config_validation(
            target_nodes=target,
            growth_steps=steps,
            seed=s
        )
        hist = out["history"]
        final_deff = hist[-1]["d_eff_proxy"]
        k = min(8, len(hist))
        late_deff = statistics.mean(h["d_eff_proxy"] for h in hist[-k:])
        peak_deff = max(h["d_eff_proxy"] for h in hist)

        results.append({
            "seed": s,
            "final": final_deff,
            "late": round(late_deff, 3),
            "peak": peak_deff,
        })
        print(f"  -> final={final_deff:.3f}  late={late_deff:.3f}  peak={peak_deff:.3f}")

    late_vals = [r["late"] for r in results]
    return {
        "name": name,
        "mean_late": round(statistics.mean(late_vals), 3),
        "std_late": round(statistics.stdev(late_vals), 3) if len(late_vals) > 1 else 0.0,
        "mean_final": round(statistics.mean(r["final"] for r in results), 3),
    }


def main():
    seeds = [20260521, 20260522, 20260523]

    # Ablation conditions for the 7 torsion types (old 4-torsion as baseline reference)
    conditions = [
        # Old 4-torsion baseline (approximated with the strongest previous weights)
        ("Old 4-torsion baseline (strong ortho + 3 planes)", 0.95, 1.05, 0.90, 1.30, 0.0, 0.0, 0.0),

        # Full new 7-torsion
        ("Full 7-torsion (all active)", 0.95, 1.05, 0.90, 1.30, 0.75, 0.65, 0.85),

        # 7-torsion without the 3 new types (signed, diversity, magnitude)
        ("7-torsion minus new 3 (only original-style 4)", 0.95, 1.05, 0.90, 1.30, 0.0, 0.0, 0.0),

        # Only the 3 new torsion types (signed + diversity + magnitude) + strong ortho
        ("New 3 torsion types + ortho only", 0.0, 0.0, 0.0, 1.30, 0.75, 0.65, 0.85),

        # High weight on signed torque
        ("High signed torque emphasis", 0.8, 0.9, 0.8, 1.0, 1.6, 0.5, 0.7),

        # High diversity (anti-alignment)
        ("High orientation diversity emphasis", 0.8, 0.9, 0.8, 1.0, 0.5, 1.6, 0.7),

        # No orthogonal (the previous strongest)
        ("7-torsion without orthogonal", 0.95, 1.05, 0.90, 0.0, 0.75, 0.65, 0.85),

        # Minimal torsion (almost off)
        ("Minimal torsion (near zero)", 0.2, 0.2, 0.2, 0.3, 0.2, 0.2, 0.2),
    ]

    print("=" * 95)
    print("PHASE 6 ABLATION — 7 Torsion Types (with 7-repulsion locked)")
    print("Comparing old 4-torsion vs expanded 7-torsion family")
    print("Target ~1400 nodes | 3 seeds | Late-time = avg last 8 steps")
    print("=" * 95)

    all_results = []
    for name, t0, t1, t2, to, ts, td, tm in conditions:
        res = run_torsion_ablation(name, t0, t1, t2, to, ts, td, tm, seeds)
        all_results.append(res)

    # Summary
    print("\n" + "=" * 95)
    print("7-TORSION ABLATION SUMMARY (sorted by mean late-time d_eff)")
    print("=" * 95)
    print(f"{'Condition':<55} {'Mean Late':>12} {'±std':>8} {'Mean Final':>12}")
    print("-" * 95)

    for r in sorted(all_results, key=lambda x: -x["mean_late"]):
        print(f"{r['name']:<55} {r['mean_late']:>12.3f} {r['std_late']:>8.3f} {r['mean_final']:>12.3f}")

    print("-" * 95)
    print("\nKey questions this ablation answers:")
    print("  • Does expanding to 7 torsion types improve over the old 4?")
    print("  • Which of the 3 new torsion signals (signed, diversity, magnitude) are useful?")
    print("  • Is the orthogonal term still the most important, or do the new ones matter more with 7-repulsion?")
    print("  • Are any of the 7 actively hindering?")


if __name__ == "__main__":
    main()