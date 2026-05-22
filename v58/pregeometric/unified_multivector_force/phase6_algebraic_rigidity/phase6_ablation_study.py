#!/usr/bin/env python3
"""
phase6_ablation_study.py

Small post-hoc ablation study on the best coefficients discovered by the
intensive Differential Evolution search.

We take the best coefficient vector and systematically zero out the strongest
terms (individually and in small groups) to measure how much each contributes
to the loss.

This gives a quantitative ranking of importance and helps decide which terms
are genuinely useful versus which were marginal.

The study is intentionally lightweight so it can be run quickly before the
full validation runs on the discrete simulator.
"""

import sys
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PHASE1_DIR = (SCRIPT_DIR.parent / "phase1_minimal_relational_models").resolve()
if str(PHASE1_DIR) not in sys.path:
    sys.path.insert(0, str(PHASE1_DIR))

from phase6_backprop import (
    ALL_TERMS,
    TERM_NAMES,
    evaluate_coefficients,
)

# Best coefficients from the intensive search (2026-05-20)
BEST_COEFFS = np.array([
    2.1491,   # term_rank_deviation
    1.2325,   # term_rank_abs
   -1.8352,   # term_closure_defect
    0.4442,   # term_balance
    0.6928,   # term_rho_rank_interaction
    0.7384,   # term_rho_closure_interaction
   -0.5153,   # term_biv_magnitude
    2.2414,   # term_rank_biv_interaction
])

def ablation_study(n_seeds: int = 1000, growth_steps: int = 100, perturbation_strength: float = 0.08):
    """
    Very intensive post-hoc ablation study with live progress writing to disk.

    High statistical power version:
    - n_seeds=1000 (or 30,000 if you have the patience)
    - growth_steps=100 (100x finer than original tiny studies)
    """
    from joblib import Parallel, delayed
    import numpy as np
    import json
    from datetime import datetime

    progress_file = SCRIPT_DIR / "phase6_ablation_live_progress.log"
    results_file = SCRIPT_DIR / "phase6_ablation_results_intensive.json"

    def log(msg):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {msg}"
        print(line, flush=True)
        with open(progress_file, "a") as f:
            f.write(line + "\n")

    def _eval_one(coeffs, seed):
        return evaluate_coefficients(
            coeffs,
            seeds=[seed],
            n_nodes=180,
            growth_steps=growth_steps,
            perturbation_strength=perturbation_strength
        )

    def evaluate_with_high_precision(coeffs, n_seeds, label=""):
        """Run many independent seeds in parallel, write partial progress every 200 seeds."""
        seeds = [20260519 + i for i in range(n_seeds)]
        all_losses = []

        chunk_size = 200
        for start in range(0, n_seeds, chunk_size):
            chunk = seeds[start:start + chunk_size]
            losses = Parallel(n_jobs=16, backend="loky", verbose=0)(
                delayed(_eval_one)(coeffs, s) for s in chunk
            )
            all_losses.extend(losses)

            current_mean = float(np.mean(all_losses))
            current_se = float(np.std(all_losses, ddof=1) / np.sqrt(len(all_losses)))
            log(f"  {label} — {len(all_losses)}/{n_seeds} seeds | mean={current_mean:.4f} ± {current_se:.4f}")

        losses = np.array(all_losses)
        mean_loss = float(np.mean(losses))
        stderr = float(np.std(losses, ddof=1) / np.sqrt(len(losses)))
        return mean_loss, stderr

    log("=== Starting VERY INTENSIVE ablation study ===")
    log(f"n_seeds per config      : {n_seeds}")
    log(f"growth_steps per eval   : {growth_steps}")
    log(f"perturbation_strength   : {perturbation_strength}")
    log(f"parallel workers        : 16")
    log("")

    # Baseline
    log("Evaluating BASELINE (full best vector)...")
    base_mean, base_se = evaluate_with_high_precision(BEST_COEFFS, n_seeds, label="BASELINE")
    log(f"Baseline: {base_mean:.4f} ± {base_se:.4f}\n")

    results = []

    # Individual term ablations
    for i, name in enumerate(TERM_NAMES):
        log(f"Evaluating ablation: ZERO {name}")
        ablated = BEST_COEFFS.copy()
        ablated[i] = 0.0
        mean_loss, se = evaluate_with_high_precision(ablated, n_seeds, label=f"ZERO {name}")
        delta = mean_loss - base_mean
        results.append({
            "term": name,
            "zeroed_indices": [i],
            "mean_loss": mean_loss,
            "stderr": se,
            "delta_from_baseline": delta
        })
        log(f"  → {mean_loss:.4f} ± {se:.4f}   (Δ = {delta:+.4f})\n")

    # Group ablation — Top 3 strongest
    strong_indices = [7, 0, 2]
    log(f"Evaluating ablation: ZERO Top 3 strongest together {[TERM_NAMES[i] for i in strong_indices]}")
    ablated = BEST_COEFFS.copy()
    ablated[strong_indices] = 0.0
    mean_loss, se = evaluate_with_high_precision(ablated, n_seeds, label="ZERO Top 3")
    delta = mean_loss - base_mean
    results.append({
        "term": "Top 3 strongest together",
        "zeroed_indices": strong_indices,
        "mean_loss": mean_loss,
        "stderr": se,
        "delta_from_baseline": delta
    })
    log(f"  → {mean_loss:.4f} ± {se:.4f}   (Δ = {delta:+.4f})\n")

    # Final results
    final_output = {
        "n_seeds": n_seeds,
        "growth_steps": growth_steps,
        "perturbation_strength": perturbation_strength,
        "baseline": {"mean": base_mean, "stderr": base_se},
        "ablations": results
    }

    with open(results_file, "w") as f:
        json.dump(final_output, f, indent=2)

    log(f"\n=== Study complete ===")
    log(f"Final results written to: {results_file}")

if __name__ == "__main__":
    # Very intensive but realistically completable run
    # 1000 seeds + 100 growth steps with 16 workers
    # (This is ~300x the original tiny study and already extremely statistically powerful)
    ablation_study(n_seeds=1000, growth_steps=100)