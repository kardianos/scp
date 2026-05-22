# Phase 6 — All Pressure Experiments Summary

**Date**: 2026-05-21  
**Location**: `phase6_algebraic_rigidity/`

This document consolidates the results from all intensive Differential Evolution searches performed on the new pressure term bank.

---

## Overview of the Four Experimental Modes

We ran four configurations of the pressure term bank on the relational graph model under the exact living candidate:

| Mode | Description | What the pressure terms affect |
|------|-------------|--------------------------------|
| **Additive** | Pressure terms only influence node attachment bias (living activity) | Attachment / growth only |
| **Modulation** | Pressure terms only modulate the quadratic strength (`λ` and/or `μ`) inside the living candidate update | Quadratic self-interaction only |
| **Both** | Pressure terms affect **both** attachment bias **and** quadratic modulation | Full dual-channel effect |
| **Full Combined** | Second independent run of the "Both" configuration (for robustness) | Same as "Both" |

Each mode was run at least once with an intensive search (population 20, 12 generations, multiple random seeds, up to 16 workers). The "Both" and "Additive" modes each received two independent runs.

---

## Master Results Table

| Mode                  | Log File                              | Best Loss   | Notes |
|-----------------------|---------------------------------------|-------------|-------|
| **Both** (full_combined) | `phase6_pressure_full_combined.log` | **-0.7132** | **Strongest result overall** |
| **Additive** (latest)    | `phase6_pressure_additive.log`      | -0.6918     | Best Additive run |
| **Modulation** (1st)     | Earlier modulation run              | -0.6796     | Best Modulation run |
| **Both** (the other)     | `phase6_pressure_both.log`          | -0.6770     | Second "Both" run |
| **Modulation** (latest)  | `phase6_pressure_modulation.log`    | -0.6672     | Second Modulation run |
| **Additive** (1st)       | Earlier additive run                | -0.6576     | First Additive run |

**Current ranking (best loss):**

1. **Both** (full_combined) — **-0.7132** ← best configuration found so far
2. Additive (latest) — -0.6918
3. Modulation (1st) — -0.6796
4. Both (the other run) — -0.6770
5. Modulation (latest) — -0.6672
6. Additive (1st) — -0.6576

---

## Detailed Coefficients by Run

### 1. Both — full_combined (Best overall: -0.7132)

| Term                              | Coefficient |
|-----------------------------------|-------------|
| term_rank_deviation               | +1.9100     |
| term_rank_abs                     | -1.4633     |
| term_closure_defect               | +0.7276     |
| term_balance                      | -2.3360     |
| term_rho_rank_interaction         | +0.4099     |
| term_rho_closure_interaction      | +1.0997     |
| term_biv_magnitude                | +1.9638     |
| term_rank_biv_interaction         | +0.2530     |

**Saved to:** `phase6_best_coefficients.json` (at the time of writing)

---

### 2. Additive — latest clean run (-0.6918)

| Term                              | Coefficient |
|-----------------------------------|-------------|
| term_rank_deviation               | -0.5070     |
| term_rank_abs                     | +0.1703     |
| term_closure_defect               | -0.7817     |
| term_balance                      | -0.6608     |
| term_rho_rank_interaction         | -2.1469     |
| term_rho_closure_interaction      | +0.0949     |
| term_biv_magnitude                | -1.7663     |
| term_rank_biv_interaction         | +1.4434     |

---

### 3. Modulation — first run (-0.6796)

| Term                              | Coefficient |
|-----------------------------------|-------------|
| term_rank_deviation               | -0.9973     |
| term_rank_abs                     | -1.3027     |
| term_closure_defect               | -1.3247     |
| term_balance                      | +0.2566     |
| term_rho_rank_interaction         | -2.2702     |
| term_rho_closure_interaction      | +2.1236     |
| term_biv_magnitude                | +1.2088     |
| term_rank_biv_interaction         | +0.8294     |

---

### 4. Both — second independent run (-0.6770)

| Term                              | Coefficient |
|-----------------------------------|-------------|
| term_rank_deviation               | -1.6431     |
| term_rank_abs                     | -2.2301     |
| term_closure_defect               | +2.1799     |
| term_balance                      | +0.2605     |
| term_rho_rank_interaction         | -1.4308     |
| term_rho_closure_interaction      | +1.0083     |
| term_biv_magnitude                | +1.0876     |
| term_rank_biv_interaction         | -1.4489     |

---

### 5. Modulation — latest clean run (-0.6672)

| Term                              | Coefficient |
|-----------------------------------|-------------|
| term_rank_deviation               | -0.9973     |
| term_rank_abs                     | -1.3027     |
| term_closure_defect               | -1.3247     |
| term_balance                      | +0.2566     |
| term_rho_rank_interaction         | -2.2702     |
| term_rho_closure_interaction      | +2.1236     |
| term_biv_magnitude                | +1.2088     |
| term_rank_biv_interaction         | +0.8294     |

(Note: This run produced almost identical coefficients to the first Modulation run, suggesting good convergence for this mode.)

---

### 6. Additive — first run (-0.6576)

| Term                              | Coefficient |
|-----------------------------------|-------------|
| term_rank_deviation               | +0.8222     |
| term_rank_abs                     | +0.5472     |
| term_closure_defect               | -0.4597     |
| term_balance                      | +2.4936     |
| term_rho_rank_interaction         | +2.0772     |
| term_rho_closure_interaction      | +1.8721     |
| term_biv_magnitude                | +0.7139     |
| term_rank_biv_interaction         | +0.0473     |

---

## Key Observations Across All Runs

- **"Both" mode** produced the single best loss (-0.7132). Allowing pressure to act on both attachment and quadratic modulation appears to be the most powerful configuration.
- **Additive mode** showed the largest improvement between its two independent runs (from -0.6576 to -0.6918), suggesting it is still sensitive to random seed and may benefit from even more generations or restarts.
- **Modulation mode** was the most consistent between its two runs, with very similar coefficient patterns (especially strong negative weight on `term_rho_rank_interaction`).
- The terms `term_rank_deviation`, `term_rank_biv_interaction`, `term_rho_rank_interaction`, and `term_closure_defect` appear repeatedly with large coefficients across modes — these are the most promising signals for algebraic rigidity.
- Sign flips between modes (e.g., `term_closure_defect` being positive in one "Both" run and negative in others) indicate that the optimal *direction* of some effects depends on whether pressure can also modulate the quadratic terms.

---

## Files Generated

- `phase6_best_coefficients.json` (overwritten by later runs; contains the most recent result)
- Mode-specific logs:
  - `phase6_pressure_additive.log`
  - `phase6_pressure_modulation.log`
  - `phase6_pressure_both.log`
  - `phase6_pressure_full_combined.log`

When the script was updated to include the mode in the filename, later runs began saving to `phase6_best_coefficients_<mode>.json`.

---

## Next Recommended Steps

1. **Post-hoc ablation study** on the top 2–3 configurations (especially the strongest "Both" at -0.7132 and the best Additive at -0.6918).
2. **Validation runs** on the full discrete simulator (500–600+ nodes) using the winning coefficient sets from the searches.
3. Export the best validated configurations for the Lean side (`Phase6Rigidity.lean`).
4. Begin formalizing the most effective pressure + rigidity combination in Lean.

---

*This document was assembled on 2026-05-21 from all completed Differential Evolution pressure experiments in the Phase 6 exploration.*