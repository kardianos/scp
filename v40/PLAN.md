# V40 Execution Plan: Evolutionary Composite Particle Search

## Architecture

Two agents run in parallel:

1. **Search Agent**: Runs the evolutionary loop — generates candidates, launches
   simulations (GPU or CPU), analyzes SFA results, selects survivors, generates
   next generation.

2. **Monitor Agent**: Watches progress, provides high-level feedback, identifies
   when the search is stuck, suggests strategy changes, and alerts on resource issues.

Communication: both agents read/write to `v40/state/` directory.

## Directory Structure

```
v40/
  POSTULATES.md          — Theory and requirements
  PLAN.md                — This file
  state/
    generation.json      — Current generation number, candidates, scores
    best.json            — Best candidate so far (params + score + SFA path)
    monitor_notes.md     — Monitor agent's observations and suggestions
    search_log.md        — Search agent's decision log
  gen_NNN/
    candidate_MMM/
      config.cfg         — Simulation config
      output.sfa         — Simulation output
      diag.tsv           — Diagnostics
      analysis.json      — Stability metrics
  tools/
    analyze_sfa.c        — SFA stability analyzer (C, fast)
    modify_sfa.c         — SFA field modifier (add structures, perturb)
    score.sh             — Score a candidate from its analysis.json
```

## Evolutionary Loop (Search Agent)

### Generation 0: Braid Baseline

1. Run standard V34 braid at N=128, L=15, T=200, f32
2. Analyze SFA: extract S(t), energy profiles, phase structure
3. This establishes the baseline S_braid ≈ 0.6

### Generation 1: Single-Axis Perturbations

Generate 8-12 candidates by modifying the braid:
- Add a second braid at 90° rotation (test binding between non-parallel braids)
- Add a counter-rotating braid (opposite k_z)
- Modulate amplitude envelope (Gaussian → toroidal)
- Vary phase offsets around the CMA-ES optimum
- Add a localized oscillon near the braid core
- Increase/decrease ellipticity
- Test with/without background field

For each: run T=50, measure S(50), rank.

### Generation N+1: Refinement

Take top-4 from generation N. For each:
- Generate 3 variants (±10% in key parameters)
- Run T=100
- Measure S(100)
- Select top-4 overall

### Convergence Criteria

Stop when:
- A candidate survives T=500 with S > 0.7 (success)
- 10 consecutive generations show no S improvement > 0.01 (stuck)
- Total GPU-hours exceed budget (currently: 20 V100-hours)

### SFA Modification Strategy

Instead of always starting fresh, MODIFY existing SFA fields:

1. Load a stable SFA frame (e.g., from a surviving candidate at T=100)
2. Add a perturbation to the field data:
   - Overlay a rotated copy of the braid pattern
   - Add a localized energy pulse
   - Modify the phase structure in a sub-region
3. Write the modified frame as a new seed SFA
4. Resume simulation from this modified state
5. This tests whether the modification is STABLE — does the system absorb it
   or reject it?

This is the "evolutionary" approach: build up complexity incrementally,
testing each addition against the stability criterion.

## Monitor Agent Protocol

The monitor agent runs every 10 minutes and:

1. Reads `state/generation.json` — checks progress
2. Reads the latest `gen_NNN/candidate_MMM/analysis.json` files
3. Identifies trends:
   - Is S(t) improving across generations?
   - Are certain parameter regions consistently better?
   - Are there failure modes that repeat?
4. Writes observations to `state/monitor_notes.md`
5. If the search is stuck (no improvement for 3 generations):
   - Suggests a strategy change (different modification type, larger perturbation)
   - May suggest switching from refinement to exploration

## GPU Resource Strategy

- **Exploratory runs** (T=50, N=128): Use CPU (OMP_NUM_THREADS=16, ~10 min each)
- **Refinement runs** (T=100-200, N=128): Use CPU or single V100 ($0.30/hr)
- **Composite assembly** (T=200, N=256-512): Require V100 (~$0.30/hr, 1-4 hours each)
- **Budget**: Start with CPU for Gen 0-3, spin up V100 for promising candidates

The search agent should batch GPU jobs: spin up one V100, run 4-8 candidates
sequentially, then destroy the instance.

## Scoring Function

```
S(t) = 0.4 × energy_score + 0.3 × coherence_score + 0.2 × theta_score + 0.1 × shape_score

energy_score = max(0, -E_pot / E_pot_braid)       # normalized to braid baseline
coherence_score = P_int(t) / P_int(0)              # triple product retention
theta_score = min(1, theta_rms / theta_rms_eq)     # theta equilibration
shape_score = 1 / (1 + aspect_spread)              # penalize extreme anisotropy
```

Where E_pot_braid ≈ -80 (from V34/V39 braid runs).

A score of S > 1.0 means the candidate is MORE bound than a single braid.
This is the target for composite particles.

## Phase 1 Execution (Immediate)

1. Build `tools/analyze_sfa.c` — reads an SFA file, computes S(t) for each frame
2. Build `tools/modify_sfa.c` — reads SFA, adds perturbation, writes new SFA
3. Run braid baseline (Gen 0) on CPU
4. Analyze baseline, establish S_braid
5. Generate Gen 1 candidates (8 configs)
6. Run Gen 1 on CPU (parallel where possible)
7. Score and select top-4
8. Begin refinement loop
