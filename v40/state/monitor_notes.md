# Monitor Notes — 2026-03-23 20:50

## Current Status

**Generation**: 0 (baseline only)
**Gen 1 launched**: NO — `gen_001/` does not exist yet
**Running simulations**: NONE (`ps aux | grep scp_sim` returns nothing)
**Search agent**: Likely the Claude process on pts/1 (PID 17732, running since 09:39)

## Timeline

All v40 files were created in the last 2-3 minutes (20:46-20:49):
- Tools built: `analyze_sfa` and `modify_sfa` compiled at 20:47
- Baseline analyzed: `gen_000/baseline.json` written at 20:48
- `generation.json` updated at 20:48
- `seed.sfa` (4.8 MB) written at 20:49

The search agent appears to have JUST finished Gen 0 setup. Gen 1 candidate generation should begin imminently.

## Gen 0 Baseline Results

| Metric | Value |
|--------|-------|
| S_final | 0.807 |
| S_mean | 0.618 |
| E_pot_final | -86.6 |
| P_int_retention | 0.862 |
| alive_frames | 41/42 (98%) |
| phi_max range | 0.19 - 1.39 |
| theta_rms | ~0.01 (very small) |

The baseline braid shows strong oscillatory behavior — S(t) ranges from 0.08 to 1.58 across the 200-unit run. The object is alive almost the entire time but has large amplitude oscillations in binding energy and coherence.

## Observations

1. **No concerns yet** — the search agent is still in initial setup. This is expected given it had to build tools, run the baseline analysis, and create seed.sfa.

2. **Baseline is healthy**: S_final=0.807 matches the reference value (0.81). The baseline SFA source at `/home/d/code/scp/v39/data/braid_control.sfa` (1.6 GB) exists and was used correctly.

3. **Tools are ready**: Both `analyze_sfa` and `modify_sfa` are compiled and in `tools/`.

4. **Seed SFA created**: A 4.8 MB `seed.sfa` is ready — likely a single frame extracted for modification.

## What to Expect Next

The search agent should now:
1. Create `gen_001/` directory
2. Generate 8 candidate configs using `modify_sfa` (different perturbation types)
3. Launch simulations (T=50, N=128, CPU with OMP)
4. Each run should take ~5-10 min on 16 cores

Expected completion of Gen 1: ~21:30-21:50 (40-80 min from now)

## Suggestions

- None yet — too early to judge. Will re-check after Gen 1 candidates start completing.
- If no `gen_001/` appears within 10 minutes, the search agent may be stalled.
