# Gen 4 Progress

## Phase A: Remote Setup
- [x] Upload source files
- [x] Fix includes and compile
- Status: DONE

## Phase B: Seed Generation + Simulations
- [x] CB15 seed generated
- [x] CB15 sim T=200 complete -- SURVIVED (42/42 alive, S_final=0.816, S_mean=2.987)
- [x] CB15 analysis complete
- [x] S20 seed generated
- [x] S20 sim T=200 complete -- SURVIVED (42/42 alive, S_final=2.341, S_mean=2.960)
- [x] S20 analysis complete
- [x] UDD_R4 seed generated
- [x] UDD_R4 sim T=200 complete -- SURVIVED (42/42 alive, S_final=2.300, S_mean=2.121)
- [x] UDD_R4 analysis complete

## Phase C: Results Summary
- [x] Summary printed (see below)

## Phase D: Download
- [x] CB15_output.sfa downloaded (10.9 GB)
- [x] S20_output.sfa downloaded (10.9 GB)
- [x] UDD_R4_output.sfa downloaded (10.9 GB)
- [x] All JSONs downloaded
- [x] All diags downloaded

Note: Remote large SFA outputs were deleted by concurrent agent activity.
All local copies appear complete (~10.9 GB each, consistent sizes).

## Results Summary

ALL 3 candidates survived T=200 (42/42 frames alive each).

| Candidate | S_final | S_mean | Alive | E_pot(final) | P_int retention | Best for Step 2? |
|-----------|---------|--------|-------|-------------|----------------|-----------------|
| **S20** | **2.341** | **2.960** | 42/42 | -430.7 | 40.2% | **BEST** |
| UDD_R4 | 2.300 | 2.121 | 42/42 | -391.9 | 75.0% | SECOND |
| CB15 | 0.816 | 2.987 | 42/42 | -133.9 | 16.8% | THIRD |

### Analysis:
- **CB15**: Highest S_mean (2.99) but worst S_final (0.82). Strong oscillations.
  The counter-braid structure pulses strongly -- high peaks but deep troughs.
  At T=200 it caught a trough (E_pot only -134). Low P_int retention (17%).
  The structure is alive but not stable -- it may be breathing/oscillating.

- **S20**: Best S_final (2.34), strong E_pot (-431 at T=200), decent retention (40%).
  More sustained binding throughout. Best candidate for Step 2 modification.

- **UDD_R4**: Highest P_int retention (75%!). S_final=2.30 very close to S20.
  Most structurally compact (aspect ratio stays near 1.1-1.5 vs 2-5 for others).
  Genuinely 3D structure, very robust. Strong second candidate.

### Recommendation for Step 2:
1. **S20** — best overall survival score at T=200
2. **UDD_R4** — best retention and structural stability, strong runner-up
3. CB15 — drop from Step 2 (oscillating, low final score)
