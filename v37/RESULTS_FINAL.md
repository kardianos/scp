# V37 Final Results and Session Summary

## What Was Accomplished

### Tools Built
| Tool | File | Purpose |
|------|------|---------|
| SFA Fragmentation Analyzer | `src/sfa_frag.c` | Connected-component cluster detection on SFA files |
| SFA Structure Analyzer | `src/sfa_structure.c` | Radial profiles, phase coherence, energy decomposition |
| Corrected Crossed Braid v2 | `src/v37_crossed_v2.c` | 7-bug fix (chirality, velocity, BC, SFA, fragmentation) |
| Knot Simulator | `src/v37_knot.c` | Trefoil + figure-eight knot geometries |
| Seed Runner | `src/v37_seedrun.c` | Load numpy seed, run C simulation with absorbing BC |
| Evolutionary Search | `backprop/field_evolve.py` | JAX evolutionary field optimization (GPU) |
| Backprop Search | `backprop/field_search.py` | JAX gradient-based field optimization (GPU) |
| Candidate Analyzer | `backprop/analyze_candidates.py` | Energy/spatial/phase analysis of evolved candidates |
| Seed Converter | `backprop/seed_to_c_init.py` | Numpy candidate → C binary for simulation |
| SFA Viewer Update | `sfa/viewer/volview.c` | Added keyboard help printout + depth-tested wireframe |

### Tools NOT Built (gaps)
- **CUDA seed runner**: v37_seedrun.c runs on CPU only. MUST be ported to CUDA for remote use.
- **CUDA evolutionary search**: field_evolve.py uses JAX (GPU), but the C validation runs on CPU.
- **Automated SFA download pipeline**: manual scp, no resume/checksum support built in.

### CRITICAL LESSON: Remote GPU Usage
The N=128 C simulation was run on a V100 but used CPU (OpenMP), NOT the GPU.
This wasted 61 minutes of GPU rental for CPU-equivalent performance. **All future
remote simulations MUST use CUDA or JAX to actually utilize the GPU.**

## Key Findings

### 1. All Previous "Surviving" Structures Were Fake
The fragmentation analyzer (sfa_frag) proved that ALL original crossed braid
experiments fragmented into 10-76 disconnected clusters by t=60. The old analysis
tools (total E_pot, P_int) failed to detect this because they tracked global
metrics, not structural coherence.

### 2. Braid3 Also Fragments
The braid3(z) — previously thought stable — fragments into 2-5 clusters and
drifts along z. It only appears to survive because periodic BC wraps the pieces
back together. With absorbing BC, it dissolves.

### 3. Evolutionary Search Finds Compact Binding From Random Seeds
The Game-of-Life-style evolutionary approach (selection + mutation + crossover on
raw 3D field values) successfully discovered field configurations with genuine
V(P) binding from pure random initialization:

| Search | T_eval | Best E_pot | R_rms | Pop | Gens |
|--------|--------|-----------|-------|-----|------|
| N=32 T=10 pop=128 | 10 | **-30.9** | 4.77 | 128 | 268 |
| N=32 T=10 pop=64 | 10 | -5.72 | 4.32 | 64 | 200 |
| N=32 T=20 pop=64 | 20 | -4.17 | 7.26 | 64 | 500 |
| N=32 T=30 pop=64 | 30 | -0.49 | 8.36 | 64 | 500 |

### 4. Best Candidate Survives T=200 at N=128
The best evolutionary candidate (T=10, pop=128, gen 268) was upscaled to N=128
and run in full C simulation with absorbing BC:
- t=0-30: Coalesces from 10 fragments into 1 cluster
- t=30-120: Coherent oscillon, E_pot oscillating -0.1 to -24 (breathing mode)
- t=120-200: Binding slowly weakens as absorbing BC drains energy
- t=200+: Binding fades to zero, structure dissolves

This is the longest any compact structure survived with absorbing BC.

### 5. Theoretical Framework: Oscillons
The Researcher agent established that this system supports oscillons:
- Derrick's theorem evaded via time dependence
- Breathing frequency ω≈0.3 far below mass gap m=1.5 → exponentially suppressed radiation
- The 90% energy loss is initial transient, not steady-state radiation
- Phase offsets δ={0, 3.0, 4.4} may be resonant tuning (cf. Gleiser-Krackow 2019)
- A well-formed oscillon could survive essentially indefinitely

### 6. SFA f32 Reader Bug Fixed
sfa_frag.c hardcoded `sizeof(double)` for column offsets, breaking f32 SFA files.
Fixed to use `sfa_dtype_size[dtype]`.

### 7. SFA fixup_index Bug Fixed
sfa_fixup_index() used first 4 bytes of compressed data as "checksum" instead of
proper CRC32 on decompressed data. Fixed to decompress, reverse BSS, and compute
real checksum.

## Candidate Files
All evolutionary search candidates stored in `backprop/candidates/`:
- `best_T10_pop128_gen268.npz` — strongest binding (-24.8)
- `best_T20_pop64_gen355.npz` — best phase coherence (0.938)
- Plus 7 intermediate checkpoints showing evolution trajectory

## SFA Files
- `data/seed_oscillon.sfa` — N=128 T=200 validation (f32, 68 frames, 2.5 GB)
  Viewable: `/home/d/code/scp/sfa/viewer/volview data/seed_oscillon.sfa`

## Documents Created
- `RESULTS_crossed_v2.md` — Bug fixes and fragmentation analysis
- `RESULTS_knots.md` — Trefoil + figure-eight results (all failed)
- `RESULTS_backprop.md` — Backprop approach results and insights
- `RESULTS_evolutionary_search.md` — Evolutionary search results
- `ANALYSIS_structure_contrast.md` — What works vs what fails (5 critical features)
- `SKEPTICISM_REPORT.md` — Skeptic agent concerns
- `CRITICAL_TESTS.md` — Proposed validation tests
- `backprop/PLAN.md` — Full backprop operational plan
- `backprop/DESIGN_search_v2.md` — Designer agent's CMA-ES proposal
