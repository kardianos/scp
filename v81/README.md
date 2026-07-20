# v81 — Stage-3 redesign (leave multi-fab product)

**Status:** OP executed; **P1 kernel port landed** (CPU); sandboxes P1/P2/P3 complete  
**Kernel auth:** accepted; `n_locks=0` leaves nuclear Cosserat path unchanged

## Why

Multi-fab B1: **F PASS**, **O FAIL** (long soft evaporates; hard shreds L).  
v80 thesis: free substrate + locks, not multiplet-on-grid.

## Result

**P1 wins** for Stage-3. In-kernel: free U(1) + typed locks; opposite attract; durable carriers.  
See `FINDINGS.md`, `path1_locks/KERNEL_PORT.md`.

## Docs

| File | Role |
|------|------|
| `OP.md` | Operational plan (executed) |
| `FINDINGS.md` | Top-level scorecard |
| `path1_locks/KERNEL_PORT.md` | Kernel land + K-L\* results |
| `path1_locks/KERNEL_DELTA.md` | Port contract (accepted) |
| `path1_locks/kernel_runs/` | K-L0/L1/L2 configs + outputs |
| `coordination/COMPARE.md` | Path rank |
| `council/SYNTHESIS.md` | Design council |

## Next

1. Medium-sourced bag (anti-lock deposits grid bag) or sequestration; \(v_t\) match for multi-rev.  
2. Optional GPU locks.  
3. Do **not** resume multi-fab product.

**Done:** J-unit Gauss fix; anti-lock lock-Higgs; N=256 campaign (`path1_locks/kernel_runs/n256/RESULTS.md`).
