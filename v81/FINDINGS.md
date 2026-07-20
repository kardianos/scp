# v81 FINDINGS — Stage-3 parallel redesign + P1 kernel port

**Date:** 2026-07-19  
**OP:** `v81/OP.md`  
**Status:** three sandboxes complete; **P1 kernel port landed** (CPU)

---

## Executive summary

Left multi-fab product (O FAIL). Three independent sandboxes, then user-accepted P1 kernel delta:

| Path | Status | Headline |
|------|--------|----------|
| **P1** monist gauge + typed locks | **DONE / KERNEL LIVE** | Sandbox L0–L2 PASS; in-kernel K-L0 PASS, K-L1 attract, K-L2 soft PASS\* |
| **P2** hyperbolic \(\ell\) + locks | **DONE** | \(\ell\) hyperbolic, non-decorative, GRIN-clear; no multi-rev park |
| **P3** token / hop-cap CA | **DONE** | Exact ledger + durable vortices; **no** EM-like opposite attract |

**Winner for atom delivery:** **P1** (sandbox + kernel).  
**Kernel docs:** `path1_locks/KERNEL_PORT.md`.

---

## Kernel scorecard (P1)

| Gate | Result |
|------|--------|
| `n_locks=0` regression | PASS (vacuum zero) |
| K-L0 pinned + Gauss floor | PASS (\(\sim10^{-13}\)) |
| K-L1 opposite attract | PASS (sign); mono soft |
| K-L2 free motion + Gauss | **PASS** — Gauss floor holds for full \(T=200\) after J-unit fix |
| K-L3 soft+\(v_t\) | PARTIAL — no multi-rev park yet |
| Pairwise Coulomb / multi-fab | Absent (soft core is optional short-range regularization only) |

**N=256 + anti-lock:** large box clean; anti-lock reduces \(\Delta\mathrm{sep}\) vs baseline; multi-rev park still open (`path1_locks/kernel_runs/n256/RESULTS.md`). Cosserat Higgs **not** used.

**Next:** medium-sourced bag / \(v_t\) match / sequestration; optional GPU.

---

## Reproduce

```bash
# Sandboxes
make -C v81/path1_locks && ./v81/path1_locks/bin/locks_pic all v81/path1_locks/results
make -C v81/path2_ell/src && ./v81/path2_ell/src/path2_ell all v81/path2_ell/out
cd v81/path3_token && python3 src/run_all.py

# Kernel P1
make -C sfa sim/scp_sim && cp sfa/sim/scp_sim bin/
OMP_NUM_THREADS=8 bin/scp_sim v81/path1_locks/kernel_runs/kl0/kl0.cfg
```

Detail: `coordination/COMPARE.md`, `path1_locks/KERNEL_PORT.md`.
