# P1 kernel port — implemented + Gauss fix

**Date:** 2026-07-19  
**Status:** CPU kernel live; **moving-lock Gauss at ~1e-13 floor** after J-unit fix  
**GPU:** deferred

---

## Code

| File | Role |
|------|------|
| `sfa/sim/scp_locks.h` | CIC ρ, VB zigzag J (correct units), Boris, soft core, tracks |
| `sfa/sim/scp_sim.c` | Grid lock fields; Gauss includes `lock_rho`; push → Ampère |
| `sfa/sim/scp_config.h` | `n_locks`, `lock0`…, `locks_file`, `locks_track`, `lock_soft_r/k` |
| `sfa/Makefile` | depends on `scp_locks.h` |

**Default `n_locks=0`:** nuclear Cosserat path unchanged.

---

## Physics conventions

| Item | Convention |
|------|------------|
| Gauss | \(\nabla\cdot E = g\,(\rho_{\mathrm{matter}}+\rho_{\mathrm{lock}})\) |
| Lock ρ | CIC of \(q\); \(\int\rho\,dV=\sum q\) |
| Ampère | After push: \(E \mathrel{-}= dt\, g\, j_{\mathrm{phys}}\) (\(J_{\mathrm{lat}}=-j\)) |
| J units | \(J = q\cdot\Delta\mathrm{cell}/(dt\,dx^2)\) on faces (fixed: was off by \(dx\)) |
| Lorentz | \(F = -(g q)\,E\) so opposite locks attract |
| Soft core | Optional form-factor: \(k(1/r-1/r_0)\) for \(r<r_0\) (`lock_soft_r`, `lock_soft_k`) |

### Critical fix (session continuation)

Zigzag face current had an **extra factor of \(dx\)** in the numerator → Ampère over-kicked \(E\) whenever locks moved → `gauss_max` grew to \(O(1)\). Correct \(J = q\,\mathrm{dseg}/(dt\,dx^2)\) restores **machine-floor Gauss** under free motion.

---

## Config example

```text
complex_phi=1
complex_gauge=1
g_gauge=1.0
A=0
init=oscillon
n_locks=2
lock0=1,12,4,0,0,0,0.08,0,0
lock1=-1,12,-4,0,0,0,-0.08,0,0
lock_soft_r=2.0
lock_soft_k=0.15
locks_track=locks_track.tsv
```

---

## Validation (`kernel_runs/`)

| Exp | Result | Key numbers (N=48, g=1, vacuum Φ) |
|-----|--------|-------------------------------------|
| **K-L0** pinned single | **PASS** | `gauss~9e-14`; \(E_{\mathrm{em}}\) holds; \(Q_{\mathrm{flux}}\approx0.99\) |
| **K-L1** ± pinned | **PASS (sign)** | Opposite attract; mono soft (self-force/images) |
| **K-L2** free \(T=200\) | **PASS** | `gauss_max=1.2e-13` **entire run**; alive=2; \(E_\star\) exact; \(E_{\mathrm{em}}\) not null; sep soft approach |
| **K-L3** soft+\(v_t\) | **PARTIAL** | Gauss floor holds; ~0.7 joining-line rev (not multi-rev park); soft core prevents singularity |
| **reg0** `n_locks=0` | **PASS** | Vacuum zero |

### K-L2 (post J-fix)

```
gauss max over run = 1.201e-13
E* = 24 exact, alive = 2
sep 8.00 → 8.31 (min 7.16)
E_em holds (grows with dynamics, never nulls)
```

---

## Anti-lock (lock-Higgs) — not Cosserat Higgs

| Key | Meaning |
|-----|---------|
| `type=1` in lock CSV | Anti-lock: \(q=0\), no EM ρ/J |
| `lock_bag_r`, `lock_bag_k` | Attract charge locks toward anti-locks for \(r<r_{\mathrm{bag}}\) |
| `locks_medium_only=1` | Skip multiplet force (vacuum Φ + Maxwell + locks) — needed for N=256 CPU |

N=256 results: **`kernel_runs/n256/RESULTS.md`**. Anti-lock slows expansion vs baseline; multi-rev park still open.

## Honest limits / next

1. **Multi-rev park** — N=256 gives room; anti-lock helps mildly; need medium-sourced bag or sequestration / \(v_t\) match.  
2. **Self-force** — still contaminates pinned \(a_{\mathrm{rel}}(D)\).  
3. **`locks_medium_only` Gauss** — \(O(10^{-3})\) free residual (not Cosserat floor).  
4. **GPU** — not ported.

---

## Reproduce

```bash
make -C sfa sim/scp_sim && cp sfa/sim/scp_sim bin/
OMP_NUM_THREADS=8 bin/scp_sim v81/path1_locks/kernel_runs/kl0/kl0.cfg
OMP_NUM_THREADS=8 bin/scp_sim v81/path1_locks/kernel_runs/kl2/kl2.cfg
OMP_NUM_THREADS=8 bin/scp_sim v81/path1_locks/kernel_runs/kl3/kl3.cfg
```
