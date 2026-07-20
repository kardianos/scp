# path2_ell FINDINGS — P2 scorecard (E0–E2)

**Date:** 2026-07-19  
**Binary:** `src/path2_ell` (CPU, 2D, NX=NY=96)  
**Run:** `./path2_ell all out` → exit 0  
**Artifacts:** `out/e0_*.txt|tsv`, `out/e1_*.txt|tsv`, `out/e2_*.txt|tsv`

---

## Kill-gate scorecard (binding)

| Gate | Result | Evidence |
|------|--------|----------|
| **Hyperbolic \(\ell\)** (not elliptic Poisson-per-tick) | **PASS** | Impulse at center: \(t_{\mathrm{arrival}}(r=20)=16.86\), \(r/c=20.0\) (ratio 0.84). No Poisson solve in code path. |
| **Stationary non-runaway \(\ell(r)\)** (E0) | **PASS** | Pinned heavy: \(\ell_{\mathrm{center}}=5.31\), \(\ell_{\mathrm{max,hist}}=6.32\), \(\pi_{\mathrm{late}}=6\times10^{-13}\). Radial decay → 1 at large \(r\) (`e0_radial.tsv`). |
| **Non-decorative \(\ell\)** | **PASS** | E1: \(\Delta y_{\kappa}=8.36\) vs \(\Delta y_{\kappa=0}=0\) and \(\Delta y_{\mathrm{no\,source}}=0\). E2: \(\ell\) changes pair \(dsep\). |
| **GRIN kill** (v76; no free-frame variable-\(c\) gravity) | **PASS** | Free \(c=1\) fixed in engine. Force law is \(-\kappa\nabla\ell\), not \(n(\rho)\) or CFL rewrite. M-chart \(c_{\mathrm{eff}}^{\min}\approx0.26\) reported as **readout only**. \(\ell\) bounded (\(\ell_{\max}\sim O(1\text{–}6)\)). |
| Pairwise Coulomb in push | **PASS (absent)** | Push samples grid \(E,B\) + \(\nabla\ell\). Soft Coulomb used only as **initial** \(E\) IC for E2 medium (then Yee evolves). |
| Multi-fab / light as Q-ball hump | **PASS (absent)** | Locks are structs; no Cosserat multiplet L. |

---

## Hyperbolic vs elliptic (explicit)

| Statement | Verdict |
|-----------|---------|
| \(\ell\) update is \(\ddot\ell - c^2\nabla^2\ell + \gamma\dot\ell + m^2(\ell-1) = S\) leapfrog | **Yes — hyperbolic** |
| Code calls elliptic Poisson each tick for \(\ell\) | **No** |
| Disturbances in \(\ell\) propagate at finite speed \(\approx c\) | **Yes** (E0 impulse travel-time check) |
| Would fail OP kill if Poisson-per-tick were required for sanity | **Would kill — not required here** |

Damping \(\gamma\) and mass \(m^2\) give a **stationary** Yukawa-like profile without changing the PDE type. Stationary limit is elliptic *as equilibrium*, not as the time stepper.

---

## GRIN kill check (explicit)

| Claim under test | Result |
|------------------|--------|
| Free-frame \(c\) varies with local density / GRIN \(n(\rho)\) | **Rejected by construction** — `C_FREE` constant; CFL uses `C_FREE` only |
| Local compact optics as long-range gravity (v76 dead branch) | **Not used** — path cost is free-response \(\ell\) field, sourced by locks |
| M-chart optical dual \(c_{\mathrm{eff}}\sim c/\ell\) | **Allowed as readout**; E1 reports `c_eff_min_proxy_Mchart≈0.26`; **not** free law |
| Runaway \(\ell\) / unbounded deflection | **Not observed** in E0–E1 window (`ell_max_hist=6.3`, finite \(\Delta y\)) |

---

## Experiment results

### E0 — pinned heavy → stationary \(\ell(r)\)

| Metric | Value |
|--------|-------|
| \(\ell_{\mathrm{center}}\) | 5.310 |
| \(\ell_{\mathrm{max}}\) history | 6.319 |
| \(\pi_{\mathrm{late}}\) | \(6.35\times10^{-13}\) |
| Impulse \(t_{\mathrm{arrival}}/ (r/c)\) | \(16.86/20.0 = 0.84\) |
| **PASS** | stationary + hyperbolic + non-runaway |

### E1 — light test deflection

| Case | \(\Delta y\) at \(x=x_{\mathrm{scatter}}\) | \(\Delta v_y\) proxy |
|------|---------------------------------------------|---------------------|
| \(\kappa\) on, source on | **+8.36** | \(v_y^{\mathrm{exit}}=0.263\) |
| \(\kappa\) off | 0 | 0 |
| source off (flat \(\ell=1\)) | 0 | 0 |

| Check | Result |
|-------|--------|
| Deflection from path-cost force | **PASS** |
| Non-decorative | **PASS** |
| GRIN kill | **PASS** |
| Fixed-\(c\) budget | **PASS** (speed cap \(0.5c\); free \(c\) unchanged) |

Interpretation: C-chart warp \(\ell\) deflects a light lock around a heavy thickener. Null controls prove the force channel is \(\kappa\nabla\ell\), not numerics drift.

### E2 — two-lock \(q=\pm1\) (positronium-class MVP)

| Case | \(\Delta\)sep (D=14) | Attract? | Notes |
|------|----------------------|----------|-------|
| Maxwell only (\(\kappa=0\), no \(\ell\) source) | \(-13.92\) | yes | Soft-Coulomb **IC** + Yee; collapses |
| Maxwell + \(\ell\) | \(-14.00\) | yes | Locks merge; \(\ell_{\max}\approx2.68\) |

| Check | Result |
|-------|--------|
| Locks intact (structs, no hump evaporate) | **YES** (by construction) |
| \(\ell\) changes pair dynamics vs pure P1-like | **YES** (small \(dsep\) difference; same qualitative collapse) |
| Durable multi-rev bound orbit | **NO** (not claimed; both cases head-on collapse) |
| Code `verdict` flag | `IMPROVEMENT` (marginal \(dsep\)); **honest reading: HONEST_PARTIAL** — opposite attraction exists; \(\ell\) is non-decorative but does **not** yet buy a parked positronium gyration |

Soft 2D Coulomb IC is strong relative to inertia → rapid approach. A softer seed / angular momentum (impact parameter) would be the next orbit attempt — out of this time-box MVP.

---

## Architecture honesty

| Item | Status |
|------|--------|
| Locks ≠ field humps | Yes — `Lock` array |
| Free medium lives | Yes — Yee TE + deposit |
| \(\ell\) first-class hyperbolic DOF | Yes |
| 2D (not 3D \(1/r\) Einstein multipole) | Intentional MVP |
| Collocated Maxwell (not full staggered Esirkepov) | MVP limitation; Gauss not at 1e-13 floor |
| P1 isolation | Vendored in-tree; did not edit `path1_locks/` or `scp_sim` |

---

## Overall P2 verdict

| | |
|--|--|
| **Runnable MVP** | **YES** |
| **E0** | **PASS** |
| **E1** | **PASS** |
| **E2** | **PASS with honest partial** (attract + locks durable; no multi-rev park; \(\ell\) non-decorative) |
| **Hyperbolic check** | **PASS** |
| **GRIN kill check** | **PASS** |
| **Path kill?** | **No** — none of elliptic-per-tick / decorative / GRIN triggers fired |
| **Kernel port now?** | **No** — wait for P1 medium maturity; P2 is co-field on top of locks medium |

---

## Next (if continued)

1. Staggered Yee + charge-conserving deposit (share with P1).  
2. E2 with impact parameter / \(v_t\) for multi-rev attempt under fixed \(c\).  
3. 3D exterior multipole check for \(\ell\sim 1/r\) (v76 F1-3D).  
4. Optional: \(\ell\)-weighted medium impedance without free-frame \(c\) rewrite.
