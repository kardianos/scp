# v79 RESULTS — Multi-fabric Z=6 atom long-T

**Instance:** Vast 45259142, Tesla V100-SXM2-32GB, N=192 L=48, g=0.05, multi-fab B2 (C/Q/L), η=0  
**Seeds:** light nucleon profile ω=1.46 via `gen_pn_core` (Z octa R=8, N R=5.5, L R=22)  
**Date:** 2026-07-18/19

---

## A. Primary: `z6_a_n6_T800` (atom = C core + Q + L shells)

| Item | Value |
|------|--------|
| Status | **COMPLETE** |
| Wall | 6160.5 s (102.7 min), 96.76 ms/step |
| Steps | 63667 (T=800, dt_factor=0.025) |
| SFA | 10 frames, ~4.2 GB → `/space/scp/v79/z6_a_n6_T800.sfa` |
| Diag | `results/z6_a_n6_T800_diag.tsv` (401 rows) |
| Gauss CG | 689 iters → gauss_max(0)=9.91e-14 |

### Bookkeeping

| Metric | INIT (t=0) | MID (t≈400) | END (t≈799) | Δ END/INIT |
|--------|------------|-------------|-------------|------------|
| Q_phi | 3540.4 | 3040.7 | 3002.1 | **−15.2%** |
| Q_flux | 455.5 | 209.3 | 145.4 | **−68%** |
| Q_core | 2500 | 2610 | 2588 | +3.5% |
| r_core | 7.29 | 7.36 | 6.60 | shrink |
| E_total | 6141 | 6725 | 7104 | **+15.7% drift** |
| E_em | 6.88 | 3.78 | 0.36 | **collapse** |
| gauss_max | 9.91e-14 | 9.91e-14 | 9.91e-14 | floor OK |

---

## B. Control: `z6_n6_nuc_T800` (Z=6 N=6 core, L=0)

| Item | Value |
|------|--------|
| Status | **COMPLETE** |
| Wall | 6128.1 s (102.1 min), 96.25 ms/step |
| SFA | 10 frames, ~2.1 GB → `/space/scp/v79/z6_n6_nuc_T800.sfa` |
| Diag | `results/z6_n6_nuc_T800_diag.tsv` (401 rows) |
| Gauss CG | 687 iters → gauss_max(0)=9.32e-14 |

### Bookkeeping

| Metric | INIT (t=0) | MID (t≈400) | END (t≈799) | Δ END/INIT |
|--------|------------|-------------|-------------|------------|
| Q_phi | 3540.4 | 3040.7 | 3002.1 | **−15.2%** |
| Q_flux | 455.5 | 533.4 | 533.4 | **+17%** |
| Q_core | 2500 | 2610 | 2588 | +3.5% |
| r_core | 7.29 | 7.36 | 6.60 | shrink |
| E_total | 5100 | 5350 | 7973 | **+56.5% drift** |
| E_em | 6.36 | 6.02 | 6.13 | **stable ~6** |
| gauss_max | 9.32e-14 | 9.33e-14 | 9.33e-14 | floor OK |

INIT E lower than atom (no L matter energy). END runlog: ω_core≈1.345 (atom printed 0 — diag/core clock differ).

---

## C. Atom vs core — side-by-side (same C/Q seed physics)

| Metric | Atom END | Nuc END | Notes |
|--------|----------|---------|--------|
| Q_phi | **3002.1** | **3002.1** | **Identical** (to printed precision) |
| Q_core | **2588** | **2588** | **Identical** |
| r_core | **6.60** | **6.60** | **Identical** |
| Q_flux | 145.4 | 533.4 | Atom −388 vs nuc |
| E_em | 0.36 | 6.13 | Atom EM collapse; nuc holds Coulomb-like energy |
| E_total drift | +15.7% | +56.5% | Nuc worse energy bookkeeping |
| gauss_max | 1e-13 floor | 1e-13 floor | Both clean |

Early Q_phi trajectory (t=50…400) is **digit-identical** between runs. Global charge shed is therefore a **C/Q nuclear / radiation channel**, not L-driven.

### Interpretation

1. **Gauge integrity OK** on both — multi-fab Gauss holds T=800.
2. **No cold multi-clock atom.** Adding L shells does **not** stabilize Q_phi or freeze a bound atom; it **cancels net flux** and **kills E_em** while C/Q inventory tracks the bare-nucleus control.
3. **L sector acts as a neutralizing / free-radiation bath** relative to Q-charged matter: flux and E_em diverge strongly; matter Q does not.
4. **Core inventory** (Q_core ~2585, r~6.6) is the robust object; free Q_phi still sheds ~15% into the absorbing BC for both.
5. **Energy drift** is large on both (worse without L). Not a stationary bound-state ledger; absorbing BC + unresolved radiation / reconfiguration.

### Verdict

| Claim | Result |
|-------|--------|
| Multi-fab T=800 is numerically healthy (Gauss) | **PASS** |
| Hand-placed Z=6 + L=6 forms a cold parked atom | **FAIL** (flux collapse, E_em death, Q shed) |
| L sector binds opposite charge into atom-like multipole | **FAIL / null** — L destroys E_em vs L=0 control |
| C/Q core survives as Q-ball-like droplet | **PARTIAL** — Q_core holds; global Q and E not stationary |

---

## D. Products

| Path | Content |
|------|---------|
| `v79/results/z6_a_n6_T800_{diag.tsv,runlog,cfg}` | Atom diagnostics |
| `v79/results/z6_n6_nuc_T800_{diag.tsv,runlog,cfg}` | Core control diagnostics |
| `/space/scp/v79/z6_a_n6_T800.sfa` | Atom frames (4.2 GB) |
| `/space/scp/v79/z6_n6_nuc_T800.sfa` | Core frames (~2.1 GB) |

Optional next: volview charge/flavor snapshots; orbit/capture designs that force n=2 binding; kernel multi-fab only if Step-2 orbit graph still blocked.

## E. Instance

Vast **45259142** still up after both completes. Teardown after SFA downloads verified (do not destroy while rsync open).
