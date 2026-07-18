# v77 Phase 2 Convergence — M1 ∧ RC1

**Date:** 2026-07-19  
**Status:** **STOP Phase 2 (A-depth)** — agreement discovered  
**Not:** NULL / DISPROVE  
**Map:** [`CAMPAIGN_MAP.md`](CAMPAIGN_MAP.md)  
**Score:** [`work/TU/PHASE2_SCORE.md`](work/TU/PHASE2_SCORE.md)

---

## 1. Terminal result

| Goal | Result |
|------|--------|
| Comprehensive dynamic full Maxwell (M1) | **MET** — `m1_claim=true` |
| Close R-compose (RC1) | **MET** — `rc1_claim=true` |
| Cohesive multi-agent agreement | **MET** — mandatory CP path all ADOPTED |
| Critical null / disprove | **Not found** |

**Mandatory path completed:**

```text
CP0 → CP-M1-SPEC → CP-M1-NUM → CP-RC1-SPEC → CP-RC1-NUM → CP-U-FINAL
 ✓        ✓            ✓            ✓             ✓            ✓
```

---

## 2. What was demonstrated

### M1 — true 2D dynamic Maxwell (`work/NE/sandbox_m1_2d.py`)

| Gate | Result (parent `--quick`) |
|------|---------------------------|
| R0 M0 regression | PASS |
| G1 vacuum | PASS |
| G2 **true 2D** beam (not 1D TEM) | PASS \(v/c\sim 1.04\) |
| G3 off-unit \(c=0.5\) | PASS |
| G4 energy drift | PASS \(\sim 0.5\%\) |
| G5 dynamic Gauss + Cont sources | PASS |
| G6 div B | PASS \(\sim 10^{-16}\) |
| G7 Faraday | PASS |
| G8 incomplete Ampère adversary | PASS (full waves, adv \(v=0\)) |
| G9 BC honesty | PASS |
| **m1_claim** | **true** |

Also: `Maxwell2D.step(rho_Q, J…)` API for RC1.

### RC1 — co-field closes R-compose (`work/NM/sandbox_rc1_cofield.py`)

| Item | Result |
|------|--------|
| Same 2D grid | free-capacity \(\psi\) (F1 Laplace) + dynamical Yee \(\mathbf{E},\mathbf{B}\) |
| Maxwell source | **import NE `Maxwell2D`** — not Φ-only |
| Fixed multi-locks | 2 Gaussians |
| Configs vacuum / neutral / same / opposite | RG0–6 **all True** |
| Shared \(c\) | \(c=1/\sqrt{\varepsilon\mu}=1\) |
| TE-IA1 \(\psi\neq\) EM | neutral: \(\psi\neq0\), \(\|E\|=0\) |
| **rc1_claim** | **true** |

```bash
cd v77/work/NE && python3 sandbox_m1_2d.py --quick
cd v77/work/NM && python3 run_rc1.py
```

---

## 3. Residuals (honest, non-blocking for Phase 2)

1. Co-field \(\psi\) is **2D** free Laplace (log multipole class). Seed 3D \(1/r\) path-cost remains v76 / separate diagnostic — not co-evolved in RC1 3D.
2. Locks are **fixed** (RC2 moving locks not required for this STOP).
3. Force-sign hierarchy soft; full Lorentz \(v\times B\) dynamics deferred to RC2.
4. M2 3D Yee Maxwell optional stretch only.
5. Poisson isomorphism Occam residual unchanged.

---

## 4. Kill criteria

K-M1 / K-RC1 **did not fire**. No program DISPROVE.

---

## 5. Relation to Phase 1

| Phase | STOP meaning |
|-------|----------------|
| **Phase 1** `CONVERGENCE.md` | Composition-tier unity (separate Maxwell + dual lite) |
| **Phase 2** this file | **Depth:** true 2D dynamic Maxwell + **one-grid co-field** R-compose close |

Phase 1 remains valid. Phase 2 **deepens** examination as mapped.

---

## 6. Orchestrator call

**STOP Phase 2 (A-depth).** Campaign complete under CAMPAIGN_MAP mandatory path. Optional later: M2, RC2 — not required for this close.
