# work/NE — Numeric EM

**Agent:** NE (Numeric — EM)

## Phase 2 — M1 green + RC1 API freeze (current)

| Item | Status |
|------|--------|
| **CP-M1-SPEC** | **ADOPT** |
| **CP-M1-NUM** | **ADOPT** — parent `m1_claim=true` |
| **CP-RC1-SPEC** | **ADOPT** — API freeze |
| M1 result | [`outputs/m1_result.json`](outputs/m1_result.json) |
| RC1 API | [`maxwell_api_rc1.md`](maxwell_api_rc1.md) |
| Sandbox | [`sandbox_m1_2d.py`](sandbox_m1_2d.py) |

```bash
cd /home/d/code/scp/v77/work/NE
python3 sandbox_m1_2d.py --quick   # regression
# NM: from sandbox_m1_2d import Maxwell2D
```

`m1_claim=true` (R0∧G1–G9). NM RC1 must call `Maxwell2D.step` (not Φ-only).

## Round 3 — M0 claim confirmed

Parent-fixed `sandbox_full_maxwell_r2.py`: **`full_maxwell_claim=true`**.  
Export: `outputs/r2_result.json`. Re-run: `bash run_r3_confirm.sh`.

## Round 2 — Full Maxwell

| Item | Status |
|------|--------|
| Scheme | **Yee staggered 2D TE+TM** (`sandbox_full_maxwell_r2.py`) |
| **full_maxwell_claim** | **true** (FM1–FM7 / KG-F1–F5 PASS) |
| O-003/004 | Full Maxwell LIVE; compose with NM dual for (A) |

### Gates

| Gate | Result | Number |
|------|--------|--------|
| FM1 vacuum | **PASS** | max\|E\|=\|H\|=0 |
| FM2 wave unit | **PASS** | v/c_th=**1.0** (c_th=1) |
| FM2 wave off-unit | **PASS** | v/c_th=**1.0** (c_th=**0.5**, ε=4) |
| FM3 div B | **PASS** | divB_max=**0** |
| FM4 Faraday | **PASS** | discrete residual=**0** |
| FM5 Gauss static TE | **PASS** | continuum projection |
| FM6 continuity | **PASS** | dQ_rel=**0** |
| FM7 Coulomb 3D (R1) | **PASS** | Q=**11.48153604** |

### Run

```bash
cd /home/d/code/scp/v77/work/NE
python3 sandbox_full_maxwell_r2.py          # full dynamical Yee
python3 sandbox_full_maxwell_r2.py --quick  # smaller grids
python3 offline_r2_full_maxwell.py
```

### Files

| Path | Role |
|------|------|
| [`design_r2_full_maxwell.md`](design_r2_full_maxwell.md) | Scheme + gate design |
| [`sandbox_full_maxwell_r2.py`](sandbox_full_maxwell_r2.py) | Dynamical Yee TE+TM |
| [`outputs/r2_result.json`](outputs/r2_result.json) | Gate export |
| [`outputs/r2_summary.txt`](outputs/r2_summary.txt) | Summary |

### Tags (TE R2)

```text
sector=1
E_origin=free_maxwell_full
em_solver=free_maxwell_yee
gravity_solver=none
embedding_dim_dynamics=2
embedding_dim_coulomb=3
te_equation_match_r1=true
```

**full_maxwell_claim bar:** KG-F1 ∧ F2 ∧ F3 ∧ F4 ∧ F5 (TE `for_ne_kill_gates_r2.md` §4)

**Scope:** 2D full Maxwell (TE+TM) + 3D quasistatic Coulomb recovery.  
**Not claimed:** 3D Yee radiation multipoles; dual-channel ψ joint (see NM R2); real α_EM.

---

## Round 1 — Maxwell-lite (retained)

KG1–4 PASS; `sandbox_ne_r1_em.py`; `outputs/r1_*`.  
See R1 log NE-001…004. Lite demos still LIVE; R2 upgrades wave to dynamical E+H.
