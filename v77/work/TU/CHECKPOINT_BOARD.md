# Phase 2 Checkpoint Board — FINAL

**Maintainer:** TU  
**Map:** `v77/CAMPAIGN_MAP.md`  
**Phase:** **P2-R5 FINAL** · **Updated:** 2026-07-19 (TU-035+)

| ID | Name | Status | Stamps in | Notes |
|----|------|--------|-----------|-------|
| **CP0** | Phase-2 kickoff | **ADOPTED** | O, TU, TE, NE | |
| **CP-M1-SPEC** | M1 gates frozen | **ADOPTED** | TE, NE | `m1_gates_v0.md` |
| **CP-M1-NUM** | M1 numeric green | **ADOPTED** | O-011; TE-017; NE; TU | `m1_claim=true` |
| **CP-RC1-SPEC** | RC1 joint API | **ADOPTED** | TE, NE, TM, TD(ψ), NM-011 | joint state + Maxwell2D API |
| **CP-RC1-NUM** | RC1 co-field green | **ADOPTED** | measured `rc1_claim=true`; NM export stamps ADOPT; TU-035 | RG0–6 True |
| CP-M2-SPEC / NUM | 3D Maxwell | PENDING | — | **optional stretch** (not required) |
| CP-RC2-SPEC / NUM | Moving locks | PENDING | — | **optional stretch** |
| CP-M3 | Charged matter | PENDING | — | stretch |
| **CP-U-FINAL** | Phase-2 unity | **ADOPTED** | **TU-036 ADOPT**; FOR_O | **STOP Phase 2 (A-depth)** |

**Status values:** PENDING | IN_PROGRESS | ADOPTED | DEFERRED | REJECTED | KILLED

---

## Mandatory path — COMPLETE

```text
CP0 → CP-M1-SPEC → CP-M1-NUM → CP-RC1-SPEC → CP-RC1-NUM → CP-U-FINAL
 ✓        ✓            ✓            ✓             ✓            ✓
```

**Phase-2 terminal bar:** **M1 ∧ RC1** green with co-agreement → **MET**.

---

## CP-M1-NUM — ADOPTED (recap)

| Item | Value |
|------|--------|
| Export | `work/NE/outputs/m1_result.json` |
| **m1_claim** | **true** |
| Gates | R0 + G1…G9 all PASS |
| True-2D | G2 `true_2d_packet`; Ampère adversary G8 PASS |
| API | `Maxwell2D.step` / `fields`; `rc1_ready=true` |

---

## CP-RC1-SPEC — ADOPTED

| Party | Vote | Evidence |
|-------|------|----------|
| TE | ADOPT | `rc1_em_interface_v0.md` |
| NE | ADOPT | `maxwell_api_rc1.md` + Maxwell2D freeze |
| TM | ADOPT | `rc1_joint_state_v0.md` STATE-RC1; RG0–RG6 |
| TD | ADOPT (ψ) | `psi_time_dep_rc1_v0.md` T0 default |
| NM | ADOPT | NM-011 |
| TU | ADOPT (board) | TU-035 |

Hard rule held: dynamical \((\mathbf{E},\mathbf{B})\), not Φ-only for claim.

---

## CP-RC1-NUM — ADOPTED

| Item | Value |
|------|--------|
| Export | `work/NM/outputs/rc1_result.json` |
| Summary | `rc1_summary.txt`; `rc1_tm_gates.json` |
| **rc1_claim** | **true** |
| RG0…RG6 | **all True** |
| Maxwell | `import_NE_sandbox_m1_2d.Maxwell2D`; `used_NE_Maxwell2D=true` |
| ψ | 2D free Laplace on **same grid** |
| Locks | fixed multi-lock (N=2); vacuum/neutral/same/opposite |
| Shared \(c\) | \(C_{\mathrm{LOCAL}}=1/\sqrt{\varepsilon\mu}=1\); JC1 |
| Not Φ-only | `rc1_not_phi_only`; `E_origin=free_maxwell_full` |
| m1 inherited | `m1_claim_inherited=true` |
| JSON stamps | CP-RC1-SPEC=ADOPT; CP-RC1-NUM=ADOPT |

NM-015 DEFER was pre-run; superseded by measured `rc1_claim=true` (+ parent run).

---

## CP-U-FINAL — ADOPTED

| Item | Content |
|------|---------|
| **TU STAMP** | **ADOPT** (TU-036) |
| Verdict | **STOP Phase 2 (A-depth)** — not DISPROVE |
| Package | `PHASE2_SCORE.md` **FINAL** |
| Residual | R-compose **closed at RC1**; RC2/M2 optional stretch only |
| FOR_O | Confirm STOP Phase 2; optional authorize M2/RC2 later |

---

## Stamp log (final summary)

| CP | Outcome |
|----|---------|
| CP0 | ADOPTED |
| CP-M1-SPEC / NUM | ADOPTED |
| CP-RC1-SPEC | ADOPTED (TE+NE+TM+TD+NM) |
| CP-RC1-NUM | ADOPTED (`rc1_claim=true`, RG0–6) |
| **CP-U-FINAL** | **ADOPTED** — Phase-2 STOP (A-depth) recommended |

---

## Reproduce

```bash
cd /home/d/code/scp/v77/work/NE && python3 sandbox_m1_2d.py --quick
cd /home/d/code/scp/v77/work/NM && python3 sandbox_rc1_cofield.py
# or: python3 run_rc1.py
```
