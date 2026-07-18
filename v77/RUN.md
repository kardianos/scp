# v77 Multi-Agent Run Protocol

**Date**: 2026-07-18  
**Status**: active  
**Problem**: [`PROBLEM.md`](PROBLEM.md)  
**Foci**: [`FOCI.md`](FOCI.md)

---

## 0. Ops shape (same spirit as v76, different agent graph)

| Keep from v76 | Change for v77 |
|---------------|----------------|
| Append-only per-agent logs | Seven agents: TE,NE,TD,ND,TM,NM,TU + O |
| Private work folders | Focus pairs (theory/numeric) not A/B/C/D |
| Read all logs; write only own | Seeded by v76 success, not blank slate |
| `FOR_X` cross-ideas | Unification agent TU + demo registry |
| Finite rounds; parent re-spawns | Goals G1–G3 multi-focus |

**Rejected:** worktrees; agents writing others’ logs; infinite daemons; kernel edits without human OK.

---

## 1. Directory map

```text
v77/
  PROBLEM.md
  FOCI.md
  RUN.md
  README.md
  logs/
    O_orchestrator.log
    TE_theory_em.log
    NE_numeric_em.log
    TD_theory_dynamics.log
    ND_numeric_dynamics.log
    TM_theory_matter.log
    NM_numeric_matter.log
    TU_theory_unification.log
  work/
    TE/ NE/ TD/ ND/ TM/ NM/ TU/
```

---

## 2. Ownership

| ID | Log | Work | Type |
|----|-----|------|------|
| **O** | `logs/O_orchestrator.log` | (parent notes optional in work/ORCHESTRATION.md) | orchestrator |
| **TE** | `logs/TE_theory_em.log` | `work/TE/` | theory EM |
| **NE** | `logs/NE_numeric_em.log` | `work/NE/` | numeric EM |
| **TD** | `logs/TD_theory_dynamics.log` | `work/TD/` | theory dynamics |
| **ND** | `logs/ND_numeric_dynamics.log` | `work/ND/` | numeric dynamics |
| **TM** | `logs/TM_theory_matter.log` | `work/TM/` | theory matter |
| **NM** | `logs/NM_numeric_matter.log` | `work/NM/` | numeric matter |
| **TU** | `logs/TU_theory_unification.log` | `work/TU/` | theory unification |

---

## 3. Append-only rules

1. Owner only appends to their log.
2. Never rewrite history.
3. Prefer:

```bash
cat >> /home/d/code/scp/v77/logs/<LOG> << 'EOF'
===== ENTRY <stamp> | <ID>-00N | <type> =====
Agent: <ID>
Round: <n>
Tags: ...
---
...
FOR_TE: ...
FOR_NE: ...
===== END =====
EOF
```

4. **Types:** `bootstrap` | `proposal` | `finding` | `attempt` | `numeric` | `theory` | `checkin` | `blocker` | `handoff` | `demo` | `kill` | `meta`

5. **short-id:** `TE-001`, `NE-012`, `O-003`, …

6. **Check-in required** each round: read other logs (at least partner + TU + O + one other focus); respond to `FOR_<self>`.

---

## 4. Pair protocol

Within a focus (EM, Dynamics, Matter):

- Theory agent writes laws first-ish; numeric implements kill-gates.
- Numeric may proceed from v76 seed + provisional equations if theory late; tag `provisional`.
- Both tag demos with IDs for TU registry (`D-EM-coulomb`, …).

Unification (TU) reads everyone; updates `WORLD.md` / `DEMOS.md` each round.

---

## 5. Constraints

- No `scp_sim` / `sfa.h` changes without explicit human authorization.
- New sandboxes only under `v77/work/{NE,ND,NM}/` (and theory notes under T*).
- May **read and cite** v76 code; prefer copy/adapt into v77 work dirs, not edit v76 logs.
- Vast/GPU: only if O authorizes for large 3D; default local Python.

---

## 6. Round lifecycle

```text
O: summary entry → spawn TE,NE,TD,ND,TM,NM,TU
agents: read → work → append → checkin → handoff
O: read all logs → judgment → next round or STOP
```

---

## 7. Stop (see PROBLEM.md §5)

**Phase 1 (done):** composition-tier unity — see `CONVERGENCE.md`.

**Phase 2 (active map):** `CAMPAIGN_MAP.md` — comprehensive dynamic Maxwell (M1+)
and R-compose close (RC1+). Mandatory path:

`CP0 → CP-M1-SPEC → CP-M1-NUM → CP-RC1-SPEC → CP-RC1-NUM → CP-U-FINAL`

Per-agent cards: `work/{TE,NE,TD,ND,TM,NM,TU}/CAMPAIGN.md`.  
Board: `work/TU/CHECKPOINT_BOARD.md`. Stamps: ADOPT | DEFER | REJECT in owner logs.

Prefer **V77-2** (dual-channel free medium) as first hard milestone.  
Terminal Phase 1: **V77-4** unified world or **V77-K** program kill.  
Terminal Phase 2: **M1 ∧ RC1** green with co-agreement, or DISPROVE.
