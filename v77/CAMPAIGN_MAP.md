# v77 Phase 2 — Campaign Map (Fullness: Maxwell + R-compose)

**Date:** 2026-07-19  
**Status:** planning map (no runs required to adopt)  
**Parent stop (A) composition-tier remains valid history** — Phase 2 does **not** re-open v76 path-cost  
**Seed:** `CONVERGENCE.md`, `work/TU/UNIFIED_PACKAGE_v1.md`, TE full Maxwell freeze  
**Why not v78:** same monist world and agent graph; this is depth of examination, not a new ontology era. v78 only if Phase 2 closes or program-kills.

---

## 0. Phase 2 terminal goals

| ID | Goal | Done when |
|----|------|-----------|
| **P2-M** | Comprehensive **dynamic** full Maxwell | Layers **M1** mandatory; **M2** preferred; **M3** stretch |
| **P2-RC** | Close **R-compose** | Tier **RC1** mandatory; **RC2** preferred; **RC3** stretch |
| **P2-U** | Cohesive agreement on the **joint** model | TU declares composition residual **closed at RC1+** (or DISPROVE) |

**STOP Phase 2** when: **M1 ∧ RC1** green with co-agreement stamps, **or** documented DISPROVE/NULL on a hard kill.

Optional stretch before any later v78: **M2 ∧ RC2**.

---

## 1. Layer map (from fullness audit)

### Maxwell layers

| Layer | Name | Content | Owner pair |
|-------|------|---------|------------|
| **M0** | Baseline (done) | Theory freeze; 2D Yee skeleton; 1D TEM wave; div B; static Coulomb recovery | TE / NE |
| **M1** | True 2D dynamic Maxwell | 2D beams/packets; energy/Poynting; Gauss+sources dynamic; radiation; PML/box; incomplete-Ampère adversary | TE gates · **NE** implement |
| **M2** | 3D dynamic Maxwell | 3D Yee; Coulomb+radiation same solver; multipoles | TE scope · **NE** |
| **M3** | Maxwell + free charged matter | Self-consistent \(\rho_Q,J_Q\); Lorentz; scatter/orbit lite | TE/TM laws · **NE** (+ ND inertia) |

### R-compose tiers

| Tier | Name | Content | Owner pair |
|------|------|---------|------------|
| **RC0** | Composition (done) | Separate sandboxes + shared ontology | NM / NE / TU |
| **RC1** | Co-field | **One grid:** dynamical Maxwell + \(\psi\); **fixed** multi-locks | TM design · **NM** lead · NE EM module · ND \(c\) |
| **RC2** | Co-evolution | Locks **move** under \(F^\psi+F^{\mathrm{EM}}\); Cont + energy | TM/TD · **NM** · NE · ND |
| **RC3** | 3D co-evolution | RC2 in 3D | All numeric · TU score |

### Dependency (do not skip)

```text
M0 (done)
  └─► M1  ──mandatory──►  RC1  ──mandatory──►  Phase-2 STOP eligible
         │                  │
         ├─► M2 (preferred) ├─► RC2 (preferred)
         └─► M3 (stretch)   └─► RC3 (stretch, needs M2)
```

**Rule:** RC2 should not start until **M1** is green (avoid debugging lock motion on thin Maxwell).  
**Rule:** RC3 should not start until **M2** is green (2D log vs 3D \(1/r\) trap).

---

## 2. Co-coordination protocol

### 2.1 Agreement stamps

When a checkpoint requires co-agreement, **each listed agent** appends to **their own log**:

```text
===== ENTRY ... | <ID>-... | checkin =====
Tags: stamp, <CHECKPOINT_ID>
---
STAMP <CHECKPOINT_ID>: ADOPT | DEFER | REJECT
Evidence: <paths / gate names>
Blockers: <none | list>
FOR_<partner>: ...
===== END =====
```

TU maintains a live table in `work/TU/CHECKPOINT_BOARD.md` (create on first Phase-2 round).  
O records round themes and STOP/CONTINUE in `logs/O_orchestrator.log`.

### 2.2 Stamp vocabulary

| Vote | Meaning |
|------|---------|
| **ADOPT** | Gate/docs meet bar; partner may proceed downstream |
| **DEFER** | Almost; need listed fix before ADOPT |
| **REJECT** | Fails monism, shared \(c\), or kill-gate; redesign |

Downstream work that depends on a REJECT/DEFER checkpoint is **out of order** (log as `blocker`).

### 2.3 Cross-read minimum each round

| Agent | Must read |
|-------|-----------|
| Every agent | O latest; TU board; **partner** log (T↔N within focus) |
| NE | TE gates; ND shared-\(c\) notes when changing \(c\) |
| NM | TE dual-channel; NE Maxwell module API; TM dual_channel |
| ND | TD J5-β; NE wave \(c\); NM force mass tags |
| TU | All focus handoffs; update board |

---

## 3. Checkpoint register (global)

| ID | Name | Depends on | Co-agreement required | Downstream unlocks |
|----|------|------------|------------------------|--------------------|
| **CP0** | Phase-2 kickoff | CONVERGENCE (A) | O + TU | All Phase-2 work |
| **CP-M1-SPEC** | M1 gate list frozen | M0 | **TE + NE** (+ TU note) | NE M1 implementation |
| **CP-M1-NUM** | M1 numeric green | CP-M1-SPEC | **NE + TE** | RC1, M2, M3 |
| **CP-RC1-SPEC** | RC1 joint state/API | CP-M1-NUM (or provisional M1-lite only if O allows) | **TM + TE + NM + NE** | NM RC1 code |
| **CP-RC1-NUM** | RC1 co-field green | CP-RC1-SPEC, CP-M1-NUM | **NM + NE + TM + TE** (+ TU) | Phase-2 STOP eligible; RC2 |
| **CP-M2-SPEC** | 3D Maxwell scope | CP-M1-NUM | TE + NE | NE M2 |
| **CP-M2-NUM** | 3D Maxwell green | CP-M2-SPEC | NE + TE | RC3 |
| **CP-RC2-SPEC** | Moving-lock laws | CP-RC1-NUM | **TM + TD + NM + ND** | NM/ND RC2 |
| **CP-RC2-NUM** | Co-evolution green | CP-RC2-SPEC, CP-M1-NUM | **NM + ND + TM + TD** (+ NE force fields) | strong unity |
| **CP-M3** | Charged matter Maxwell | CP-M1-NUM | TE + TM + NE | optional |
| **CP-U-FINAL** | Phase-2 unity / kill | CP-RC1-NUM min | **TU + O** (+ all ADOPT or explicit residual) | STOP Phase 2 |

**Mandatory path to Phase-2 STOP:**  
`CP0 → CP-M1-SPEC → CP-M1-NUM → CP-RC1-SPEC → CP-RC1-NUM → CP-U-FINAL`

---

## 4. Per-agent campaign maps

---

### 4.1 TE — Theory EM

**Mission:** Keep Maxwell monist, complete, and **gate-authoritative** for deeper exams.

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| T0 | Maintain freeze `full_maxwell_monist_v0.md` (no ontology fork) | CP0 |
| T1 | `work/TE/m1_gates_v0.md` — M1 kill-gates (2D beam, energy, dynamic Gauss, radiation, Ampère adversary, PML) | **CP-M1-SPEC** with NE |
| T2 | Stamp NE M1 results; amend gates only if kill | **CP-M1-NUM** |
| T3 | RC1 interface note: how dynamical \(\mathbf{E},\mathbf{B}\) attach to dual-source locks (no \(\psi\equiv\Phi\)) | **CP-RC1-SPEC** |
| T4 | Stamp RC1 EM channel | **CP-RC1-NUM** |
| T5 | Optional M2/M3 theory scope notes | CP-M2-SPEC / CP-M3 |

**Co-agree with:** NE (all M*), TM/NM (RC1 dual-source), TU (WORLD consistency).  
**Must not:** Reopen local-GRIN or hand-Poisson-as-monism.

---

### 4.2 NE — Numeric EM

**Mission:** Comprehensive **dynamic** full Maxwell numerics (close thin 1D-only wave story).

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| N0 | Keep `sandbox_full_maxwell_r2.py` as M0 regression | CP0 |
| N1 | Implement M1 suite under `work/NE/sandbox_m1_*.py` | after CP-M1-SPEC |
| N2 | Export `outputs/m1_result.json`; all M1 gates | **CP-M1-NUM** with TE |
| N3 | Expose **module API** for NM: step Maxwell given \(\rho_Q,J_Q\); read \(\mathbf{E},\mathbf{B}\) | **CP-RC1-SPEC** |
| N4 | Support NM RC1 integration; stamp EM half of co-field | **CP-RC1-NUM** |
| N5 | Optional M2 3D Yee + multipoles | CP-M2-* |
| N6 | Optional M3 charged sources self-consistent | CP-M3 |

**M1 gate minimum (coordinate exact list with TE at CP-M1-SPEC):**

1. Vacuum  
2. **2D** wave packet / beam \(v\approx c\) (not only 1D TEM)  
3. Off-unit \(\varepsilon,\mu\)  
4. Energy / Poynting vacuum  
5. Dynamic Gauss residual with prescribed Cont-safe sources  
6. div B  
7. Faraday discrete  
8. Incomplete-Ampère adversary fails free wave  
9. Absorbing or large-box honesty note  

**Co-agree with:** TE (gates), NM (API for RC1), ND (shared \(c\) numbers).

---

### 4.3 TD — Theory Dynamics

**Mission:** Time-dep free capacity + inertia law for co-evolution; no SR overclaim.

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| D0 | J5-β remains default | CP0 |
| D1 | `work/TD/psi_time_dep_rc1_v0.md` — T0/T1/T2 choice for RC1 (quasistatic F1 vs relaxational) | **CP-RC1-SPEC** |
| D2 | `work/TD/lock_motion_rc2_v0.md` — \(m\) in \(F=ma\) under J5-β; energy | **CP-RC2-SPEC** |
| D3 | Stamp ND numerics | CP-RC1-NUM (shared \(c\)); CP-RC2-NUM |

**Co-agree with:** ND, TM (forces), NM (motion).  
**Critical co-agree:** At **CP-RC2-SPEC**, TD+TM+ND must agree which mass enters lock acceleration.

---

### 4.4 ND — Numeric Dynamics

**Mission:** Shared \(c\), free-capacity dynamics, inertia tags for moving locks.

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| Dyn0 | Keep dual-channel-\(c\) + J5-β exports | CP0 |
| Dyn1 | Free-capacity time-dep probe consistent with TD T1/T2 | after CP-RC1-SPEC |
| Dyn2 | Shared-\(c\) regression when NE M1 changes media | **CP-M1-NUM** (observe); **CP-RC1-NUM** (stamp) |
| Dyn3 | RC2: non-tautological mass measurement under motion | **CP-RC2-NUM** |

**Co-agree with:** TD, NE (\(c\)), NM (forces on locks).

---

### 4.5 TM — Theory Matter

**Mission:** Dual-source locks and force laws for co-field / co-evolution.

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| Mat0 | dual_channel_v0 remains baseline | CP0 |
| Mat1 | `work/TM/rc1_joint_state_v0.md` — full state vector; Supp; \(\psi\) solver; Maxwell dynamical (not Φ-only) | **CP-RC1-SPEC** with TE, NE, NM |
| Mat2 | Force law RC1: \(F^\psi + q\mathbf{E}\) (fixed locks: force diagnostics only) | CP-RC1-SPEC |
| Mat3 | RC2 motion + hierarchy + Lorentz tier B | **CP-RC2-SPEC** |
| Mat4 | Stamp NM results | CP-RC1-NUM, CP-RC2-NUM |

**Co-agree with:** TE (no channel collapse), NM (implementability), TD (inertia).

---

### 4.6 NM — Numeric Matter

**Mission:** **Lead RC1/RC2** implementation — where R-compose dies or lives.

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| MatN0 | Keep dual-channel R2 (static lite) as regression | CP0 |
| MatN1 | After CP-RC1-SPEC: `sandbox_rc1_cofield.py` — fixed multi-locks; F1 \(\psi\) + **NE dynamical Maxwell** (not Φ-only) | — |
| MatN2 | Gates: neutral / same / opposite; shared \(c\); energy split; vacuum; sibling independence | **CP-RC1-NUM** |
| MatN3 | Optional RC2 moving locks | CP-RC2-NUM |
| MatN4 | Export unified JSON for TU | each CP-*-NUM |

**Co-agree with:** NE (must consume Maxwell step API), TM, TE, TU.  
**Hard rule:** RC1 **REJECT** if still only Maxwell-lite \(\Phi\) without dynamical \(\mathbf{E},\mathbf{B}\) path (unless O explicitly weakens bar — default: no).

---

### 4.7 TU — Theory Unification

**Mission:** Board, WORLD, kill criteria, Phase-2 STOP recommendation.

| Stage | Deliverable | Checkpoint |
|-------|-------------|------------|
| U0 | `work/TU/CHECKPOINT_BOARD.md` initialized from this map | **CP0** |
| U1 | Update WORLD if M1/RC1 change constitutive story (unlikely) | as needed |
| U2 | After each CP-*-NUM: board + DEMOS | continuous |
| U3 | `work/TU/PHASE2_SCORE.md` when CP-RC1-NUM green | **CP-U-FINAL** |
| U4 | Recommend STOP Phase 2 (A-depth) or DISPROVE or CONTINUE stretch | CP-U-FINAL |

**Co-agree with:** everyone at CP-U-FINAL; O final call.

---

### 4.8 O — Orchestrator

| Stage | Action |
|-------|--------|
| CP0 | This map + spawn Phase-2 Round 1 theme **M1-SPEC** |
| Each round | Read stamps; unblock or hold; no agent writes O log |
| After CP-RC1-NUM | Prefer STOP Phase 2 or authorize RC2/M2 stretch |
| DISPROVE | Only if hard kill in KILL_CRITERIA met |

---

## 5. Round plan (suggested sequence)

| Round | Theme | Primary agents | Checkpoint target |
|-------|--------|----------------|-------------------|
| **P2-R1** | Freeze M1 gates + NE design | TE, NE, TU | **CP-M1-SPEC** |
| **P2-R2** | Implement M1 numerics | NE, TE, ND (\(c\)) | **CP-M1-NUM** |
| **P2-R3** | RC1 joint state + Maxwell API | TM, TE, NE, NM, TD | **CP-RC1-SPEC** |
| **P2-R4** | RC1 co-field sandbox | NM, NE, TM, TE, TU | **CP-RC1-NUM** |
| **P2-R5** | Score Phase-2 STOP | TU, O, all stamp | **CP-U-FINAL** |
| **P2-R6+** | Optional M2 / RC2 | as authorized | CP-M2 / CP-RC2 |

Agents may parallelize **within** a round only on tasks that do not skip a REJECT dependency.

---

## 6. Co-agreement matrix (who must stamp what)

| Checkpoint | TE | NE | TD | ND | TM | NM | TU | O |
|------------|----|----|----|----|----|----|----|---|
| CP0 | · | · | · | · | · | · | ✓ | ✓ |
| CP-M1-SPEC | **✓** | **✓** | · | · | · | · | note | · |
| CP-M1-NUM | **✓** | **✓** | · | note | · | · | note | · |
| CP-RC1-SPEC | **✓** | **✓** | note | · | **✓** | **✓** | note | · |
| CP-RC1-NUM | **✓** | **✓** | · | note | **✓** | **✓** | **✓** | · |
| CP-RC2-SPEC | · | note | **✓** | **✓** | **✓** | **✓** | note | · |
| CP-RC2-NUM | · | note | **✓** | **✓** | **✓** | **✓** | **✓** | · |
| CP-U-FINAL | note | note | note | note | note | note | **✓** | **✓** |

✓ = ADOPT/DEFER/REJECT stamp required · note = optional checkin · · = no stamp

---

## 7. Kill / NULL triggers (Phase 2)

From `work/TU/KILL_CRITERIA.md` plus:

| ID | Trigger |
|----|---------|
| **K-M1** | Cannot get 2D dynamical Maxwell energy+wave without dualist stage that breaks free/bound |
| **K-RC1** | Co-field forces \(\psi\equiv\Phi\) or shared \(c\) fails on joint grid |
| **K-RC2** | Moving locks require non-monist inertia that contradicts J5-β and cannot be renormed honestly |

If K-* fires with multi-agent ADOPT of kill: TU proposes **DISPROVE**; O confirms.

---

## 8. Explicit non-goals (Phase 2)

- Real \(G\), \(\alpha_{\mathrm{EM}}\), particle data fits  
- Full GR / BH / galactic DM  
- scp_sim / sfa.h edits without human OK  
- Replacing v76 F1-3D path-cost success  
- Declaring Phase-2 done on M0+RC0 alone (that is already CONVERGENCE composition-tier)

---

## 9. File pointers agents should create/update

| Agent | Campaign artifacts |
|-------|-------------------|
| TE | `m1_gates_v0.md`, stamps in TE log |
| NE | `sandbox_m1_*.py`, `outputs/m1_result.json`, Maxwell **API** for NM |
| TD | `psi_time_dep_rc1_v0.md`, `lock_motion_rc2_v0.md` |
| ND | `outputs/` shared-c + RC2 inertia |
| TM | `rc1_joint_state_v0.md`, RC2 force/motion |
| NM | `sandbox_rc1_cofield.py`, `outputs/rc1_*.json` |
| TU | `CHECKPOINT_BOARD.md`, `PHASE2_SCORE.md`, DEMOS updates |
| O | `logs/O_orchestrator.log` only |

---

## 10. One-page agent cheat sheet

```text
TE:  freeze Maxwell → write M1 gates → stamp NE → RC1 EM interface → stamp RC1
NE:  implement M1 → export API → help NM RC1 → optional M2/M3
TD:  ψ time-dep for RC1 → mass law for RC2 → stamp ND
ND:  guard shared c → RC2 inertia numerics
TM:  RC1 joint state+forces → RC2 motion laws → stamp NM
NM:  BUILD RC1 co-field (lead) → optional RC2 move locks
TU:  board + score + STOP/KILL recommend
O:   sequence rounds; final STOP
```

**First action when Phase-2 starts:** O logs CP0; spawn **P2-R1** (TE+NE+TU primary) targeting **CP-M1-SPEC**.
