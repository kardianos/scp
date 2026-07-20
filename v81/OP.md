# v81 OP — Stage-3 parallel redesign (kernel authorized)

**Status:** **EXECUTED** 2026-07-19 — three parallel path agents done; see `FINDINGS.md` + `coordination/COMPARE.md`  
**Written:** 2026-07-19  
**Purpose:** After context compaction, a new session can execute this file without prior chat.  
**Post-exec note:** Kernel port **not** applied; P1 `KERNEL_DELTA.md` awaits user accept.

---

## 0. Authorization (binding)

### Kernel / format

| Scope | Authorization |
|-------|----------------|
| `sfa/sim/scp_sim.c`, `sfa/sim/scp_sim.cu` | **AUTHORIZED** for Stage-3 redesign work described here |
| `sfa/format/sfa.h` | **AUTHORIZED** only if a path needs new export columns/chunks; prefer sandbox-only SFA export first |
| Multi-fab product experiments (C/Q/L atom parks, soft \(v_t\) orbit grids, Z6+L6) | **OUT OF SCOPE** — do not resume |
| Cosserat nuclear stack (v74 Q-balls, Z-carbon) | **READ-ONLY asset** — do not break; paths may *later* couple a bag-lock to a Cosserat core |

User statement (binding): step away from multi-fab product; three independent Stage-3 redesigns; **kernel updates approved**; run three paths **in parallel in different sub-agents**.

### Standing goal

Carbon atom structural analog from fabric alone. Stages 0–2A (nuclei) stand.  
**This OP only opens Stage 3:** light opposite-charge / positronium-class bound system.  
**Not this OP:** carbon atom, multi-L shell, Stage 5 abundance.

### Evidence that closed multi-fab product

| Fact | Ref |
|------|-----|
| F PASS short rest C↔L attract | `v80/campaign/tracks/TRACK_RESULTS.md` |
| O FAIL long soft evaporates; hard shreds L | `v80/campaign/tracks/orbit2/ORBIT2_FAIL.md` |
| v79 L nulls \(E_{\mathrm{em}}\) | `v79/` |
| v80 free+locks thesis | `v80/REPRESENTATION.md`, `SHAPES.md` |
| Council consensus | `v81/council/SYNTHESIS.md` + three `RECOMMENDATION.md` |

**Verdict already recorded:** multi-fab product ladder = `STOP_AND_RECONSIDER`. Do not re-litigate.

---

## 1. Thesis constraints (all three paths must obey)

```text
REAL  = FreeSubstrate + Locks   (or pure free medium with topological locks)
POV   = Charts C / M / E as re-partitions of one state
EXPORT = optional grid/SFA projection — never the home ontology
```

**Hard kill gates (any path dies if):**

1. Primary state reduces to multiplet \(\phi(x)\) (or N copies) with no free-structure dynamics.  
2. “Particle” is only a field hump / Q-ball profile by definition.  
3. Binding uses **inserted pairwise Coulomb** / springs / damping (medium is decorative).  
4. Light sector is defined by hand-placed multiplet L / multi-fab.  
5. Free-frame \(c\) becomes GRIN gravity (v76 kill); \(c_{\mathrm{eff}}\) only as M-chart readout.  
6. First metric is carbon spectroscopy instead of **durability + ledger honesty**.

**First success metric (all paths):** durable opposite “charges” + closed energy/charge ledger on long \(T\).  
**Second metric:** emergent attraction (F-like) and/or multi-rev / capture without evaporation.

---

## 2. Three independent paths (parallel)

| ID | Directory | Name | v80 shapes | Kernel touch? |
|----|-----------|------|------------|---------------|
| **P1** | `v81/path1_locks/` | Locks on gauge medium (PIC-monist) | 4 + 10 + 6 | Sandbox first; **may** port `scp_locks` into kernel if green |
| **P2** | `v81/path2_ell/` | Path-measure \(\ell\) + locks | 3 + 4 + 10 + 6 | Sandbox first; kernel only if P1 medium reuse + \(\ell\) proven hyperbolic |
| **P3** | `v81/path3_token/` | Token / update-budget CA | 6 + 4 (pure) | **Never** into `scp_sim` in this OP — pure sandbox spike |

Council rank was P1 > P2 > P3 for atom delivery; **this OP still runs all three in parallel** for independent evidence. Do not merge implementations.

### P1 — Locks on gauge medium (primary atom candidate)

**State**

- Free medium: Yee lattice \(E,B\) (or link \(A\)) U(1), fixed free \(c\) (CFL = engine law).  
- Locks: array of structs `{id, type, q∈{±1}, E_star, p, x, footprint, …}` — **not** fields.  
- Optional later: capacity \(n_{\mathrm{free}}/n_{\mathrm{seq}}\) for sequestration rest energy.

**Step**

1. Maxwell/Yee update.  
2. Interpolate fields → lock force (Boris-class push).  
3. Deposit current with **charge-conserving** scheme (Esirkepov-class) → Gauss at machine floor.  
4. Energy ledger: \(\Delta E_\star + \Delta KE + \Delta E_{\mathrm{field}} +\) boundary flux ≈ 0.

**Forbidden:** pairwise \(1/r^2\) force; light = `radial_qball` profile; multi-fab L.

**First experiments (CPU)**

| Exp | Setup | Pass |
|-----|--------|------|
| L0 | Single lock, static | Gauss floor; exterior field sensible |
| L1 | Rest pair \(q=\pm1\), \(D\in\{8,12,16,20,24\}\) | Monotone \(a_{\mathrm{rel}}(D)\) from **medium only** |
| L2 | Soft long \(T\gtrsim 2000\) | Locks intact; \(E_{\mathrm{em}}\) does not null; ledger drift small |
| L3 | Hard \(v_t\) / orbit attempt | No shred-by-definition; ≥ few revs or honest RR inspiral |

**Kernel port (only if L0–L2 green):** new module e.g. `sfa/sim/scp_locks.h` + hooks after E-update in `scp_sim`; Cosserat nuclear path untouched; config keys `n_locks`, per-lock init. Document delta in `v81/path1_locks/KERNEL_DELTA.md` before editing.

**Deliverables**

- `v81/path1_locks/README.md` — state/step/build  
- `v81/path1_locks/src/` — standalone C (OpenMP OK)  
- `v81/path1_locks/FINDINGS.md` — scorecard L0–L3  
- Optional: `KERNEL_DELTA.md` + PR-style kernel patch if green  

### P2 — Path-measure \(\ell\) + locks

**State**

- Everything in P1 free medium + locks, **plus** scalar \(\ell(x,t)\) (path cost / capacity).  
- Locks source \(\ell\) (thicken). Propagation on \(\ell\)-weighted structure at **fixed coordinate \(c\)**.

**Step**

- \(\ell\) must obey a **hyperbolic** sourced law (e.g. \(\ddot\ell - c^2\nabla^2\ell = -S(\rho)\)).  
- **Kill if** elliptic Poisson-per-tick is required for sanity (side-Poisson collapse).  
- **Kill if** dynamics unchanged when \(\ell\) removed (decorative).  
- **Kill if** runaway GRIN-like “gravity” (v76).

**First experiments**

| Exp | Setup | Pass |
|-----|--------|------|
| E0 | One heavy lock pinned; \(\ell\) relaxes hyperbolically | Stationary non-runaway \(\ell(r)\) |
| E1 | Light test lock deflection | Deflection consistent with fixed-\(c\) budget, not free GRIN claim |
| E2 | Two-lock positronium-class with \(\ell\) on | Improvement vs pure P1 **or** honest null |

**Dependency note:** P2 may **copy or link** P1 medium/lock numerics as a library **only via explicit files under `path2_ell/`** (vendored copy preferred for isolation). Do not block P2 waiting for P1 green; if P1 incomplete, implement minimal shared medium inside P2.

**Deliverables**

- `v81/path2_ell/README.md`, `src/`, `FINDINGS.md`  
- Explicit statement: hyperbolic vs elliptic; GRIN kill check  

### P3 — Token / update-budget CA (research)

**State**

- Token density + link fluxes; hard per-tick hop cap \(c\).  
- Charge = topological circulation / vortex; mass ≈ sequestered tokens.  
- **No** Maxwell field required for MVP.

**Step**

- Conservative lattice-gas (or similar) transfers; total tokens exact.  
- Nucleate vortex–antivortex; measure interaction.

**First experiments (2D)**

| Exp | Setup | Pass |
|-----|--------|------|
| T0 | Token conservation + hop cap | Exact token total; \(c\) binds |
| T1 | Vortex–antivortex pair | Opposite attract (or clear signed force) |
| T2 | Long-lived bound gyration | Patterns do not evaporate like multiplet humps |

**Kill at two weeks** if no clean opposite interaction or patterns evaporate under own “radiation.”

**Deliverables**

- `v81/path3_token/README.md`, `src/` (Python or C), `FINDINGS.md`  
- **No** `scp_sim` edits from P3 in this OP  

---

## 3. Parallel sub-agent execution plan

### Orchestrator (parent agent)

1. Read this entire `OP.md`.  
2. Create directories if missing: `v81/path1_locks/`, `v81/path2_ell/`, `v81/path3_token/`, `v81/coordination/`.  
3. Spawn **three sub-agents in parallel** (worktree isolation preferred if available).  
4. Do **not** serialize “P1 then P2 then P3” unless a sub-agent fails to start.  
5. After all three report, write `v81/coordination/COMPARE.md` (ledger honesty, durability, F-like force, implementability, kernel readiness).  
6. Recommend **one** path for kernel integration (expected: P1 if green).  

### Sub-agent spawn template

Use the environment’s parallel agent tool (e.g. `spawn_subagent` / Task) with:

| Field | Value |
|-------|--------|
| **subagent_type** | `general-purpose` |
| **isolation** | `worktree` if supported (prevents cross-path file fights); else separate `cwd` under `v81/pathN_*/` |
| **capability** | read-write + execute (compile/run tests) |
| **background** | true (all three at once) |

### Prompt: Path 1 agent

```text
You are the P1 agent for SCP v81. Read and obey v81/OP.md sections 0–2 (P1 only).

Workspace root: repo root. Work ONLY under v81/path1_locks/ (and read-only refs:
v81/OP.md, v81/council/SYNTHESIS.md, v80/REPRESENTATION.md, v80/SHAPES.md,
v77 monist Maxwell notes if present).

Goal: implement standalone CPU sandbox: free U(1)/Yee medium + typed locks
(charge-conserving deposit, Boris push). No multi-fab. No pairwise Coulomb.
No Cosserat light L. Run L0–L2 experiments; document FINDINGS.md with kill gates.

Kernel: do NOT edit scp_sim until L0–L2 are green AND FINDINGS says GREEN for port.
If green, write KERNEL_DELTA.md proposing scp_locks hooks; implement port only if
KERNEL_DELTA is complete and tests still pass.

Deliver: build instructions, runnable binary/scripts, FINDINGS.md scorecard.
When done, write v81/path1_locks/STATUS.md = DONE | BLOCKED | KILLED with one-line reason.
```

### Prompt: Path 2 agent

```text
You are the P2 agent for SCP v81. Read and obey v81/OP.md sections 0–2 (P2 only).

Work ONLY under v81/path2_ell/. Read-only: v81/OP.md, v81/council/SYNTHESIS.md,
v80/REPRESENTATION.md, v80/SHAPES.md, v76 F1-3D path-cost notes if present.

Goal: free medium + locks + hyperbolic path-measure ℓ(x,t). Fixed free c.
Kill if elliptic-per-tick or decorative ℓ or GRIN gravity claim.
Run E0–E2; FINDINGS.md must state hyperbolic check and GRIN kill check.

May vendor a minimal copy of medium/lock code inside path2_ell/src — do not
edit path1_locks/ or scp_sim unless documenting a later merge.

Deliver: build, runs, FINDINGS.md, STATUS.md = DONE | BLOCKED | KILLED.
```

### Prompt: Path 3 agent

```text
You are the P3 agent for SCP v81. Read and obey v81/OP.md sections 0–2 (P3 only).

Work ONLY under v81/path3_token/. Read-only: v81/OP.md, v81/council/SYNTHESIS.md,
v80/SHAPES.md (shape 6).

Goal: 2D token/update-budget cellular substrate; hop cap = c; charge as
circulation/vortex. No scp_sim edits. No Maxwell required for MVP.
Run T0–T2; hard two-week kill if no opposite attract or patterns evaporate.

Deliver: build/run (Python or C), FINDINGS.md, STATUS.md = DONE | BLOCKED | KILLED.
```

### Coordination rules

| Rule | Detail |
|------|--------|
| File isolation | Each path owns only its directory |
| Shared theory | Read `OP.md` / council / v80 — do not rewrite thesis mid-flight |
| No multi-fab | Any agent proposing multi-fab L re-run is wrong; ignore |
| Commits | Prefer one commit per path when green, or orchestrator single commit after compare |
| GPU | Not required for this OP; CPU only unless a path explicitly needs GPU and user approves |
| Time box | Aim for runnable MVP + FINDINGS within one focused session per path; do not boil ocean |

### Orchestrator after-merge checklist

- [ ] All three `STATUS.md` exist  
- [ ] `coordination/COMPARE.md` ranks paths on durability, ledger, F-like force, cost  
- [ ] Kernel edit (if any) limited to green path + documented `KERNEL_DELTA.md`  
- [ ] `v81/FINDINGS.md` top-level summary  
- [ ] Update `v81/README.md` status line  

---

## 4. Explicit non-goals

- Multi-fab C/Q/L atom parks, soft/hard orbit product campaigns, Z6+L6  
- Weighted composite \(G\) scoring  
- Stage 4 carbon atom packaging  
- Rewriting CONCEPT.md as lab notebook  
- “Just add another fabric array”  
- Starting GPU Vast instances without user ask  

---

## 5. How to start (user or post-compaction agent)

When the user says to **execute OP** / **start v81 parallel paths**:

```text
1. Confirm OP.md is current.
2. Spawn three sub-agents with the prompts in §3 (parallel, isolated).
3. Monitor STATUS.md files; do not block agents on each other.
4. On completion: COMPARE.md + top-level FINDINGS.md.
5. Only then consider kernel port for the winning path.
```

**Do not start this OP as part of writing this file.**  
Writing this file ≠ execution.

---

## 6. Quick reference links

| Doc | Path |
|-----|------|
| This OP | `v81/OP.md` |
| Council synthesis | `v81/council/SYNTHESIS.md` |
| Design brief | `v81/DESIGN_BRIEF.md` |
| O FAIL | `v80/campaign/tracks/orbit2/ORBIT2_FAIL.md` |
| Scoreboard | `v80/campaign/GOAL.md` |
| Representation | `v80/REPRESENTATION.md` |
| Shapes | `v80/SHAPES.md` |
