# v77 — Foci, Open Problems, and Agent Mandates

**Date**: 2026-07-18  
**Depends on**: [`PROBLEM.md`](PROBLEM.md), v76 seed package

---

## 0. Shared seed (every agent)

Read before Round 1 work:

1. `v76/CONVERGENCE.md`
2. `v76/work/A/THEORETICAL_PACKAGE_v1.md`
3. `v76/APPROACHES.md` §1–2 (primitives: field, energy, mass, \(c\), warp)
4. `v77/PROBLEM.md`, `v77/RUN.md`

**Core slogan (seed hinge):**  
> Mass is field locked; \(c\) is free-field locality, constant in our frame;  
> energy is the continuum ledger; warp is what constant local \(c\) looks like around locks.

**Path-cost refinement (from human + O):** local free **interaction rate** fixed (\(=c\) in frame); near locks free channels thinner; long-range path cost requires free **response** \(\psi\), not only local free density.

---

## 1. Focus map

```text
                    ┌─────────────────┐
                    │  TU Unification │
                    │  G1 G2 G3       │
                    └────────┬────────┘
           ┌─────────────────┼─────────────────┐
           ▼                 ▼                 ▼
      ┌─────────┐       ┌─────────┐       ┌─────────┐
      │   EM    │       │ Dynamics│       │ Matter  │
      │ TE + NE │       │ TD + ND │       │ TM + NM │
      └─────────┘       └─────────┘       └─────────┘
           │                 │                 │
           └──────── seed: v76 F1-3D ──────────┘
```

---

## 2. Focus EM — Maxwell / free gauge

### Open problems

1. State Maxwell (or Maxwell-lite) as **free continuum linear laws**, not fields on a dualist stage.
2. Identify **charge** with lock–gauge coupling (ledger), consistent with free/bound budget.
3. Relate \(c = 1/\sqrt{\varepsilon\mu}\) to free locality (same \(c\) as path-cost locality).
4. Keep **sibling channels**: free capacity \(\psi\) (path cost) vs free gauge \(A\) or \(\mathbf{E},\mathbf{B}\) — shared medium, different constitutive laws.
5. Numeric: Gauss law, Coulomb \(1/r^2\), wave propagation at \(c\), vacuum control.

### Agents

| ID | Mandate |
|----|---------|
| **TE** | Equations, ontology, interface to \(\psi\), constitutive table, kill conditions |
| **NE** | Sandbox(es) under `work/NE/`; dualist EM baseline for Occam; export maps |

### Demo success (EM)

- Theory note: monist Maxwell-lite  
- Numeric: Coulomb + wave at \(c\) **or** documented monist no-go for EM channel  

---

## 3. Focus Dynamics — inertia, time, frames

### Open problems

1. Close **J5** (non-tautological \(m_{\mathrm{inertial}}\) vs \(E_\star/c^2\)) or prove renorm necessity.
2. Time-dependent free-capacity / free-gauge dynamics (beyond quasistatic F1).
3. Rods/clocks from free medium (C1 unfinished).
4. Lorentz / relativity structure from locality-\(c\) (sketch or kill overclaim).
5. Free radiation (waves of free medium) vs path-cost statics.

### Agents

| ID | Mandate |
|----|---------|
| **TD** | Dynamics laws, inertia targets, frame operationalization, kill conditions |
| **ND** | Numeric inertia/boost/wave tests; inherit v76 `sandbox_r4_inertia` lessons |

### Demo success (Dynamics)

- Theory: time-dep or inertia law statement  
- Numeric: J5 pass **or** clean form-factor/renorm result with kill of naïve equality  

---

## 4. Focus Matter — locks, multi-body, dual forces

### Open problems

1. What **is** a lock (stability, formation criteria) in monist language.
2. Multi-lock: path-cost attraction/repulsion vs Coulomb.
3. **Dual source**: same \(\rho_b\) (or related ledgers) sourcing \(\psi\) and gauge.
4. Hierarchy: EM strong at small scale, path-cost weak at large (constitutive, not hand-waved).
5. Bridge *vocabulary* to SCP Q-balls / charge — **without** claiming fixed-grid Q-balls prove monism.

### Agents

| ID | Mandate |
|----|---------|
| **TM** | Lock ontology, multi-lock force taxonomy, dual-coupling design |
| **NM** | Multi-lock numeric media; dual-source demos; simple bound pairs if feasible |

### Demo success (Matter)

- Theory: lock + dual-force package  
- Numeric: ≥2 locks with measurable free response (and/or Coulomb if EM ready)  

---

## 5. Focus Unification — world map and falsification

### Open problems

1. Write the **widely applicable idea** in ≤1 page (update each round).
2. Score whether EM + Dynamics + Matter share one core (G1–G2).
3. Maintain **program kill criteria** (when to declare V77-K).
4. Maintain **demo registry**: which demos live, dead, blocked.
5. Propose next-round themes to O via log (FOR_O).

### Agent

| ID | Mandate |
|----|---------|
| **TU** | Theory only; read all logs; `work/TU/WORLD.md`, `DEMOS.md`, `KILL_CRITERIA.md` |

TU does **not** own numeric sandboxes; may request NE/ND/NM tests via `FOR_NE` etc.

---

## 6. Cross-messages (same as v76)

- Write only own log + own work folder.
- Tag `FOR_TE`, `FOR_NE`, … `FOR_TU`, `FOR_O` in **your** log.
- Check in every round: read other logs; adopt/defer/reject `FOR_<self>`.

---

## 7. Demo registry template (TU maintains)

| Demo ID | Focus | Theory doc | Numeric | Status | Blocks |
|---------|-------|------------|---------|--------|--------|
| D-ψ-3D | seed | v76 package | v76 B | LIVE | residuals |
| D-EM-… | EM | … | … | … | … |
| … | … | … | … | … | … |
