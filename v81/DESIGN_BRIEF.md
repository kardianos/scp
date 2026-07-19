# v81 design brief — Stage-3 redesign (leave multi-fab product)

**Date:** 2026-07-19  
**Audience:** Independent design advisors (claude fable, glm-5.2, kimi k3).  
**Do not invent simulation results.** Use only this brief + cited v80 docs.

---

## Standing goal

**Carbon atom structural analog from space-time fabric alone** — Stages 0–2A (nuclei) done; **Stage 3 light opposite-charge sector** is the atom blocker; Stage 4 = nucleus + bound light charges.

---

## What just failed (product multi-fab)

Multi-fabric B1 (C/Q/L on **fixed N³**, three Cosserat multiplet copies, shared U(1)):

| Result | Detail |
|--------|--------|
| **F PASS** | Opposite C↔L rest pair: attract; \(a_{\mathrm{rel}}(D)\) monotone for \(D\in\{8..24\}\) |
| **O FAIL** | Soft \(T=400\): mild inspiral, \(E_{\mathrm{em}}\) held. Soft \(T\sim4000\): after \(t\sim1500\text{–}2000\), **\(Q\) and \(E_{\mathrm{em}}\) collapse** (evaporation). Hard \(v_t=0.5\): **L shredded** in \(\lesssim0.5\) rev |
| **v79** | Hand-placed Z6+L6: L **nulls \(E_{\mathrm{em}}\)** (radiation bath) |
| **S4** | Minimal C+L short/medium hold was real (COM window artifact on fixed-center \(Q_{\mathrm{core}}\)) |

**Interpretation:** force without durable bound product. Multi-fab on fixed grid **emulates** dual charge but does not implement monist free/bound sequestration. Further multi-fab orbit/multi-L parks are **out of scope**.

---

## v80 thesis (must inform redesign)

Primary state should **not** be “another multiplet of \(\phi(x)\) on fixed Eulerian \(N^3\)” (or N copies of that).

```text
REAL  = FreeSubstrate + Locks
POV   = Charts { C, M, E } over the same state
EXPORT = optional φ-grid / SFA / volview projections
```

Focus shapes (v80/SHAPES.md):

| Id | Idea |
|----|------|
| **3** | Path-measure / free capacity as primary medium |
| **4** | Dual substrate: free medium + **typed locks** (particle = lock, not hump) |
| **6** | \(c\) as **update/causal budget** per tick (literal step law) |
| **10** | Gauge/defects primary for charge (not first-class \(\Phi\) only) |

**Blend (working):** 4 architecture + 3 free content + 6 step law + 10 charge.

Kill gates for candidates that collapse back to multiplet-on-grid with no free-structure dynamics; or define “electron” only as hand-placed light Q-ball.

Docs: `v80/REPRESENTATION.md`, `SHAPES.md`, `DUAL_CHARTS.md`, `CONTEXT.md`.

---

## Authorization (this phase)

- **Kernel modification of `scp_sim` / `sfa.h` is authorized** by the user for Stage-3 redesign.
- **Step away from multi-fab product** as the path to atoms.
- Prefer designs that can become **numeric “real” simulation** (not pure philosophy).

---

## What we need from you

Propose **up to three independent Stage-3 redesigns** (or rank three paths) that:

1. Target **positronium / hydrogenoid first** (light opposite-charge bound system), not carbon atom.
2. **Capture the v80 free/bound + locks thesis** as far as is implementable in a real sim.
3. Specify **state vector**, **step law**, **conserved quantities**, **success metrics**, and a **minimum numeric experiment** (grid size, T, observables).
4. Say what goes **into the kernel** vs a **standalone sandbox** first.
5. Give **kill gates** and an order of implementation (week-scale, not multi-year).
6. Be honest about risk of “still just \(\phi(x)\) with extra arrays.”

**Deliverable format (markdown):**

```text
## Timeline reading (5–10 bullets)
## Three redesign options (or ranked paths)
  For each: state / step / how it maps v80 shapes / kernel delta / first experiment / kill gate
## Recommended primary path + first 2 weeks of work
## What NOT to do
```

Do **not** recommend re-running multi-fab Z6+L6 or another soft \(v_t\) orbit grid.
