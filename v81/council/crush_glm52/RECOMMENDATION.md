## Timeline reading (10 bullets)

1. **Stages 0–2A are done in the Q-ball/Cosserat kernel** (`v74/CARBON_MAP.md`, `v74/RESULTS.md`): Z-carbon parks (Q→650), A=12 super-critical. That kernel is the *particle-engineering track* and stays — it is not the atom blocker and should not be rewritten.
2. **v75 multi-fab C/Q/L** was the first atom attempt: three Cosserat multiplet copies on one fixed grid sharing a U(1). It *emulates* dual charge but never had a distinct substrate type — "same kind of object × N" (v80/CONTEXT.md §3).
3. **v76 F1-3D** established the monist free-capacity math: free medium ψ, path-cost `ℓ = ℓ₀ + γψ`, mass `m = E⋆/c²`, exterior ~1/r in 3D. **KILL:** local GRIN `n(ρ_free)` as gravity. **Numeric PASS:** free-deficit ~0.13–0.17, 1/r multipole (R²~0.998). The only monist substrate with working numerics.
4. **v77** added full Maxwell (Yee) free gauge, dual-source locks (`supp|ρ_Q| ⊆ supp ρ_b`), J5-β inertia. Explicitly names the **R-compose residual:** no production binary co-evolves Yee + ψ + locks. That residual *is* the v80 gap, stated ahead of v80.
5. **v78** froze R1–R10 as a *scoring* matrix (recipes on paper), honestly noting monist ψ is optional at nuclear T. Atom product declared **BLOCKED** — not because the equations are wrong, but because no simulator represents the ontology.
6. **v79 long-T GPU atom test:** hand-placed L on multi-fab. Atom **NULL** — `Q_phi` matches L=0 control, `E_em` → 0, `Q_flux` collapses. L never opened a free-substrate sequestration chart; it only exported a third multiplet hump.
7. **v80 diagnosis:** *all* existing fields (incl. multi-fab) live on a fixed grid as primary state. The code is always `φ(x)` floats-at-sites. Multi-fab is parallel cosplay. Required primary state: `FreeSubstrate + Locks`, with C/M/E as re-partition charts over one state.
8. **v81 mandate:** leave multi-fab product; kernel mod authorized; capture the v80 thesis as far as implementable; positronium first.
9. **The honest core problem:** v76/v77 produced the monist *equations* and *sandbox* numerics; what was never built is the **production co-evolving binary** where locks are first-class discrete objects sourcing a live free medium. That is the natural 2-week target and maps ~1:1 onto the v77 R-compose residual.
10. **Dominant risk:** any "lock" re-encoded as a φ hump "to stabilize it" instantly collapses back to v79. The entire ballgame is the discipline: *locks are discrete registry objects; forces route through the live free medium; no direct lock-lock potential.* Everything below is machinery around that discipline.

## Three redesign options

Three genuinely independent paths, differing in **how matter/charge is represented**. All reuse v76 ψ and/or v77 Yee gauge as the *free medium*; they differ in the lock ontology. Ranked by near-term positronium deliverability, not philosophical purity.

### Option 1 — Lock-registry co-evolution (PIC-monist)  *(recommended primary)*

| Aspect | Specification |
|--------|---------------|
| **Shapes** | 4 (architecture) primary; 3 (path-measure) via ψ; partial 6, 10 |
| **State vector** | **Free medium on grid:** ψ(x) free-capacity scalar (v76 F1-3D); (E,B) free gauge 6-vector (v77 Yee). **No matter φ field.** **Lock registry (discrete):** `{id, species ∈ {heavy-bag, light-plus, light-minus}, x_p, v_p, ω_p, Q_p, E⋆_p, M_p = E⋆_p/c²}`. |
| **Step law (C-chart, fixed c)** | (1) Free-capacity solve each tick: `-∇·(σ∇ψ) = s·Σ ρ_b(p)` (SOR/FFT). (2) Yee Maxwell with `ρ_Q(p)`, `J_Q = Σ Q_p v_p δ(x−x_p)` (charge-conserving deposition). (3) Boris/leapfrog lock push: `F_p = q_p(E + v_p×B) + F_path(∇ψ)`. `c = 1/√(εμ)` is the literal free-gauge hop speed — first-class causal law, not just CFL. |
| **Conserved** | `Σ Q_p` (Gauss exact under charge-conserving deposition); total energy = gauge `∫(εE²+μB²)/2` + ψ + `Σ(M_p c² + ½M_p v²)`; free/bound budget `ρ_f + ρ_b = ρ_tot`. |
| **Charts** | **C**: c fixed, read ψ warp + E,B + ledgers. **M**: hold `M_p` fixed, redraw ψ as `c_eff(x)` (optical dual, *readout only*). **E**: region budget + boundary radiation flux (v79 `E_total` drift becomes honest readout). |
| **Kernel delta** | **Standalone first:** `monist_pic` (e.g. `v81/sim/monist_pic.c`) reusing v76 ψ solver + v77 Yee + lock registry + Esirkepov deposition + Boris pusher. **`scp_sim.c/.cu` NOT modified.** **`sfa.h`:** add a "lock frame" alongside the field frame (per-lock id/species/x/v/Q/E⋆), minimal delta at export. |
| **First experiment** | `N=128, L=32 (Δx=0.25), T=2000`, dt = CFL/2. Two light locks (light-plus Q=+1, light-minus Q=−1, small `E⋆`, small core) at `D=8`, small `v_rel`. Observables: `x_p(t)` orbit, `D(t)`, `E_em(t)`, ψ at cores, total charge (=0), total energy, E-chart boundary flux. Control: same locks with **gauge disabled** — must show no binding. |
| **Kill gate** | Dies if: (1) locks get re-encoded as φ humps → pivot to Option 2; (2) binding needs a direct lock-lock potential; (3) orbit collapses/evaporates <2 rev like v79 L; (4) `E_em` drifts to 0; (5) M-chart/E-chart cannot be defined on the same state (v80 kill #3). |

### Option 2 — Defect-primary gauge monism  *(deeper track / pivot target)*

| Aspect | Specification |
|--------|---------------|
| **Shapes** | 10 (gauge/defect) primary; 4 (typed locks); partial 6, 3 |
| **State vector** | **Primary:** U(1) gauge links `U_μ(x) ∈ U(1)` (or Yee E,B as Lie algebra). *This is the free medium; no matter φ.* **Defect registry:** `{id, type ∈ {plus-vortex, minus-vortex, bag-composite}, x_p, winding w_p ∈ ℤ, r_c, E⋆_p}`. ψ either read off gauge stress or co-evolved (F1-3D with `ρ_b` = defect energy density). |
| **Step law** | Free gauge Yee/wave. Defects transported under Lorentz-like force from gauge gradients at the core; winding topologically conserved. **Gauss: `∇·E = Σ w_p δ(x−x_p)`** — charge *is* winding. Positronium = +1 and −1 vortex bound by mutual gauge flux (shape 10 fully realized). |
| **Conserved** | `Σ w_p` (topological charge); Gauss exact; energy (gauge + core). |
| **Charts** | C/M/E on gauge stress + defect ledgers (cleanest of the three for v80 kill #3). |
| **Kernel delta** | New standalone binary, lattice-gauge + defect tracker; does **not** reuse φ infra. **Hard R&D:** defect-core regularization (Villain / Dirac string / compact core) — the schedule risk. |
| **First experiment** | `N=64, L=16, T=2000`. +1 and −1 vortex pair at `D=4`, small `v_rel`. Pure gauge/defect — no φ anywhere. PASS = stable bound pair. |
| **Kill gate** | Dies if: (1) defect core collapses to a φ-hump; (2) ± defects don't bind without a hand potential; (3) winding not conserved numerically. **Deferred** (not killed) if not reachable in 6 weeks — the right long-track answer to "still just φ(x)." |

### Option 3 — Path-measure / causal-budget monism  *(research note)*

| Aspect | Specification |
|--------|---------------|
| **Shapes** | 6 + 3 primary; 4 partial; 10 weak |
| **State vector** | **Primary:** one scalar `b(x)` = per-site update budget / travel-time / path-cost (shape 6's "coin" literal). Optional graph variant (shape 1). **Lock registry:** `{id, type, x_p, sequestered-budget B_p, sign s_p ∈ {+,−,0}, E⋆_p}`. No separate gauge: charge sign couples to vorticity of b-current. |
| **Step law** | `b` evolves hyperbolically at c (literal causal budget/tick). Locks sequester budget (`b` dips at lock) and move under `∇b`. ± signs source/sink b-current vorticity; opposites attract via measure torsion. `c` = max budget hop/tick. |
| **Conserved** | `∫b + Σ B_p`; charge-sign sum; energy (budget + sequestered). |
| **Charts** | C/M/E on b + ledgers. Shape-10 coupling is the weak link. |
| **Kernel delta** | New standalone binary. Cheapest substrate, but charge coupling is speculative. |
| **First experiment** | `N=64, T=1000`. Two opposite-sign locks in pure `b`; observe attraction through budget field alone. |
| **Kill gate** | Dies if: (1) `b(x)` behaves like a scalar φ with no free/bound or causal-budget semantics; (2) charge coupling must be a hand potential; (3) too metaphorical for clean conservation laws (v80 kill #5). |

## Recommended primary path + first 2 weeks of work

**Primary: Option 1 (lock-registry co-evolution / PIC-monist).**

**Why:** The only path that delivers a numeric positronium on a 2-week schedule (PIC numerics are textbook; reuses v76 ψ + v77 Yee). It directly closes the **v77 R-compose residual** — the exact gap v80 names. It satisfies the v80 kill gates *if disciplined*: locks first-class discrete (not φ humps), forces route through the live medium (not direct potentials), C/M/E charts are re-partitions of one state. Its own kill gate #1 (collapses to φ-hump) is the natural trigger to pivot to Option 2.

**Honest residual risk:** ψ and (E,B) are still fields on a fixed grid. The defense is v80/CONTEXT.md §49 — fixed grid is permitted as a *view*; the kill is when primary state *reduces to φ(x) with no free-structure dynamics*. Here ψ is the **free** medium (path-cost ontology, F1-3D equation), not matter; matter is the discrete lock registry. Full shape 6 (literal causal coin) and full shape 10 (topological defect charge) are *not* captured by Option 1 — they are what Options 3 and 2 are for.

**Two-week plan (week-scale, kernel mod deferred):**

*Week 1 — standalone `monist_pic`, no `scp_sim`/`sfa.h` changes:*
- **D1–D2:** Extract v76 ψ solver + v77 Yee into a standalone C file. Validate each against known v76/v77 results (free-deficit ~0.13–0.17; vacuum wave at c).
- **D3–D4:** Add lock registry, charge-conserving deposition (Esirkepov), Boris pusher. Single-lock-in-box: ψ well forms, Gauss holds to ~1e-12.
- **D5:** Two same-sign locks → confirm repulsion *through the gauge*. Two opposite-sign → confirm attraction (the F-PASS analog, now monist). **Gate to Week 2** — triage here if it fails.
- **D6–D7:** Chart readouts — M-chart (`c_eff(x)` at fixed `M_p`), E-chart (region budget + boundary flux). SFA export (field frame + lock frame).

*Week 2 — positronium + decision:*
- **D8–D9:** Tune light-lock species (small `E⋆`, small core) so orbit resolves on grid. ± locks at `D=8`, small `v_rel`. Run `T=2000`.
- **D10–D11:** Analyze orbit graph, `E_em(t)`, energy drift on E-chart boundary. Compare to gauge-disabled control.
- **D12:** Hard-`v_t` stress test (the v79 shred case). Does the light lock survive? Diagnose core resolution / deposition noise / budget leak.
- **D13:** **Decision.** PASS → `v81/FINDINGS.md` + Stage-4 atom path. FAIL → triage vs kill gates; if #1 or #3 triggered, open Option 2 pivot stub (4-week follow-on).
- **D14:** Writeup + SFA archive.

## What NOT to do

- **Do NOT re-run multi-fab Z6+L6 or another soft-`v_t` orbit grid** (brief-explicit; v79 proved the recipe null).
- **Do NOT represent locks as φ humps "to stabilize them."** Kill gate #1; the exact v79 failure mode. Locks are discrete registry objects, full stop.
- **Do NOT use a direct lock-lock Coulomb potential.** Force MUST route through the live free gauge (Yee). A direct potential is death by kill #2 and makes the medium decorative.
- **Do NOT modify `scp_sim.c/.cu` or `sfa.h` for Option 1.** Standalone until positronium is shown. Kernel-mod authorization is for the Option 2 pivot or Stage-4 merge — not Week 1.
- **Do NOT promote local GRIN / `c_eff(x)` as gravity.** v76 kill holds. `c_eff` is an **M-chart readout**, not the free-frame law.
- **Do NOT target carbon atom directly.** Positronium first (brief-explicit). Carbon = Stage 4.
- **Do NOT conflate "lock species" with "multi-fab fabric."** One free medium, typed locks. Three Cosserat copies is the thing we are leaving.
- **Do NOT skip the chart readouts.** M-chart and E-chart must be demonstrable D6–D7 — otherwise the v80 thesis isn't captured and the work collapses to "PIC with extra bookkeeping."
- **Do NOT let `E_em` drift silently to 0.** That is the v79 signature; instrument E-chart boundary flux from day one.
- **Do NOT block positronium on Option 3.** Most literal read of 6+3 but most metaphorical; keep it as a research note, not a schedule item.
