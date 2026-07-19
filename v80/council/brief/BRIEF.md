# Shared brief — next-step recommendations for SCP (v80 era)

**Date:** 2026-07-19  
**Audience:** Independent advisors. You did **not** write the campaign. Recommend **next steps only**.  
**Do not** invent simulation results. Use only this brief.

---

## 1. Standing long-term goal

Produce a **carbon atom structural analog from the space-time fabric alone** (gauged complex Cosserat field content) — not by importing chemistry or external particle species.

Carbon is a **scale stack** (do not skip):

| Stage | Target | Status |
|-------|--------|--------|
| 0 | Carbon mapping spec (Z, A, observables via conserved bookkeeping) | Done — `v74/CARBON_MAP.md` |
| 1 | Nucleon template (gauged Q-ball) | Largely done |
| **2A** | Liquid-drop Z-carbon nucleus | **Done** — Z-carbon parks; free A=12 super-critical at g=0.05 |
| 2B | Multi-center carbon (optional) | After 2A |
| **3** | Light opposite-charge sector (electron analog); positronium first | **Blocking for atoms** |
| **4** | Atom = C-nucleus + 6 opposite light charges, Coulomb-bound | After 2+3 |
| 5 | Spontaneous production / abundance peak near carbon | After engineered seeds work |

**Success (structural):** conserved-quantity bookkeeping at machine precision; objects on Q-ball/fusion branch; closed fabric ledger where process-form applies; no fatal radiation. Parallels to real carbon are structural, not quantitative.

**Default theory era:** complex Cosserat + gauged U(1); multi-fabric C/Q/L (B2: p=C+Q, n=C-only, q_C=0,q_Q=+1,q_L=−1). Kernel (`scp_sim`) must **not** be modified without explicit user authorization.

---

## 2. Recent timeline (condensed)

### v74–v75
- Z-carbon Stage 2A done; multi-fabric architecture toward atoms designed/authorized.
- E-lite force (same fabric ±ω): opposite attract, same repel measured (F8–F9).
- Orbit/capture: partial (classical-ish at large D; head-on merges).

### v76–v77
- Monist free-capacity path-cost (F1-3D); local GRIN killed.
- Full Maxwell + dual-channel RC1 co-field (ψ ≠ EM); shared c.

### v78
- R1–R10 physics relations freeze; particle ladder + recipes.
- Campaign STOP: packages ADOPT; **ideal C₁₂ atom product BLOCKED**.
- Residual #3: long-T multi-fab GPU atom matrix.

### v79 — multi-fab Z=6 atom long-T (V100)
Primary `z6_a_n6_T800` (Z6N6+L6) vs control `z6_n6_nuc_T800` (L=0), N=192 L=48 T=800:

| Metric | Atom END | Core L=0 END |
|--------|----------|--------------|
| Q_phi | 3002 (−15%) | **same 3002** |
| Q_core / r_core | 2588 / 6.60 | **identical** |
| Q_flux | 145 (−68%) | 533 (+17%) |
| E_em | **0.36 collapse** | **~6.13 stable** |
| gauss_max | 1e-13 floor | 1e-13 floor |

**v79 verdict:** Gauge OK; **no cold atom**. L shell **nulls E_em and flux** while C/Q inventory tracks bare nucleus. Hand-placed multi-L is a neutralizing/radiation bath, not a parked multipole.

### v80 — representation thesis + overnight product campaign

**Thesis (design, not proven by GPU night):** multi-fab and fixed-grid multiplet fields are an **emulation** of free/bound continuum. v78 has numeric expression (R1–R10) but lacks a **primary simulation state** that *is* free substrate + locks (path measure, dual free+defect, c-as-update-budget, gauge defects). External kimi k3 review: grid is innocent; density-only state is the problem; recommend 2D free-medium + typed locks toy before kernel rewrite.

**Overnight V100 product campaign** (force/orbit graph, **not** Z6 re-park): 16/16 jobs OK, ~4.4 h wall, then instance torn down.

Goal function:
```
G = 0.25*S_force + 0.25*S_pair + 0.20*S_orbit
  + 0.15*S_Lhold + 0.10*S_ledger + 0.05*S_morph
```
**Final G = 0.62** (threshold ≥0.55 = continue product).

| Component | Score | Evidence |
|-----------|------:|----------|
| S_force | 0.55 | Force jobs complete; elite opp Q≈0, same Q~230; E_em(same)>E_em(opp) at D=12; **no a_rel pair tracks** (SFAs pruned) |
| S_pair | 0.65 | C–L rest T=400: Q and E_em held; head-on E_em 0.74→0.31; rest E drift +14% |
| S_orbit | 0.50 | 3 vt complete; low vt mild; **vt=0.08 E drift −27%**; multi-rev not measured |
| S_Lhold | 0.70 | S4 hydrogenoid-class N=192 T=800: **E_em held 0.71→0.75** (≠ v79 collapse); Q_phi held; Q_flux 114→4.7 and Q_core 114→0.46 (likely COM/window, not charge death) |
| S_ledger | 1.00 | Gauss floor all runs |
| S_morph | 0.25 | No volview |

**Key contrast v79 vs v80 S4:** fat Z6+L6 kills E_em; **minimal C+L hydrogenoid-class holds E_em** over long T.

**Not claimed:** measured Coulomb a_rel; multi-rev orbit; cold multi-electron atom; monist free-substrate implemented.

---

## 3. Constraints for recommendations

1. Do **not** modify `scp_sim` / `sfa.h` without explicit user auth.  
2. Prefer experiments via configs + seed generators (`gen_qball_multi`, `gen_mf_pair_boost`, `gen_mf_shell_orbit`, `gen_pn_core`, `mf_pair_track`).  
3. Large data → `/space/scp/`; keep diags local.  
4. Two parallel program lines are legitimate: **(A) product multi-fab mechanics**, **(B) representation substrate (v80 design)**.  
5. Soft-kill multi-fab atoms only if G≪0.3 with clear force death — **not** the case now.

---

## 4. Your task

Write a markdown report with these sections:

### A. Timeline reading (5–10 bullets)
What the last few versions actually established vs still open.

### B. Long-term goal status
Where Stage 3–5 stand after v79+v80. What would count as real progress toward carbon atom.

### C. Recommended next step (primary)
**One** primary next action (or a tightly ordered 2–3 step sequence if inseparable). Include:
- What to run/build
- Why this before alternatives
- Success / kill criteria
- Estimated effort (CPU toy / short GPU / overnight GPU)
- Kernel auth needed? (yes/no)

### D. Explicitly defer or kill
What **not** to do next (and why).

### E. Product vs representation
How you would split effort between multi-fab product ladder and free-substrate representation work over the next 1–2 weeks of calendar time.

### F. Confidence
Self-score 1–5 on recommendation quality given the evidence gaps (missing pair tracks, etc.).

Be blunt. Prefer concrete experiments over slogans. Do not re-run Z6+L6 park as the main next experiment unless you give a strong new mechanism.

Output **only** the markdown report (sections A–F).
