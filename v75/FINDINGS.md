# v75 FINDINGS — Setup, baseline, and multi-fabric design log

**Updated**: 2026-07-13  
**Docs**: `PROPOSAL.md` (architecture), `FIRST_STEPS.md` (actions), this file (record)

---

## 1. Setup (where things live)

### Code / theory baseline

| Item | Location / value |
|------|------------------|
| Kernel | `sfa/sim/scp_sim.c`, `scp_sim.cu` — multi-fabric **authorized** (2026-07-12) but **not yet used**; Step 1 ran E-lite |
| SFA | `sfa/format/sfa.h` |
| Standard U(1) package | `complex_phi=1`, `complex_gauge=1`, `g_gauge=0.05`, `m_theta=1.6`, `eta=0` |
| Potential | \(V_t(s)=(\mu/2)s/(1+\kappa s)\), \(s=\prod\|\Phi_a\|^2\), \(m^2=2.25\), \(\mu=-41.345\), \(\kappa=50\) |
| Gauged Q-ball window | \(\omega\in(1.406,1.5)\), **\(Q_{\max}=921\)** at \(g=0.05\) |
| Branch data | `v69/theory/gscan.tsv` |
| Multi-ball seeds | `bin/gen_qball_multi` (neg. \(\omega\) = opposite charge) |
| Nested grids (design) | `sfa/MULTI_RESOLUTION.md` |
| Runner | `scp-runner` MCP |
| Large data | `/space/scp/` |

### Standing goal

CLAUDE.md: **carbon atom from fabric**. Stages 0–2A nuclear done (v74);
Stage 3 = light opposite-charge + multi-scale binding → v75.

### v74 nuclear baseline (measured)

| Run | Recipe | Seed Q | Final Q (T=300) | vs \(Q_{\max}\) | Morphology |
|-----|--------|--------|-----------------|-----------------|------------|
| **c6_light** | 6× ω=1.46 octahedron D=10 | 907.8 | **649.8** | 0.71× | One droplet, mid-branch |
| **c12_light** | 12× icosahedron D=10 | 1959.9 | **1410.5** | 1.53× super | One droplet, evaporating |

### Step 1 force campaign setup (2026-07-12 → 13)

| Item | Value |
|------|--------|
| Grid | N=192, L=48, T=100, dx≈0.503 |
| Balls | light ω=±1.46 profile `f_w146_g005`, per-ball Q≈114.1, M≈173.3 |
| Separations | D=16 and D=20 (centers ±D/2 on x) |
| Runs | `pm_*` opposite (±ω), `pp_*` same (++) co-phase |
| GPU | Vast V100-SXM2-16GB (`v75step1`), ~13.5–14 min/run, ~102–105 ms/step |
| Analysis | `sfa_qball_track` (threshold 0.25) → D(t); quadratic fit \(D=D_0+v_0 t+\frac12 a_{\mathrm{rel}}t^2\), t≥20 |
| Naive Coulomb | \(a_{\mathrm{rel}}^{\mathrm{naive}}=2g^2 Q^2/(4\pi D^2 M)\) |
| Estimator | \(a_C=(a_{\mathrm{same}}-a_{\mathrm{opp}})/2\) cancels systematics (v70 method) |
| Data | `/space/scp/v75/results/` (diags, tracks, partial SFAs); tracks = analysis source of truth |
| Instance | torn down 2026-07-13 after downloads stalled/finished |

---

## 2. Findings (physics + design)

### F1 — One same-sign fabric → nuclei, not atoms [v74]

Liquid-drop fusion; A=12 always super-critical at g=0.05; Z-carbon (A=6) parks.

### F2 — Q-ball alone is not an atom [design]

Needs light opposite-charge cloud + bind without short-range fusion of core+cloud.

### F3–F6 — Multi-fabric architecture + fidelity ladder [design]

Three fabrics C/Q/L primary; L0–L3 fidelity ladder; engagement from first principles.
See `PROPOSAL.md`.

### F7 — Opposite charge is seedable without kernel change [tooling]

`gen_qball_multi` ±ω works. Step 1 = Option **E-lite**.

### F8 — Step 1 L0 force law: **PASS** [measured 2026-07-12/13]

**Question:** Do opposite-charge light Q-balls attract under Coulomb and stay as
**two** clusters over T=100?

#### Results (from track TSVs)

| Run | D₀ → D_end | ΔD | a_rel (t≥20) | a/naive | Clusters |
|-----|------------|-----|--------------|---------|----------|
| **pm D16** (±) | 16.04 → 15.43 | **−0.61** | **−1.60×10⁻⁴** | −1.37 | **2–2** |
| **pp D16** (++) | 16.04 → 15.80 | −0.23 | −2.61×10⁻⁵ | −0.22 | **2–2** |
| **pm D20** (±) | 20.02 → 19.68 | **−0.34** | **−5.20×10⁻⁵** | −0.70 | **2–2** |
| **pp D20** (++) | 20.02 → 20.22 | **+0.20** | **+4.26×10⁻⁵** | +0.57 | **2–2** |

Coulomb estimator \(a_C=(a_{\mathrm{pp}}-a_{\mathrm{pm}})/2\):

| D | a_C | naive | ratio |
|---|-----|-------|-------|
| 16 | 6.69×10⁻⁵ | 1.17×10⁻⁴ | **0.57** |
| 20 | 4.73×10⁻⁵ | 7.48×10⁻⁵ | **0.63** |

**Gauss / charge bookkeeping**

- pm: \(Q_{\mathrm{tot}}\approx 0\) (machine noise ~1e-10…1e-13); gauss_max ~1e-13  
- pp: \(Q_{\mathrm{tot}}\approx 230\)–232 ≈ 2× light Q; gauss floor held  
- Energy drift ~0.01% over T=100  

**Interpretation**

1. **Opposite charges attract** (pm: a_rel < 0 at both D).  
2. **Same charges repel at D=20** (pp a_rel > 0) — clean Coulomb sign.  
3. **D=16 same-charge still net approaches** — co-phase **contact residual** dominates
   short-range (light balls, larger relative size; known v70/EXISTENCE clock interference).  
4. **Two distinct clusters for entire T=100** on all four runs — no instant merge
   or annihilation on this timescale.  
5. Absolute \(a_C\) is ~0.6× naive vacuum Coulomb — same class of **box / sponge
   systematics** as v70 (L=48, D=16–20; wall distance not ≫ D). Force **sign**
   and relative (same vs opp) are trustworthy; exponent-from-pair-dynamics still
   needs L ≳ 2D + sponge or halo method.  
6. Consistent with v70 flo/fl campaigns (heavier ω=1.42 balls); light ω=1.46
   reproduces the same Coulomb story.

**SFA download integrity (secondary)**

| File | Local status after long rsync |
|------|-------------------------------|
| pp_force_D16.sfa | **OK** 52 frames t=0–100 |
| pm_force_D20.sfa | **OK** 52 frames t=0–100 |
| pm_force_D16.sfa | Corrupt index (repair → t=0–0.63 only) |
| pp_force_D20.sfa | Bad (24-col seed-like header) |

**Analysis does not depend on the bad SFAs** — tracks were produced on remote
from complete files. Prefer re-download only if volume renders needed for all four.

### F9 — Impact on multi-fabric process [2026-07-13 decision update]

Step 1 **does not kill** multi-fabric, and **does not make full C/Q/L urgent
tomorrow**. It **reorders** the path.

#### What Step 1 proved (E-lite = same fabric, ±ω)

| Claim | Status after Step 1 |
|-------|---------------------|
| Shared U(1) mediates long-range force | **Confirmed** (signs + two-body a_C) |
| Opposite charge can be multi-center (2 clusters) | **Confirmed** over T=100 |
| Instant annihilation / fusion of ± | **Not observed** at D=16–20, T=100 |
| Immortal bound orbit / positronium | **Not yet tested** (force at rest only) |
| e⁻ as separate light sector with private potential | **Not tested** (same \(V_t\), same mass window) |
| Nuclear bag cannot swallow leptons | **Not tested** (no N+e composite) |

#### How this changes Option B (three-fabric C/Q/L)

| Multi-fabric design element | Impact of Step 1 |
|-----------------------------|------------------|
| **Shared massless A** as engagement channel | **Validated** — keep as primary long-range bar |
| **Opposite \(g_Q g_L < 0\)** | **Validated** as the right charge assignment |
| **Private nuclear bag \(s_C\)** (L never in product potential) | **Still required for atoms**, not for ± force — same-fabric contact residual at D=16 shows short-range same-\(V\) physics is still active between co-phased cores |
| **Separate light mass window \(m_L \ll m_N\)** | **Still open** — Step 1 used nuclear-class light Q-balls (ω=1.46), not a lighter electron sector |
| **Kernel multi-fabric block (true C/Q/L split)** | **Deprioritized for next sprint** — E-lite is enough to develop Coulomb + orbit tooling; build multi-fabric when E-lite fails at orbit/atom scale or when private bag is mandatory |
| **Option A / E path** | **Strengthened** as near-term implementation: same matter multiplet + ± charge + A, before three full Cosserat copies |
| **L3 n-bar force graph** | Can now set \(V_{\mathrm{Coulomb}}(r)\) sign from L0; short-range contact term still needs a measured merge/annihilation law |

#### Revised process (what to do next)

```
[done]  F1–F8 baselines + Step 1 force (E-lite)
   │
   ├─► Step 2: ± orbit / fall dynamics (tangential boost or head-on)
   │     • If long-lived multi-rev orbit → atom path on E-lite / Option A
   │     • If capture+annihilation always → multi-fabric isolation becomes mandatory
   │
   ├─► Scale co-design: lighter L balls (higher ω, smaller Q) for m_e ≪ M_N
   │
   ├─► Only then: multi-fabric kernel (user-auth already given)
   │     Implement when Step 2 shows short-range same-V fusion/annihilation
   │     blocks hydrogenoid, OR when we need L not to feel s_C
   │
   └─► L3 n-bar after orbit+force graph measured
```

#### One-sentence process impact

> **EM selective engagement (shared A, opposite charge) is real in 3D; multi-fabric’s remaining job is isolation of short-range bag physics and mass hierarchy, not inventing Coulomb.**

---

## 3. Open nulls / risks

| Risk | Status |
|------|--------|
| ± force / two clusters | **Closed PASS** (Step 1) |
| Long-lived ± orbit (positronium) | **Open — Step 2** |
| Annihilation rate law (CONCEPT slow segregation) | Open — longer T / contact D |
| Light electron-scale Q-ball window | Open |
| Nested grid production | Design only |
| Multi-fabric kernel implementation | Authorized; **not started** — wait Step 2 gate |
| Absolute force calibration / n=2 from pairs | Need L ≳ 2D or halo (v70) |
| SFA download corruption under slow rsync | Two of four SFAs bad; tracks OK |

---

## 4. Chronology

| When | Event |
|------|-------|
| 2026-07-11 | v74 Stage 0–2A nuclear campaign |
| 2026-07-12 | Multi-fabric proposal; fidelity ladder; FIRST_STEPS |
| 2026-07-12 | Step 1 seeds + GPU force campaign (pm/pp × D16/D20) |
| 2026-07-12/13 | Track analysis: Coulomb signs; 2 clusters; estimator ~0.6× naive |
| 2026-07-13 | Download monitor timeout; 2/4 SFAs fully good; instance torn down |
| 2026-07-13 | This update: Step 1 → multi-fabric process impact |

---

## 5. Bottom line

**Nuclear (v74):** one fabric makes droplets/nuclei, not atoms.  

**Force (v75 Step 1):** opposite-charge light Q-balls **attract**, same-charge
**repel** at D=20, and **remain two centers** for T=100 — shared U(1) engagement
works in full 3D without a multi-fabric kernel.  

**Multi-fabric process:** keep C/Q/L as the **long-term architecture** for
true atoms (private bag + light sector + scale hierarchy). **Near-term work
stays E-lite:** orbits, lighter L balls, then hydrogenoid. Implement true
multi-fabric kernel when E-lite short-range physics (same \(V_t\)) blocks
stable N+e composites or when we need L outside \(s_C\).

**Next concrete step:** FIRST_STEPS Step 2 — ± dynamical orbit / capture
(not another force-at-rest).
