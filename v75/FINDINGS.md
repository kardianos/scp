# v75 FINDINGS — Setup, baseline, and multi-fabric design log

**Updated**: 2026-07-13  
**Docs**: `PROPOSAL.md`, `MULTIFABRIC_SPEC.md`, `B2_PLAN.md`, `FIRST_STEPS.md`, this file

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
| Immortal bound orbit / positronium | **Partial** (F10: classical arc + scan; multi-rev open) |
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
   ├─► [done] Step 2: ± orbit PARTIAL PASS (F10) — classical arcs; head-on merges
   │     • Multi-rev still open (Step 2b); contact merge ⇒ bag isolation for N+e
   │
   ├─► Scale co-design: lighter L balls (higher ω, smaller Q) for m_e ≪ M_N
   │
   ├─► Multi-fabric kernel (user-auth already given) for overlapping N+e
   │     Not required for D~20 orbits; required when cores share volume
   │
   └─► L3 n-bar after orbit+force graph measured
```

#### One-sentence process impact

> **EM selective engagement (shared A, opposite charge) is real in 3D; multi-fabric’s remaining job is isolation of short-range bag physics and mass hierarchy, not inventing Coulomb.**

### F10 — Step 2 L0 orbit / capture: **PARTIAL PASS** [measured 2026-07-13]

**Question:** Do opposite-charge light Q-balls execute Coulomb orbital motion
without immediate bag fusion, and what happens on a head-on contact trajectory?

#### Setup

| Item | Value |
|------|--------|
| Grid | N=192, L=48, **T=400**, dx≈0.503, g=0.05, m_θ=1.6, η=0 |
| Balls | light ω=±1.46, same profile as Step 1 (track Q≈±64–65) |
| Seed tool | `bin/gen_qball_pair_boost` (±ω + independent boosts) |
| Separation | D₀=20 on x; equal-and-opposite velocities |
| Orbit seeds | **vc** vt=0.0193 (naive circular), **v05** vt=0.010 (sub), **v15** vt=0.029 (super) |
| Head-on | vr=0.01 each (closing) |
| GPU | Vast V100-SXM2-32GB (`v75step2`), ~31 min/run, ~59 ms/step |
| Analysis | remote `sfa_qball_track` thr=0.25 → D(t), φ₊(t); head-on end-state from run log |
| Data | tracks+diags `/space/scp/v75/results/`; reduced D(t) in `v75/analysis/*_Dt.tsv` |

Naive circular used \(v_{\mathrm{rel}}=\sqrt{a_{\mathrm{naive}}D}\) with
\(a_{\mathrm{naive}}=7.48\times10^{-5}\) → vt=0.0193. Measured \(a_C\) (Step 1)
gives true circular **vt≈0.0154** — so vc is slightly super-circular.

#### Orbit results (always **2 clusters**, Q≈±64.5 stable)

| Run | vt | D₀ → D_end | ΔD | φ₊ sweep | revs in T=400 | T_orb (⟨ω⟩) | Class |
|-----|-----|------------|-----|----------|---------------|-------------|--------|
| **vc** | 0.0193 | 20.02 → **20.46** | +0.44 | 0°→43.1° | **0.120** | ~3345 | **Near-circular** (mild outspiral) |
| **v05** | 0.010 | 20.02 → **15.95** | −4.07 | 0°→26.3° | 0.073 | ~5476 | **Inspiral** (sub-circular) |
| **v15** | 0.029 | 20.02 → **26.55** | +6.53 | 0°→53.9° | 0.150 | ~2671 | **Outspiral / unbound** |

Diagnostics (all three orbits): E drift **−0.015%**, gauss_max **~1e-13**,
\(Q_{\mathrm{tot}}\sim10^{-8}\), s_max flat ~0.047–0.049 — no merge signature.
Motion stays in the xy plane (cz≈0). Seed-expected ω(vc)=0.00193 vs measured
⟨ω⟩=0.00188 (−3%). Mild D growth on vc matches vt ~25% above measured-a_C circular.

#### Head-on results (contact channel)

Diag truncated at t≈256 by remote disk full; **run log completed** through T=400:

| Phase | t | Signature |
|-------|---|-----------|
| Approach | 0→250 | r_core 10.5→6.8, E_em slowly falling |
| Contact | ~300–377 | E drift grows; E_pot swings; E_em collapses |
| End state | 400 | **s_max=0.134**, r_core=6.72, **E_pot=−40.3**, E_em=0.068, E drift **−0.48%**, Q_tot≈10⁻⁶ |

**Class: CAPTURE / MERGE** — same-fabric bag wins at contact (b≈0).

#### Interpretation vs success criteria

| Criterion | Outcome |
|-----------|---------|
| Multiple revolutions, slow radiation | **Not yet** — T=400 ≈ 0.12 period (~3300). Need T≳3500 at vt≈0.0154 |
| Bound orbital arc without fusion | **YES** — vc D≈20±2%, 2 clusters full T |
| Sub / super scan classical | **YES** |
| Head-on capture | **YES** — merge at contact |
| Immediate disperse | **No** |

**Gate: PARTIAL PASS**

1. **E-lite supports Coulomb orbits at D~20** without annihilation. Multi-fabric
   is **not** required for non-contact two-body EM dynamics.
2. **Contact → merge** reconfirms private bag isolation for overlapping N+e
   (hydrogenoid), not for wide positronium-like orbits.
3. Multi-rev “immortal” positronium remains open (Step 2b).

#### Process after F10

```
[done]  Step 1 force PASS
[done]  Step 2 orbit PARTIAL PASS (arc + classical scan; head-on merge)
   │
   ├─► Step 2b (optional): multi-rev at vt≈0.0154, T≥3500
   ├─► Lighter L scout (higher ω / smaller Q)
   ├─► Multi-fabric kernel for overlapping N+e bag isolation
   └─► L3 n-bar / nested grids as needed
```

---

## 3. Open nulls / risks

| Risk | Status |
|------|--------|
| ± force / two clusters | **Closed PASS** (Step 1) |
| Coulomb orbit kinematics (sub/circ/super) | **Closed PARTIAL PASS** (Step 2 F10) |
| Multi-rev long-lived ± orbit | **Open** — need T≳ period (~3300) at vt≈0.0154 |
| Contact capture / merge law | E-lite merge (F10); **MF no-merge** (F11 G3) |
| Light electron-scale Q-ball window | Open |
| Nested grid production | Design only |
| Multi-fabric kernel implementation | **B1 G2/G3 PASS** (F11); B2 unlock next |
| Absolute force calibration / n=2 from pairs | Need L ≳ 2D or halo (v70) |
| SFA download integrity under slow rsync | Force: 2/4 bad; orbit SFAs frames=0 locally — **tracks are source of truth** |

---


### F11 — Multi-fabric B1 GPU: energy, charge law, G2/G3 **PASS** [2026-07-13]

**Question:** Does B1 Shape β (private bags, shared A, Φ_Q≡Φ_C) reproduce
E-lite Coulomb attract and **survive head-on contact without merge**?

#### Setup

| Item | Value |
|------|--------|
| Grid | N=192, L=48, g=0.05, m_θ=1.6, η=0, `n_fabrics=3`, `mf_lock_CQ=1` |
| Charges | q_C=0, q_Q=+1, q_L=−1 |
| Balls | ω_C=ω_L=**+1.46** (`f_w146_g005`), Q≈114.7 each, M≈174, D₀=20 |
| Seeds | `gen_mf_pair_boost` → dual 24-col C/L SFAs |
| GPU | Vast V100-32GB `v75mf`, q-fixed `scp_sim_mf_cuda` |
| Analysis | `mf_pair_track` (C vs L centroids) — **source of truth** |

#### Kernel fixes required for science

1. **Energy diag (GPU):** `reduce_energy` was C-only. Fixed with
   `reduce_energy_add_fabric_kernel` → E_total ≈ **347.6** (2× ball).
2. **Charge–seed rule:** E-lite ±ω under *one* charge; B1 has *fabric* charges.
   **opp ω + opp q double-flips L** → same-sign ρ_em, null relative force.
   **Rule: same-sign Noether ω on C and L.** (`gen_mf_pair_boost.c` comment.)
3. **q_em transport:** `q_em` only scaled Ampère current; matter always used
   \(U=\mathrm{e}^{i\theta}\). L with q_L=−1 **sourced − but felt +** → COM co-drift,
   ΔD≈0. **Fix:** transport \(U^q=\mathrm{e}^{i q\theta}\) in CPU+GPU force paths.

#### G2 force D=20, T=100 (rest)

| Variant | ΔD | COM | Q_flux (t=0) | Verdict |
|---------|-----|-----|--------------|---------|
| wrong_opp_omega (pre-fix) | +0.0004 | co-drift +0.18 | ~27 | null relative force |
| same_ω, no U^q (D20b) | +0.004 | co-drift −0.17 | ~0 | still null relative |
| **same_ω + U^q (D20c)** | **−0.336** | symmetric attract | ~0 | **PASS** |

D20c: D=20.000→19.664; C +10→+9.83, L −10→−9.83. **Matches E-lite PM D20
ΔD=−0.336.** E drift +1.04% (elevated vs pre-U^q ~0.01%; watch). Gauss floor.

#### G3 head-on D=20, vr=0.1, T=400

| Metric | Value |
|--------|--------|
| Wall | 51.4 min, ~97 ms/step, 41 frames |
| E | ~348.8, drift **−0.14%**, gauss ~1e-13 |
| D(t) | 20 → **D_min≈1.02** (t≈100) → **60.5** (recede) |
| mass C/L | both **~78** full run |
| Morphology | **pass-through**: centroids swap (C +10→−30, L −10→+30) |

**Class: NO MERGE** — private bags isolate short-range; opposite EM + contact
does **not** fuse. Contrast F10 E-lite head-on **MERGE**. E_em collapsed near
closest approach (~0.05 at t~100) then recovered.

#### Interpretation

| Criterion | Outcome |
|-----------|---------|
| Shared A mediates multi-fabric Coulomb | **PASS** (G2 = E-lite a_rel) |
| Two clusters at large D | **PASS** |
| Head-on contact without bag merge | **PASS** (G3) ★ |
| Energy bookkeeping with U^q | **Acceptable** on G3 (−0.14%); G2 +1% open |
| Absolute Coulomb calibration | Still ~0.6× naive (box; same as F8) |

#### Process impact

> **B1 Shape β is science-ready for atom-scale contact isolation.** E-lite remains
> the force-calibration tool; multi-fabric’s job (private bag) is **validated** at
> head-on. Unlock **B(2)** (Φ_Q free of Φ_C) or lighter L next.

Data: tracks+diags `/space/scp/v75/results/mf_*`; SFAs ~11–23G remain on `v75mf`
(not pulled — tracks are SoT).

### F12 — Hierarchy H1: heavy C@1.42 + light L@1.46 [2026-07-13]

**Question:** Under B1 lock, does a **mass hierarchy** (M_C/M_L ≈ 2.6) still
show multi-fabric Coulomb attract and **head-on no-merge**?

#### Setup

| Item | Value |
|------|--------|
| Package | B1 Shape β, q_C=0, q_Q=+1, q_L=−1, g=0.05, U^q transport |
| **C (heavy)** | ω=+1.42, `f_w142_g005`, Q≈311–315, E≈456, mass_ρ≈222 |
| **L (light)** | ω=+1.46, `f_w146_g005`, Q≈115, E≈173, mass_ρ≈78 |
| Grid | N=192, L=48, D₀=20 |
| Tool | dual-profile `gen_mf_pair_boost` |
| GPU | `v75mf` V100-32GB |

#### B-H1 force (rest, T=100)

| Metric | Value |
|--------|--------|
| Wall | 16.9 min |
| E_total | ~635 → 638 (≈ E_C+E_L), drift **+0.56%** |
| gauss | ~1e-13 |
| **D₀ → D_end** | **20.000 → 19.344** |
| **ΔD** | **−0.656** |
| Centroids | C +10→+9.83; L −10→**−9.52** (light falls in harder) |
| masses | C 222→223, L 79→79 stable |

Compare F11 equal-mass (1.46+1.46): ΔD=**−0.336**. Hierarchy |ΔD| ≈ **2.0×**
equal-mass, consistent with larger Q_C Q_L product and reduced-mass kinematics
(naive a_rel ~1.4×10⁻⁴ vs ~7×10⁻⁵).

**Class: PASS** — attract; two centers; unequal-mass kinematics correct.

#### B-H1b head-on (vr=0.08 each, T=400)

| Metric | Value |
|--------|--------|
| Wall | 51.3 min |
| E | ~637, drift **+0.006%**, gauss floor |
| D(t) | 20 → **D_min≈0.15** (t≈110) → **46.8** |
| massC / massL | **221→225 / 78→79** both intact |
| Morphology | **pass-through**; centroids swap sides |

**Class: PASS (no merge)** — light L **not** swallowed by heavy C bag at
deep contact. Same isolation story as equal-mass G3 (F11).

#### Interpretation

| Criterion | Outcome |
|-----------|---------|
| Hierarchy Coulomb attract | **PASS** |
| Light falls faster than heavy | **PASS** |
| Head-on isolation with M_C≠M_L | **PASS** ★ |
| Net EM charge (Q_C≠Q_L) | Q_flux ~38 on force (expected); not a fail |

> **First hierarchy rung works on B1.** Soft-branch min Q≈90 (ω≈1.485) is the
> next lightening step (H2). True m_e ≪ M_N still needs multi-L or package change
> beyond single-ball Q_min≈90 at g=0.05.

Data: `/space/scp/v75/results/mf_h1_*_{diag,track}.tsv`

### F13 — Hierarchy H2: soft-branch L@1.485 (Q≈90) [2026-07-13/14]

**Question:** Does the soft-branch minimum light ball (ω=1.485, Q≈90, E≈137)
still attract under Coulomb and survive head-on against heavy C@1.42?

#### Setup

| Item | Value |
|------|--------|
| C | ω=+1.42, `f_w142_g005`, Q≈315, mass_ρ≈222 |
| L | ω=+**1.485**, `f_w1485_g005` (continued from w146), Q≈**90**, mass_ρ≈60 |
| M_C/M_L | ≈ **3.3** (vs H1 ≈2.6) |
| Package | B1 lock, U^q, same as F12 |
| Profile | `v74/profiles/f_w1485_g005.txt` |

#### B-H2 force (rest, T=100) — **PASS**

| Metric | Value |
|--------|--------|
| Wall | 17.0 min |
| E | ~600, drift +0.23%, gauss floor |
| **ΔD** | **−0.616** (20.00→19.38) |
| masses | C 222→223, L 61→61 stable |
| Ql | 90 intact |

| Run | ΔD | Q_L |
|-----|-----|-----|
| F11 equal 1.46 | −0.336 | 115 |
| H1 C1.42+L1.46 | −0.656 | 115 |
| **H2 C1.42+L1.485** | **−0.616** | **90** |

Attract; slightly weaker than H1 (smaller Q_L), as expected.

#### B-H2b head-on (vr=0.08, T=400) — **PARTIAL PASS**

| Metric | Value |
|--------|--------|
| Wall | 51.0 min, 41 frames |
| E | 600→**572**, drift **−4.8%** (elevated late) |
| gauss | floor held |
| D(t) | 20 → **D_min≈0.35** (t≈110) → **40.8** |
| massC | 221→**225** intact |
| **massL** | 60→**47** (−22% by t=400); Ql 90→**72** late |
| Morphology | pass-through (centroids swap); L **sheds** after contact |

**vs H1b:** H1 light L held mass 78→79; H2 soft L loses ~22% after deep contact.
Two centers remain (not full bag merge), but soft-branch L is **less robust** under
hierarchy head-on.

| Criterion | Outcome |
|-----------|---------|
| Force attract at soft min Q | **PASS** |
| Head-on no full merge | **PASS** (two centers) |
| L mass/Q stability through contact | **FAIL soft** (−22% mass, −20% Ql) |
| Energy conservation | **MARGINAL** (−4.8%) |

> Soft-branch L works for **force** and **survives as a separate center**, but
> contact damages the light lump. Hierarchy ceiling at g=0.05 single-ball Q_min≈90
> is **marginally viable**; further lightening may need multi-L / softer kinematics
> rather than softer ω alone.

Data: `/space/scp/v75/results/mf_h2_*_{diag,track}.tsv`

### F14 — Atom ladder A1–A3: soft orbit, multi-L, Z6+L6 [2026-07-14]

**Question:** Can we build a multi-fabric “atom” — heavy nuclear C + light L cloud —
that **persists without bag merge** (toward carbon-atom / C12-scale goal)?

#### Ladder results

| Rung | Setup | Outcome | Key numbers |
|------|--------|---------|-------------|
| **A1** soft orbit | C@1.42 + L@1.46, D=20, vt hierarchy | **PASS** | D 20→17.3 inspiral; massL 78.5→78.8; no merge |
| **A2** dual-L | C origin + 2×L@±16 rest | **PASS** | massL 157→158; Ql 229 stable; D≈0 (geometry) |
| **A3** Z6 atom | C: 6×ω=1.42 octa D=10; L: 6×ω=1.46 shell R≈15.6 | **PARTIAL / atom FAIL** | see below |

#### A3 detail (T=400, wall 51 min)

| Metric | t=0 | t=400 | Notes |
|--------|-----|-------|-------|
| massC | 1638 | **1425** (−13%) | Nuclear multi-ball evaporates (v74-class) |
| massL | 474 | **308** (−**35%**) | L cloud erodes after nuclear dynamics |
| Qc | 2326 | 1931 (−17%) | Super-critical Q loss |
| Ql | 692 | 440 (−36%) | Light charge loss |
| E | 4418 | 3639 | drift **−17.6%** |
| s_max | — | ~0.06 end; **1.82 peak** @ t~126 | Nuclear fusion spike |
| D(C,L COM) | 0 | ~0 | Both sectors centered — not a merge diagnostic |

**Isolation:** L fabric mass remains **nonzero** (not fully absorbed into C bag) —
private bags still separate sectors.

**Atom persistence:** **FAIL** for “workable C12 atom.” Nuclear C does not park
cleanly (Q_max=921 at g=0.05; 6×311 seed is supercritical), and the L shell
loses ~1/3 of its mass concurrent with nuclear rearrangement. Not a stable
nucleus+cloud composite on this timescale.

#### Campaign conclusion (F11–F14)

| Claim | Status |
|-------|--------|
| Multi-fabric Coulomb attract | **PASS** (F11–F13) |
| Private bag no-merge (pair / dual-L / soft orbit) | **PASS** |
| Soft-branch ultra-light L under head-on | **PARTIAL** (F13 shed) |
| Hierarchy soft orbit | **PASS** (A1) |
| Multi-L rest around single C | **PASS** (A2) |
| **Z6 nuclear + L6 shell = persistent atom** | **FAIL** on this package |

> **Reasonable failure for C12 atom at g=0.05 with single-package balls.**
> Isolation works; **nuclear supercriticality + L erosion** block a stable atom.
> Paths forward (not yet tested): lower-g or sub-Q_max nuclear seed (true parked
> Z-carbon), larger L shell, multi-rev soft orbit only, or B2 unlock — not more
> head-on at softer ω.

Data: `/space/scp/v75/results/mf_a1_*`, `mf_a2_*`, `mf_a3_*` diags/tracks.

### F15 — Self-tune Stage 1: **DEFINITIVE SUCCESS** [2026-07-14]

**Living log:** `v75/SELF_TUNE_LOG.md` · ledger `/space/scp/v75/self_tune/`  
**VERDICT:** `/space/scp/v75/self_tune/VERDICT.md`  
**Design:** `v75/SELF_TUNE_C.md` · controller `v75/analysis/self_tune_controller.py`

**Question:** Can outer soft θ on frozen B1 action produce a scorecard-PASS
nucleus + light composite?

**Answer: YES** — first PASS at **B4_full** (single heavy C + L6 shell).

#### Winning θ (freeze)

| Param | Value |
|-------|--------|
| n_C | **1** (single heavy ball, ω_C=**1.42**) |
| n_L | **6** (ω_L=1.46, R_shell=**18**) |
| Package | B1 lock, g=0.05, q_Q=+1, q_L=−1, N=192 L=48 T=400 |

| Metric | t=0 | t=400 |
|--------|-----|-------|
| massC | 222 | 226 |
| massL | 472 | **474** (stable) |
| Qc | 315.4 | **315.4** (c_Q≈0) |
| Ql | 689 | **689** |
| E | 1505 | 1673 (+11%; c_E only) |
| gauss_max | — | **9e-14** floor |
| cost | — | **0.112 PASS** |

#### Full-trial ladder (T=400)

| ID | θ | cost | massL | Note |
|----|---|------|-------|------|
| B1_full | Z6 R18 | 0.259 | −16% | L shed + park |
| B1a_full | Z6 R22 | 0.251 | **0%** | R fixes L |
| B1b_full | Z6 softL R22 | 0.203 | −2% | best Z6 |
| B2_full | Z4 R18 | 0.161 | −1% | |
| B2a_full | Z4 R22 | 0.187 | **0%** | |
| B3_full | Z2 R18 | 0.109 | **0%** | PARTIAL c_Q=0.19 |
| B3a_full | Z2 R22 | 0.081 | **0%** | best multi-ball PARTIAL |
| **B4_full** | **Z1+L6 R18** | **0.112** | **0%** | **PASS** |

#### Claims

1. **Soft-θ packaging success** for single heavy center + multi-L shell under B1.
2. **L erosion is soft-θ fixable** (R↑) even when nuclear multi-ball parks.
3. **Multi-ball nuclei (Z2–Z6)** never hit PASS (seed→park c_Q); **single-C does**.
4. **Option C / K_st not required** for this scale claim.
5. **Caveats:** concentric rest (not orbit); E drift +11%; not multi-Z carbon nucleus;
   tree stopped early (B4a/B5 unrun). Next: orbit/bind test on frozen B4 θ;
   multi-Z needs park-aware cost or lower-g / sub-Q_max seeds.

#### Chronology note

| When | Event |
|------|-------|
| 2026-07-14 14:02 | Campaign start |
| 2026-07-14 22:43 | B4_full PASS → **DEFINITIVE_SUCCESS** |

## 4. Chronology

| When | Event |
|------|-------|
| 2026-07-11 | v74 Stage 0–2A nuclear campaign |
| 2026-07-12 | Multi-fabric proposal; fidelity ladder; FIRST_STEPS |
| 2026-07-12 | Step 1 seeds + GPU force campaign (pm/pp × D16/D20) |
| 2026-07-12/13 | Track analysis: Coulomb signs; 2 clusters; estimator ~0.6× naive |
| 2026-07-13 | Step 1 → multi-fabric process impact (F9); instance torn down |
| 2026-07-13 | Step 2: `gen_qball_pair_boost`; 4 GPU runs (vc/v05/v15/headon) T=400 |
| 2026-07-13 | Track analysis F10: near-circular arc; inspiral/outspiral; head-on merge |
| 2026-07-13 | Instance `v75step2` torn down after track/diag pull |
| 2026-07-13 | B1 MF kernel; energy L-add; q-transport U^q fix |
| 2026-07-13 | G2 D20c attract PASS (ΔD=−0.336); G3 head-on no-merge PASS |
| 2026-07-13 | F12 H1 hierarchy force+headon PASS (ΔD=−0.66; Dmin=0.15 no-merge) |
| 2026-07-14 | F13 H2 soft L@1.485: force PASS ΔD=−0.62; headon PARTIAL (L sheds 22%) |
| 2026-07-14 | F14 Atom ladder A1–A3: orbit+dual-L PASS; Z6+L6 atom FAIL (L −35%) |
| 2026-07-14 | Self-tune Stage 1 start (v75st); soft θ + ledger |
| 2026-07-14 | F15 interim mid-tree; then B4_full PASS → DEFINITIVE_SUCCESS |

| 2026-07-13 | B1 MF kernel CPU+GPU; energy L-add; G2 wrong_opp_omega null force |
| 2026-07-13 | Corrected same-ω seeds; G2 force D20b + G3 headon queue |

---

## 5. Bottom line

**Nuclear (v74):** one fabric makes droplets/nuclei, not atoms.

**Force (v75 Step 1):** opposite-charge light Q-balls **attract**, same-charge
**repel** at D=20, and **remain two centers** for T=100 — shared U(1) engagement
works in full 3D without a multi-fabric kernel.

**Orbit (v75 Step 2):** ± pairs with tangential boost follow **classical Coulomb
kinematics** (sub → inspiral, naive-circular → near-flat D with mild outspiral,
super → expand). **Two clusters** for full T=400 on all orbit runs. Head-on
**merges** at contact — same-fabric bag wins when cores collide. Multi-rev
immortality not yet shown (period ~3300 ≫ T=400).

**Multi-fabric B1 (F11):** G2 Coulomb attract + G3 **no-merge** pass-through —
private bags do the isolation job that E-lite head-on failed. Tracks are SoT;
full SFAs optional archive only.

**Atom ladder (F14): A1/A2 PASS; A3 Z6+L6 FAIL as persistent atom** (L −35%,
nuclear Q evaporates). Isolation holds; multi-Z packaging hard at g=0.05.

**Self-tune (F15): DEFINITIVE SUCCESS** — soft θ B4 (single heavy C + L6 shell
R=18) holds massC/L and Qc/Ql to T=400 under B1; scorecard PASS. Multi-ball
nuclei (Z2–Z6) never PASSed (seed→park). Option C not required for this claim.
