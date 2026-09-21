# v93 theory/sustain campaign — STATUS COMPLETE

**Closed:** 2026-08-09 ~07:00 CDT  
**Calling window:** 2026-08-08 → 2026-08-09 11:00 CDT (hourly loop **cancelled**)  
**Adoption judgments (theory/code/numerics for every lever):**  
→ **`FINAL_STATUS_ADOPTION.md`** (revised after Claude review  
`REVIEW_FINAL_STATUS_claude.md` — MOSTLY SOUND → fixes applied).  
**This file** = FCS paths, visualizations, scoreboard index.

---

## 0. One-line answer

**Yes — a compact condensate sustains** through T=10 000 under **unitary dense hop + law `q_detune`**, with no fuel. Radiance alone taxes (no Em\* park). Dual intake is real (quiet SLOW_REFILL) but **not adopted**. Battery **87/87 GREEN**. Shear is **not** a birth switch; noise fuel is **box-fill**, not a particle.

---

## 1. Headline scoreboard

| campaign | result | detail |
|---|---|---|
| **Sustain T=3000** (29 arms) | **DONE** | Compact HOLD: S00/S04/S14/S16 (+S00b,S10); additive/gated/qd0 fail classes |
| **R1 T=10000** (8 arms) | **6/8 PASS HOLD** | S00,S04,S14,S16,s111,amp12 flat Em~147–207; rad arms SLOW_DRAIN |
| **Mono** (k_rad ladder) | **no Em\* park** | Monotone tax vs k_rad |
| **P2 dual (bath)** | **box-fill** | Wrong apparatus for park; dual raises Em with xb |
| **P2 dual (quiet T=5k)** | **SLOW_REFILL** | d1/d0 ≈ 1.19–1.30; not a fixed point; default-off |
| **Claude qgraph** | **VERDICT YES** | Unified I/O + shells + Q=Z−N (symbolic); not law-adopted |
| **Battery** | **ALL GREEN 87** | After C Em_tot mirror + abx mask; `rad_dual=0` inert |

### Required switches for HOLD

- `amp_nat` unitary dense hop (linearized)  
- `q_detune` > 0 (law trap)  
- Quiet (no noise box-fill); workfn optional for pure-field wall  

### Rejected / negative

| claim | verdict |
|---|---|
| Shear genesis switch | **No** (G01≈G02) |
| Noise fuel = particle | **No** (box-fill Em~2k+) |
| Radiance monodispersity | **No** (SLOW_DRAIN only) |
| Dual Em\* park | **Not yet** (SLOW_REFILL) |
| Adopt `rad_dual` into V3a | **No** (default 0) |

---

## 2. Document map (read order)

| priority | path | role |
|---:|---|---|
| **1** | **`theory/FINAL_STATUS_ADOPTION.md`** | **Adoption judgments (detailed)** |
| **2** | **`theory/STATUS_COMPLETE.md`** | **this file** — FCS + viz |
| 3 | `theory/REVIEW_FINAL_STATUS_claude.md` | External review (MOSTLY SOUND) |
| 4 | `theory/CAMPAIGN_JUDGMENT.md` | Short mission answer |
| 3 | `theory/SUSTAIN_REPORT.md` | T=3000 structural scoreboard |
| 4 | `theory/R1_SCORE.md` | T=10k endurance |
| 5 | `theory/P2QUIET_T5K.md` | Quiet dual T=5k final |
| 6 | `theory/NUMERICS.md` | Full numeric trail |
| 7 | `theory/MECHANISM.md` | Isolation → mechanism |
| 8 | `theory/qgraph/docs/UNIFIED_LAW.md` | Claude symbolic law |
| 9 | `theory/qgraph/docs/VERDICT.md` | YES + queue rules |
| 10 | `theory/HOURLOG.md` | Chronological log |
| 11 | `/space/scp/v93/README.md` | Showcase FCS + WebP catalog |

---

## 3. Visualizations (WebP / PNG)

**Root:** `/space/scp/v93/`

| kind | path |
|---|---|
| Animated WebP | `/space/scp/v93/webp/*.webp` |
| Triptych stills | `/space/scp/v93/webp/*_triptych.webp` |
| Frame PNGs | `/space/scp/v93/png/<name>_em/` and `*_ee/` |
| Showcase FCS | `/space/scp/v93/fcs/*.fcs` |
| Showcase README | `/space/scp/v93/README.md` |

### Recommended films

| webp | story | source FCS (showcase) |
|---|---|---|
| `hold_S00_em.webp` | Quiet unitary HOLD (the positive result) | `fcs/hold_S00_uni_quiet.fcs` |
| `hold_S04_em.webp` | Deep trap HOLD | `fcs/hold_S04_deep.fcs` |
| `rad_tax_em.webp` / `rad_tax_ee.webp` | Radiance: Em drain, Ee rise | `fcs/rad_tax_S05.fcs` |
| `fuel_noise_em.webp` | Noise → box-fill (not particle) | `fcs/fuel_noise_S07.fcs` |
| `progressive_kx_em.webp` / `progressive_kx_phase.webp` | Phase tilt HOLD_ext | `fcs/progressive_S14_kx.fcs` |
| `genesis_wind_em.webp` | Open-door diffuse genesis | `fcs/genesis_G01_wind.fcs` |
| `collision_blob2_em.webp` | Two-blob collision | `fcs/collision_blob2.fcs` |
| `emf_absorb_noise_em.webp` | Door intake + field fuel | `fcs/emf_absorb_noise.fcs` |
| `ambient_vs_core_rad_em.webp` | Core vs ambient under tax | `fcs/ambient_vs_core_rad.fcs` |
| `core_quiet_em.webp` | Quiet control | `fcs/core_quiet.fcs` |

View with any image browser, or re-render:

```bash
# example: Em view of long-horizon HOLD
# (use campaign FCS below; showcase already has short-cadence streams)
volview -view em -px 512 /home/d/code/scp/v93/theory/runs/R1_T10000_20260808_172409/R1_S00.fcs
```

---

## 4. Campaign FCS archive (primary data)

All under `v93/theory/runs/` (symlinks `ACTIVE_*` point at latest).

### 4.1 Sustain T=3000 — 29 FCS

**Dir:** `theory/runs/sustain_T3000_20260808_104548/`  
**Score:** `theory/SUSTAIN_REPORT.md`

| stem | role |
|---|---|
| `S00_uni_base.fcs` | **Canonical quiet HOLD** |
| `S00b_uni_base_s111.fcs` | Seed-robust HOLD |
| `S04_qd36_deep.fcs` | Deep trap HOLD |
| `S14_kx_tilt.fcs` | Progressive HOLD_ext |
| `S16_seed_mw.fcs` | Winding HOLD |
| `S01_additive.fcs` | FAIL class (additive) |
| `S02_uni_gated.fcs` | STRUCT_MELT |
| `S03_qd0_melt.fcs` | DIFFUSE |
| `S05_rad.fcs` / `S06_deep_rad.fcs` | Radiance tax |
| `S07_noise.fcs` (+S08/S09/S11*) | Box-fill fuel |
| `G01`–`G05_*.fcs` | Genesis / workfn / shear |
| `R01_ring_uni.fcs` / `R02_ring_bath.fcs` | Ring |

Logs: same stem `.log`; Em series often `*.em.csv`.

### 4.2 R1 endurance T=10000 — 8 FCS

**Dir:** `theory/runs/R1_T10000_20260808_172409/`  
**Score:** `theory/R1_SCORE.md`

| FCS | Em class |
|---|---|
| `R1_S00.fcs` | HOLD ~148 |
| `R1_S04.fcs` | HOLD ~147 |
| `R1_S14.fcs` | HOLD_ext ~148 |
| `R1_S16.fcs` | HOLD ~148 |
| `R1_S00_s111.fcs` | HOLD seed-robust |
| `R1_S00_amp12.fcs` | HOLD higher mass ~207 |
| `R1_S00_rad.fcs` | SLOW_DRAIN ~105 |
| `R1_S04_rad.fcs` | SLOW_DRAIN ~114 |

### 4.3 Mono (radiance ladder) — 8 FCS

**Dir:** `theory/runs/mono_T2000_20260808_190046/`  
**Score:** `theory/MONO_FINAL.md`  
Stems: `M_qd{12,36}_kr{00,02,05,10}.fcs`

### 4.4 Dual intake P2

| suite | dir | FCS |
|---|---|---|
| Bath thermometer (box-fill) | `p2_T600_20260808_204313/` | `P2_xb*_d{0,1}.fcs` + ctrl |
| Quiet T=2000 | `p2quiet_T2000_20260808_230025/` | `Q_qd{12,36}_d{0,1}.fcs` |
| Quiet T=5000 **final** | `p2quiet_T5000_20260809_010524/` | `Q5_qd{12,36}_d{0,1}.fcs` |

Scores: `theory/P2_RESULTS.md`, `P2QUIET_SCORE.md`, `P2QUIET_T5K.md`

### 4.5 qgraph P3/P5

**Dir:** `theory/runs/qgraph_P5_20260808_175524/`  
`P3_prad{2,4}.fcs`, `P5_door{0,1}.fcs`  
Docs: `theory/qgraph/docs/P3_P5_RESULTS.md`, `P1_RESULTS.txt`

### 4.6 Early progressive probes

`theory/runs/probe_20260808_102{445,657,702}/` — shear/genesis short arms.

---

## 5. Kernel / battery artifacts

| item | path |
|---|---|
| Dual intake design | `theory/p2/DESIGN.md` |
| Pre-change baseline | `theory/p2/baseline/` |
| P2 kernel note | `theory/qgraph/docs/P2_KERNEL.md` |
| Battery GREEN log | `theory/runs/battery_p2fix2_20260809_010524/BATTERY.log` |
| Knobs (default inert) | `rad_dual=0`, `bath_xb=0`, `bath_em=0.01` |

**C + Go** both carry dual path; battery surface unchanged at defaults.

---

## 6. Symbolic package (Claude)

**Dir:** `theory/qgraph/`

| deliverable | path |
|---|---|
| Unified law | `docs/UNIFIED_LAW.md` |
| Verdict YES | `docs/VERDICT.md` |
| Predictions | `docs/SIM_PREDICTIONS.md` |
| Maxima | `maxima/*.mac` |
| Go library | `*.go` (`go test ./theory/qgraph` green) |

Not law-adopted. P1 offline: matter HOLD \(d_s\sim1.1–1.6\); foam \(\sim2.9\).

---

## 7. Next-internal (frozen recommendation)

From `CAMPAIGN_JUDGMENT.md` §5:

1. **Channel of record:** unitary dense hop + law `q_detune` (T=10k HOLD).  
2. **Do not** adopt `rad_dual` until Em\* park with restoring I′ &lt; O′.  
3. **Do not** path atoms through noise box-fill.  
4. Radiance stays far-field tax, not monodispersity engine.  
5. Genesis/coarsening is a **separate** problem from quiet sustain.

---

## 8. Quick open paths

```bash
# Judgment
less /home/d/code/scp/v93/theory/CAMPAIGN_JUDGMENT.md

# Best visualization (quiet HOLD)
xdg-open /space/scp/v93/webp/hold_S00_em.webp

# Long-horizon FCS
ls /home/d/code/scp/v93/theory/runs/R1_T10000_20260808_172409/*.fcs

# Showcase catalog
less /space/scp/v93/README.md
```

---

*Hourly campaign loop cancelled 2026-08-09 on user request. Campaign scoring complete.*
