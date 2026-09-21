# SUSTAIN REPORT — FINAL (T=3000 matrix)

**Dir:** `theory/runs/sustain_T3000_20260808_104548`  
**Kernel:** fabrun (Go) · seed=20260802 · L=16 · dt=0.02 · bath=1 · **29/29 DONE**  
**Scored:** 2026-08-08 ~18:00 CDT (structural classes — not Em-only auto-HOLD)

---

## Mission answer

| question | answer |
|---|---|
| Does a condensate **sustain**? | **Yes** — compact HOLD under **unitary linear hop + q_detune trap**, no external fuel. |
| Which switches? | **amp_nat unitary** required; **q_detune>0** required for blob class; additive / gated-τ / qd=0 fail structure. |
| Internal vs fuel? | Quiet HOLD is **internal** (trap+U). Noise is **box-fill capacity**, not a particle. |
| Genesis? | Open door (`wf_on=0`) → **diffuse plateau** Em~640, nHot~370–390, not compact HOLD. `wf_on=1` pure field → **empty**. |
| Shear? | **Not the birth switch** (G01≈G02, ΔEm ~2%). |

---

## Structural classification (last frame + last-half slope)

Criterion (honest — Em alone is insufficient):

| class | rule of thumb |
|---|---|
| **HOLD** | Em flat ≳100; nHot≳20; Em_max≳1.0; compact support |
| **HOLD_ext** | HOLD + progressive spread (larger Rg / fd) |
| **HOLD_tiny** | Em flat; small inventory (ring topology) |
| **genesis_plateau** | birth → mid Em, nHot hundreds, diffuse |
| **box-fill** | noise-fed Em≳1500, nHot≳1000 |
| **SLOW_DRAIN / DRAIN+MELT / STRUCT_MELT / DIFFUSE / empty** | as named |

### Compact HOLD (R1 candidates)

| arm | Em | max | nHot | slope | class |
|---|---:|---:|---:|---:|---|
| **S00** uni base | 148.30 | 1.93 | 50 | ~0 | **HOLD** |
| **S00b** seed=111 | 146.85 | 1.82 | 56 | ~0 | **HOLD** (seed-robust) |
| **S04** qd=3.6 deep | 147.43 | 2.51 | 51 | −0.0003 | **HOLD** (taller peaks) |
| **S10** nowf quiet | 148.30 | 1.93 | 50 | ~0 | **HOLD** ≡S00 |
| **S14** kx progressive | 147.71 | 1.84 | 56 | ~0 | **HOLD_ext** fd=+650 |
| **S16** seed_mw=2 | 148.07 | 1.90 | 51 | ~0 | **HOLD** winding OK |

### Structure fails / drains

| arm | Em | max | nHot | class |
|---|---:|---:|---:|---|
| S01 additive | 77.6 | 0.52 | 4 | **DRAIN+MELT** |
| S02 uni gated | 148.8 | 0.34 | 0 | **STRUCT_MELT** (energy holds, peaks die) |
| S03 qd=0 | 148.8 | 0.80 | 10 | **DIFFUSE** melt |
| S05 rad quiet | 109.9 | 0.89 | 25 | **SLOW_DRAIN** |
| S06 deep+rad | 119.7 | 1.11 | 40 | **SLOW_DRAIN** (deep resists somewhat) |
| S18 V3a quiet | 109.9 | 0.89 | 25 | **SLOW_DRAIN** ≡S05 |

### Box-fill (fuel — not particle)

| arm | Em | max | nHot | note |
|---|---:|---:|---:|---|
| S07 noise | 2338 | 3.41 | 1388 | capacity plateau |
| S08 noise+rad | 2207 | 2.97 | 1553 | I≈O thruput |
| S09 noise deep | 2367 | 5.66 | 1329 | |
| S11 nowf noise | 2339 | 3.32 | 1379 | |
| S11b seed111 | 2280 | 3.56 | 1363 | seed-robust fill |
| S12 nowf noise+rad | 2206 | 3.20 | 1534 | metabolic thruput |
| S13 deep noise | 2360 | 6.03 | 1304 | |
| S15 kx noise | 2317 | 3.81 | 1364 | |
| S17 add+noise | 2157 | 2.30 | 1744 | weaker peaks |
| S19 V3a noise | 2207 | 2.97 | 1553 | ≡S08 |

### Genesis long-T

| arm | Em | max | nHot | class |
|---|---:|---:|---:|---|
| G01 nowf kx | 650 | 1.19 | 392 | **genesis_plateau diffuse** |
| G02 nowf k0 | 638 | 1.30 | 370 | ≈G01 (**+2% Em** — shear weak) |
| G03 wf kx | 0 | 0 | 0 | **empty** door wall |
| G04 nowf kx+rad | 439 | 1.12 | 74 | taxed birth −33% Em |
| G05 nowf k0+rad | 443 | 1.16 | 70 | ≈G04 (shear still weak under tax) |

### Ring

| arm | Em | max | nHot | class |
|---|---:|---:|---:|---|
| R01 ring uni | 4.88 | 0.88 | 6 | **HOLD_tiny** mono CV~0.05 |
| R02 ring bath | 4.88 | 1.54 | **1** | Em flat but **nHot 6→1** — bath dephases multi-cell ring |

---

## Isolation map (what the matrix isolates)

| switch | contrast | result |
|---|---|---|
| transport uni vs add | S00 vs S01 | unitary holds; additive drains |
| gate τ | S00 vs S02 | gated → STRUCT_MELT |
| q_detune trap | S03 / S00 / S04 | qd=0 DIFFUSE; law qd HOLD; deep taller |
| radiance outtake | S00 vs S05; S04 vs S06 | quiet tax SLOW_DRAIN; deep resists |
| field fuel (noise) | S00 vs S07 | box-fill ≠ compact HOLD |
| workfn door | S00 vs S10; G01 vs G03 | quiet blob wf-insensitive; pure field door absolute |
| progressive dense | S14, S16 | HOLD_ext / winding OK |
| genesis shear | G01 vs G02; G04 vs G05 | **weak** (~2%) |
| seed robustness | S00b, S11b | HOLD and box-fill both seed-robust |

---

## Bohmian-method note

Progressive meters (Jx, fd, support geometry) — not pilot particles. S14 shows large fd with Em flat and nHot stable → guidance real without destroying inventory. Tags lie: use Em_tot / nHot / Em_max, not Em_tag alone.

---

## Next (endurance / internalize)

1. **R1 T=10000** on S00/S04/S14/S16 (+ rad / seed / amp refinements) — see `R1_ENDURANCE.md`; live under `runs/ACTIVE_R1`.
2. Keep unitary dense channel + law qd as hold path; amp_nat adoption still needs R0 safety later.
3. Monodispersity open (blob multi-hot; ring mono but tiny).
4. Do not promote noise I as path to atoms.
5. Genesis follow-on separate: can diffuse birth coarsen to compact HOLD?

---

*Auto Em-only table retired. This report is the structural scoreboard of record for the T=3000 matrix.*
