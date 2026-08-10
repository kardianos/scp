# Sustain mechanism table (T=3000 campaign) — FINAL

**Campaign:** `runs/ACTIVE_SUSTAIN` → `sustain_T3000_20260808_104548` (**29/29 DONE**)  
**As of:** 2026-08-08 ~18:00 CDT  
**Full scoreboard:** `runs/sustain_T3000_20260808_104548/SUSTAIN_REPORT.md`

---

## Working definition

| class | meaning |
|---|---|
| **HOLD** | Em flat; nHot≳20; Em_max≳1; compact clump |
| **HOLD_ext** | HOLD with progressive spread / large fd |
| **HOLD_tiny** | Em flat; small inventory (ring topology) |
| **genesis_plateau** | birth from field → mid Em, nHot hundreds, **diffuse** |
| **box-fill** | noise-fed Em~2200+, nHot≳1000 |
| **SLOW_DRAIN / STRUCT_MELT / DIFFUSE / empty** | as before |

---

## What sustains (compact)

| arm | Em | max | nHot | note |
|---|---:|---:|---:|---|
| **S00** uni quiet | 148.3 | 1.93 | 50 | baseline HOLD |
| **S00b** seed=111 | 146.9 | 1.82 | 56 | seed-robust HOLD |
| **S04** deep qd | 147.4 | 2.51 | 51 | taller peaks |
| **S10** nowf quiet | 148.3 | 1.93 | 50 | ≡S00 |
| **S14** kx | 147.7 | 1.84 | 56 | HOLD_ext fd=+650 |
| **S16** seed_mw | 148.1 | 1.90 | 51 | winding OK |
| **R01** ring | 4.88 | 0.88 | 6 | HOLD_tiny; CV~0.05 |

**Required switches:** unitary linear hop + q_detune trap (blob class).  
**Compatible:** progressive kx (extends), seed_mw, deep trap, wf on/off for quiet blob, seed change.

## What fails structure (blob class)

| arm | class | why |
|---|---|---|
| S01 additive | DRAIN+MELT | transport algebra |
| S02 gated | STRUCT_MELT | τ kills guidance |
| S03 qd=0 | DIFFUSE | no trap |
| S05/S18 rad | SLOW_DRAIN | O-channel tax |
| S06 deep+rad | SLOW_DRAIN | deep resists somewhat |
| S07–S15/S17/S19 noise* | box-fill | open I fills volume |
| S08/S12/S19 | box-fill thruput | I≈O at filled capacity ≠ particle |
| R02 ring+bath | dephase | Em flat; nHot 6→1 |

## Genesis long-T (all five)

| arm | Em_end | max | nHot | class |
|---|---:|---:|---:|---|
| **G01** nowf kx | 650 | 1.19 | 392 | genesis_plateau diffuse |
| **G02** nowf k0 | 638 | 1.30 | 370 | ≈G01 (**+2% Em** — shear not switch) |
| **G03** wf kx | **0** | 0 | 0 | **empty** door wall |
| **G04** nowf kx+rad | 439 | 1.12 | 74 | genesis taxed ~−33% Em |
| **G05** nowf k0+rad | 443 | 1.16 | 70 | ≈G04 |

Genesis under open door produces **diffuse multi-hot matter**, not S00-class compact clump. Shear remains weak for birth mass. Radiance cuts inventory. Door (wf) is absolute for pure field.

## Headline (mission answer)

> **Yes, a condensate can sustain to T=3000** under **unitary + q_detune** without external fuel (S00/S04/S16 + seed-robust S00b).  
> Mechanism is **internal**: trap (detune) + unitary guidance — not noise metabolism.  
> Noise “sustain” is **box-fill capacity**, not a particle.  
> Progressive phase holds with extended support (Bohmian fd diagnostic).  
> Pure-field genesis needs open door and yields diffuse plateau, not quiet-blob HOLD.

## Endurance (R1 FINAL) and monodispersity

- **R1 T=10000 FULL PASS (6/8 HOLD):** S00/S04/S14/S16/s111/amp12 — Em flat, structure intact.  
- **R1 rad (2/8):** SLOW_DRAIN; Em_max drops &lt;1 — **no Em\* park** under O-only to T=10k.  
- **Mono k_rad×qd:** tax monotone; deep resists; no park.  
- **P2 dual (experimental, not adopted):**  
  - bath_xb apparatus → box-fill (wrong meter).  
  - **quiet dual T=2000:** d1 recovers (113→134) vs d0 drain — first non-flood I.  
  - **quiet dual T=5000 FINAL (both qd 1.2 & 3.6):** d0 **SLOW_DRAIN**; d1 **SLOW_REFILL**
    (Em~140 vs ~108–117; ratios 1.30 / 1.19). Deep trap resists O (+9) but does
    **not** create Em\* park; dual keeps climbing. **Not adopted.**  
  - Battery ALL GREEN 87 with dual default-off.

## Next-internal-dynamics recommendation (campaign)

See **`CAMPAIGN_JUDGMENT.md`** for the full deliverable writeup.

1. **Sustain channel = unitary dense hop + law q_detune** — validated T=3000 and T=10000; seed- and mass-robust.  
2. **Do not rely on radiance for monodispersity** — O-only is SLOW_DRAIN.  
3. **Noise / bath fuel ≠ particle** — keep classes separate.  
4. **P2 dual** closed I on quiet clump confirmed; **no Em\* park to T=5k** at either
   detune → leave default 0; redesign balance (O/I matching) before adopt.  
5. **Battery ALL GREEN 87** (p2fix2) — C≡Go restored; dual not in law table.  
6. Genesis remains a separate coarsening problem (diffuse plateau → compact?).

---

*Updated 2026-08-09 ~04:05 CDT — quiet dual T=5000 full final.*
