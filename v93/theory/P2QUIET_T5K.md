# P2 quiet dual T=5000 — FULL FINAL 4/4

**Dir:** `theory/runs/p2quiet_T5000_20260809_010524` (`.done` 2026-08-09 04:02 CDT)  
**Setup:** unitary blob (`amp_nat=2`), k_rad=0.05, p_rad=4, **bath_xb=0**, noise=0,
wf_on=1, L=16, T=5000, seed=20260802; rad_dual∈{0,1} × q_detune∈{1.2, 3.6}  
**Scored:** 2026-08-09 ~04:05 CDT

## Scoreboard

| arm | Em path (0→1k→2k→3k→4k→5k) | Em_end | slope 3k→5k | Em_tag_end | class |
|---|---|---:|---:|---:|---|
| Q5_qd12_d0 | 149→121→113→110→109→**108** | **108.1** | −0.0009 | 29.5 | **SLOW_DRAIN** |
| Q5_qd12_d1 | 149→125→134→136→140→**140** | **140.4** | **+0.0023** | 35.0 | **SLOW_REFILL** |
| Q5_qd36_d0 | 149→129→122→120→118→**117** | **117.2** | −0.0012 | 44.7 | **SLOW_DRAIN** |
| Q5_qd36_d1 | 149→131→129→133→138→**139** | **139.1** | **+0.0030** | 60.2 | **SLOW_REFILL** |

| comparison | value |
|---|---:|
| qd12 d1/d0 Em ratio | **1.30** (Δ +32) |
| qd36 d1/d0 Em ratio | **1.19** (Δ +22) |
| deep d0 resist vs law d0 | +9 Em (117 vs 108) — trap slows O somewhat |
| deep d1 vs law d1 | ~same Em_end (~139–140) — I still wins |

## Reading

1. **O-only (d0):** both traps keep **SLOW_DRAIN** through T=5k; slopes small and
   still negative — no Em\* floor in window. Deep trap retains ~+9 Em vs law.
2. **Dual (d1):** both traps **recover then SLOW_REFILL** — no park, no box-fill.
   After mid-horizon recovery, I > O and inventory keeps climbing slowly.
3. **Deep trap does not fix the park:** qd=3.6 dual still refills (slope even
   slightly steeper late); Em_end ≈ law dual. Trap resists O but dual I is
   trap-insensitive at this Em.
4. **Quiet dual = real I channel on the clump** without bath flood — confirmed
   at both detunes to T=5k. **Not an Em\* monodisperse attractor** as implemented.
5. **Adoption still blocked** on park demonstration; leave `rad_dual` default 0.

## Class hygiene

| class | this suite |
|---|---|
| compact HOLD (no rad) | not tested here (see R1) |
| SLOW_DRAIN (O tax) | d0 both traps |
| SLOW_REFILL (I>O) | d1 both traps |
| Em\* park | **none** |
| box-fill | **none** (quiet apparatus correct) |

*Also written under `runs/p2quiet_T5000_20260809_010524/P2QUIET_T5K.md`.*
