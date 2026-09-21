# Campaign judgment — sustain / mechanism / next dynamics

**Deadline target:** 2026-08-09 11:00 CDT  
**As of:** 2026-08-09 ~04:05 CDT  
**Evidence:** T=3000 matrix · R1 T=10000 · mono · P2 dual T=5k full · battery GREEN

---

## 1. Does a condensate sustain?

**Yes.** Under **unitary linear dense hop (`amp_nat`) + law `q_detune` trap**, a compact multi-hot clump holds Em flat with nHot≳50 and Em_max≳1.5 from T=3000 through **T=10000**.

| evidence | result |
|---|---|
| T=3000 S00/S04/S14/S16 (+S00b,S10) | structural HOLD / HOLD_ext |
| R1 T=10000 same + s111 + amp12 | **PASS HOLD** (6/8) |
| R1 rad arms | SLOW_DRAIN — tax, not park |

**Required switches:** unitary transport + q_detune>0.  
**Fails:** additive (DRAIN+MELT), gated τ (STRUCT_MELT), qd=0 (DIFFUSE).

---

## 2. Why (mechanism isolation)

| switch | verdict |
|---|---|
| Transport algebra | Unitary holds; additive drains |
| Trap (q_detune) | Essential for blob class |
| Progressive kx / seed_mw | Compatible (HOLD_ext / winding OK) |
| Workfn door | Quiet blob insensitive; pure field absolute wall |
| Radiance O | SLOW_DRAIN; deep resists; **no Em\*** |
| Noise fuel | box-fill Em~2200+ — **not a particle** |
| Genesis open-door | diffuse Em~640 plateau; shear **weak** (~2%) |
| Shear as birth switch | **Rejected** (G01≈G02) |

**Internal vs external:** quiet HOLD is **internal** (trap+U). Noise/bath is external capacity. Keep classes separate.

---

## 3. Monodispersity / Em\* park

- k_rad ladder (mono T=2000): monotone tax, no fixed point.  
- R1 rad T=10k: still SLOW_DRAIN; peaks weaken (Em_max&lt;1).  
- **P2 dual intake** (experimental, default 0) — **T=5000 quiet FULL FINAL:**  
  | | d0 Em | d1 Em | class d1 |
  |---|---:|---:|---|
  | qd=1.2 | 108 | 140 | SLOW_REFILL (ratio 1.30) |
  | qd=3.6 | 117 | 139 | SLOW_REFILL (ratio 1.19) |
  - With bath_xb → box-fill (wrong apparatus).  
  - Quiet dual: real I on clump; deep trap resists O (+9) but **does not park**.  
- Dual **not adopted** (no Em\* fixed point at either detune through T=5k).

---

## 4. Kernel / battery

- Go+C: `rad_dual`, `bath_xb`/`bath_em`; Em_tot on both diag lines.  
- **Battery ALL GREEN 87/87** after Em_tot mirror + abx banner/credit mask (`battery_p2fix2_*`, exit 0).  
- Physics identity C≡Go holds on masked abx; dual remains off by default.

---

## 5. Next-internal-dynamics recommendation

1. **Sustain channel of record:** unitary dense hop + law q_detune — validated to T=10k; seed- and mass-robust.  
2. **Adoption path for amp_nat:** battery already green at defaults; optional R0 amp_nat≠0 panel later — do not flip law table without user sign-off.  
3. **Radiance:** keep as far-field tax; do not expect monodisperse Em\* from O alone.  
4. **Dual (P2):** T=5k closed — **SLOW_REFILL not park** at law and deep trap; leave default 0. Next design = balance O/I (not longer T alone). Adopt only after Em\* + GREEN battery, dual still default 0.  
5. **Do not** path atoms through noise box-fill.  
6. **Genesis:** separate coarsening problem (diffuse → compact?).  
7. **Bohmian-method:** progressive Jx/fd are real diagnostics; not pilot particles.

---

## 6. Deliverable checklist (CALLING)

| # | item | status |
|---|---|---|
| 1 | T=3000 matrix scored | **DONE** `SUSTAIN_REPORT.md` |
| 2 | R1 T=10000 HOLD winners | **DONE** `R1_SCORE.md` |
| 3 | Mono / rad / progressive probes | **DONE** mono + R1 rad + quiet dual T=5k full |
| 4 | MECHANISM + NUMERICS + HOURLOG | **DONE** (trail live) |
| 5 | Next-internal recommendation | **THIS FILE §5** |

*Quiet dual T=5000 FULL FINAL 04:02 CDT — SLOW_REFILL both traps; no Em\* park. Sustain yes/no already closed by R1 HOLD. Freeze writeup trail by 11:00; no new multi-hour arms after 09:00.*
