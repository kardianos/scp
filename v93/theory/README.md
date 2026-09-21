# v93/theory — analytical notes + Go numeric probes

## ★ STATUS COMPLETE (start here)

**→ [`FINAL_STATUS_ADOPTION.md`](FINAL_STATUS_ADOPTION.md)** — **adoption judgments**
(§1 high-level + §2 item-by-item theory/code/numerics for every lever).  
**→ [`STATUS_COMPLETE.md`](STATUS_COMPLETE.md)** — FCS archive + WebP/PNG paths
(`/space/scp/v93/`), scoreboard index.

Hourly loop **cancelled** 2026-08-09. Prefer **Go `fabrun`**. Kernel edits
only when measured and battery-safe at defaults.

## Campaign status (sustain 08-08 → 08-09) — CLOSED

| deliverable | path | status |
|---|---|---|
| **Status complete** | **`STATUS_COMPLETE.md`** | **DONE** — entry point |
| T=3000 matrix | `SUSTAIN_REPORT.md` · `runs/ACTIVE_SUSTAIN` | **DONE** 29/29 |
| R1 T=10000 | `R1_SCORE.md` · `runs/ACTIVE_R1` | **DONE** 6/8 HOLD, 2 rad SLOW_DRAIN |
| Mono / dual | `MONO_FINAL.md` · `P2QUIET_T5K.md` | **DONE** — no Em\* park; dual SLOW_REFILL |
| Judgment + next | **`CAMPAIGN_JUDGMENT.md`** | **DONE** — sustain YES (uni+qd internal) |
| Showcase viz | `/space/scp/v93/{webp,fcs,png}/` | **DONE** |
| Scoreboard trail | `NUMERICS.md` · `MECHANISM.md` · `HOURLOG.md` | updated |

**Headline:** quiet compact HOLD is real through T=10k under **unitary dense hop
+ law q_detune**. Radiance taxes; dual intakes but does not park; noise box-fill
≠ particle; shear ≠ birth switch; battery **87 GREEN** (dual default-off).

## Documents

| note | topic |
|---|---|
| **`STATUS_COMPLETE.md`** | **Campaign closeout + FCS + viz pointers** |
| **`CAMPAIGN_JUDGMENT.md`** | Mission answer + next-internal recommendation |
| **`SUSTAIN_REPORT.md`** | T=3000 structural scoreboard |
| **`R1_SCORE.md`** | T=10000 endurance |
| **`NUMERICS.md`** | Full campaign scoreboard |
| `MECHANISM.md` | Isolation → mechanism |
| `METABOLIC_MATTER.md` | I/O balance, T⁴ radiance, melt, DM-as-field-density |
| `T_HARMONICS_AND_WIND.md` | T as harmonic attribute; field shear ≠ T |
| `SHEAR_GENESIS_DESIGN.md` | Shear genesis design (pre-reg); see NUMERICS |
| `BOHMIAN_HYBRID.md` | Bohmian as **methodology** — morph of the cell, not tracers |
| `PROGRESSIVE_CELL.md` | Symbolic progressive-cell + early probes |

## Go probe suite

```
theory/probe/     # parallel freecell driver + progressive meters
theory/Makefile   # make probe / probe-run / probe-shear / probe-prog
theory/bin/probe  # binary
theory/runs/      # probe_<timestamp>/{logs,fcs,SUMMARY.md,results.csv}
```

```bash
cd /home/d/code/scp/v93
make -C theory probe-run J=6
```

**22 arms in parallel:** shear winds, spin, imp, high-amp kill pairs, noise,
slit, vacuum (invalid), unitary/additive/gated dense, melt/deep, ring winding.

**Latest verdict:**

| run | result |
|---|---|
| `probe_20260808_102445` (22 arms, wf_on=1) | pure-field genesis **cond≡0**; field \(J_x\) **positive**; dense fd **+284** |
| `probe_20260808_102702` (door-side C*) | **wf_on=0 unlocks birth** (cond~480); shear vs standing birth **nearly equal**; blob+noise births under V3a workfn |

Details: **`NUMERICS.md`**.
