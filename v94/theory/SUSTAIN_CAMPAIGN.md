# Sustain campaign T=3000 — plan, Bohmian link, isolation

**Status:** RUNNING **25/29** (2026-08-08 ~16:52 CDT) — final batch G05/R02/S00b/S11b  
**Live read:** `MECHANISM.md` + draft `SUSTAIN_REPORT.md` — compact HOLD yes; genesis diffuse; R1 ready
**Kernel:** Go `fabrun` (multi-process parallel; GOMAXPROCS per worker)  
**Horizon:** T=3000 (unless noted)  
**Harness:** `theory/bin/probe -suite sustain -kernel fabrun -T 3000 -j 4`  
**Live read:** `MECHANISM.md` + run-dir `INTERIM_SCORE.md` (S00/S04 HOLD; S05 tax; S07 box-fill)

---

## 1. Question

Can a condensate **sustain** (Em stays O(1) with small late-time slope), and
**which switches** enable that? If yes, map the mechanism (trap vs door
intake vs tax vs progressive phase) so it can later become *internal*
dynamics rather than external fuel.

Genesis (C1/C2) already showed: shear barely changes birth; workfn gates
pure-field birth. Sustain is a **different bar**.

---

## 2. Isolation matrix (what each arm turns)

| arm | uni | gate | qd | rad | wf | noise | other | isolates |
|---|---|---|---|---|---|---|---|---|
| S00 | ✓ lin | off | 1.2 | 0 | 1 | 0 | | **baseline hold** |
| S01 | add | — | 1.2 | 0 | 1 | 0 | | transport algebra |
| S02 | ✓ gat | on | 1.2 | 0 | 1 | 0 | | phase-gated hops |
| S03 | ✓ | | **0** | 0 | 1 | 0 | | binding off (melt) |
| S04 | ✓ | | **3.6** | 0 | 1 | 0 | | deep trap |
| S05 | ✓ | | 1.2 | **0.05** | 1 | 0 | | radiance outtake |
| S06 | ✓ | | 3.6 | **0.05** | 1 | 0 | | deep vs tax |
| S07 | ✓ | | 1.2 | 0 | 1 | **1.5** | | field fuel (I) |
| S08 | ✓ | | 1.2 | **0.05** | 1 | **1.5** | | I and O together |
| S09 | ✓ | | 3.6 | 0 | 1 | **1.5** | | deep + fuel |
| S10 | ✓ | | 1.2 | 0 | **0** | 0 | | open door, no fuel |
| S11 | ✓ | | 1.2 | 0 | **0** | **1.5** | | open door + fuel |
| S12 | ✓ | | 1.2 | **0.05** | **0** | **1.5** | | **metabolic candidate** |
| S13 | ✓ | | 3.6 | 0 | **0** | **1.5** | | deep metabolic |
| S14 | ✓ | | 1.2 | 0 | 1 | 0 | **kx** | progressive dense phase |
| S15 | ✓ | | 1.2 | 0 | 1 | **1.5** | **kx** | progressive + fuel |
| S16 | ✓ | | 1.2 | 0 | 1 | 0 | **seed_mw** | winding progression |
| S17 | add | | 1.2 | 0 | 1 | **1.5** | | add + fuel |
| S18 | ✓ | | 1.2 | **0.05** | 1 | 0 | | V3a quiet |
| S19 | ✓ | | 1.2 | **0.05** | 1 | **1.5** | | V3a + fuel |
| G01/G02 | — | | | 0 | 0 | pulse kx/k0 | genesis long shear |
| G03 | | | | 0 | 1 | pulse kx | workfn blocks birth? |
| G04/G05 | | | | **0.05** | 0 | pulse | genesis + tax |
| R01/R02 | ✓ | | | 0 | 0 | ring mw | topology long |

**Sustain pass (operational):** last-half `|dEm/dt| ≲ 0.05` and Em_end ≳ 20
→ classify **HOLD** (see `classifySustain` in probe).

---

## 3. What is “Bohmian” here (and what is not)

**Not:** pilot particles, hidden variables, config-space N-body ψ, tracers.

**Yes (methodology from `BOHMIAN_HYBRID` / `PROGRESSIVE_CELL`):**

| idea | where in campaign | what we already saw (short T) |
|---|---|---|
| Transport = guidance of amplitude | `fd` / pass U (S00,S14,S15) | B1 fd≈+284; additive uses `fl` not `fd` |
| Harmonics are line-like (J on links) | progressive meters Jx, Jd | field Jx kill pairs A6 vs A7 |
| Object ≠ cell tags | Em series, not Em_tag alone | F_tag falls while Em holds (Face C) |
| Phase progression on dense ψ | S14 kx, S16 seed_mw | ring holds W in vacuum (Face A) |
| Global R₂ is not quality | meters report R2f but do not rank by it | traveling wind R2f≈0, standing R2f high |

**What Bohmian framing predicts for sustain:**

- A **held** object is a support \(S(t)\) with stable progressive structure
  (trap and/or continuous guidance), not a fixed index set.
- External field fuel (noise) is an *open-system* \(I\); internal sustain
  would need the same books closed without noise (self-recapture) — that
  is the translation target if S11/S12 HOLD but S00 does not.

---

## 4. How to run / monitor

```bash
cd /home/d/code/scp/v93
go build -o theory/bin/probe ./theory/probe
# ~28 arms × ~1 h / 4 workers ≈ 7 h
setsid env GOMAXPROCS=4 ./theory/bin/probe \
  -suite sustain -kernel fabrun -T 3000 -j 4 -seed 20260802 \
  -snap-every 5000 \
  > theory/runs/sustain_launch.log 2>&1 &
```

Outputs per arm: `.log`, `.fcs`, `.em.csv`; campaign: `results.csv`,
`SUSTAIN_REPORT.md`, `results.json`.

---

## 5. After completion — mechanism reading

| pattern | interpretation |
|---|---|
| S00 HOLD, S03 MELT | q_detune trap necessary |
| S00 MELT, S04 HOLD | need deep trap |
| S00 MELT, S07/S11 HOLD | need external field fuel |
| S11 HOLD, S12 HOLD with rad | metabolic balance possible |
| S00 HOLD, S05 DRAIN | radiance kills quiet matter |
| S14/S16 differ from S00 | progressive phase matters for hold |
| G01 plateau vs G02 | shear helps long genesis/growth |

If only fuel-fed arms HOLD → sustain is **open-system**, not yet internal.
Next design: coherent self-recapture (Face B open arm) to internalize \(I\).
