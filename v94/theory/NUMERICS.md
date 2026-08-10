# theory/ numerics — campaign status (Go probe)

**Last parallel runs:**

| id | dir | focus |
|---|---|---|
| full matrix | `theory/runs/probe_20260808_102445` | 22 arms shear+prog |
| door-side | `theory/runs/probe_20260808_102702` | C1–C5 workfn/presence |

**Tool:** `theory/bin/probe` (pure Go; drives `v93/freecell`; meters via `fcsdump`+Go)

---

## 0. Honest headline

| claim | status |
|---|---|
| **Shear as unique genesis driver** (organized field beats equal-energy standing for birth) | **WEAK / mostly negative** — see C1 vs C2 |
| **Genesis under default workfn (wf_on=1) from pure field** | **NEGATIVE** (cond≡0 on all A/S arms) |
| **Genesis when door is not presence-gated (wf_on=0)** | **STRONG** — cond~470–490, Em~450+ |
| **Genesis needs matter presence or open door** | **CONFIRMED** — C4 quiet blob cond=0; C3 blob+noise cond=498 |
| **Progressive field structure** (signed \(J_x\)) | **POSITIVE** — A6 vs A7 kill pair |
| **Progressive dense guidance** (unitary `fd`) | **POSITIVE** — B1 ≈+284; C5 ≈+440 with noise |
| **Global \(R_2\) as quality** | **AVOID** — inverted for traveling waves |

**Root of earlier “all negative” genesis:** laws_V3a **workfn** puts condensation
demand on *bound matter* (empty cells → wf_far). Pure field shear never
opens the door. That is a **door/regime** result, not a progressive-meter failure.

### T=3000 sustain (partial — 4/29 done, 2026-08-08 11:55 CDT)

| arm | T=3000 Em | structure | verdict |
|---|---:|---|---|
| **S00 uni+qd1.2 quiet** | 148.3 flat | nHot~50, Em_max~1.9 | **STRUCTURAL HOLD** |
| S01 additive | 77 drain | nHot~4 | DRAIN |
| S02 gated | 148.8 flat | nHot→0 | STRUCT_MELT |
| S03 qd=0 | 148.8 flat | diffuse | DIFFUSE |

→ Quiet unitary + law detune can sustain without external fuel. See `MECHANISM.md`.

---

## 1. What changed in theory/ (this face)

| path | role |
|---|---|
| `probe/*.go` | parallel freecell driver + log parse + progressive meters |
| `Makefile` | `make -C theory probe-run` |
| `bin/probe` | built binary |
| `runs/probe_*/` | logs, fcs, `SUMMARY.md`, `results.csv`, `results.json` |
| docs (this file, README) | updated; older design notes retained |

No Python in this suite. Kernel untouched.

---

## 2. Shear / field matrix (genesis hunt)

All arms: `cond = 0`. Selected progressive field meters (last frame):

| arm | Ee | Jx | \|J\| | R2f | qcnvF | reading |
|---|---:|---:|---:|---:|---:|---|
| S0_pulse_k0 | 548 | **~0** | 0.08 | **0.76** | 0 | standing sterile |
| A1_kx06 | 548 | −0.16 | 0.26 | 0.12 | 0 | mild wind |
| A2_kx12 | 548 | −0.28 | 0.45 | 0.002 | 0.2 | wind (prior) |
| A3_kx20 | 548 | −0.16 | 0.43 | 0.01 | 7 | steep — Jx not monotone in k |
| A4_spin2 | 548 | ~0 | 0.22 | 0.002 | 0 | azimuthal; weak net Jx |
| A5_imp_kx | 548 | **−0.54** | 0.93 | 0.08 | 72 | wind+radial strong \|J\| |
| A6_hiamp_kx | 1645 | **−1.04** | 1.73 | 0.001 | 230 | strongest progressive |
| A7_hiamp_k0 | 1645 | **~0** | 0.34 | **0.73** | 68 | energy twin of A6 — standing |
| A8 noise+kx | 1905 | −0.12 | 1.13 | 0.006 | 114 | wind diluted by noise |
| A9 noise | 1372 | ~0 | 0.80 | 0.02 | 22 | isotropic |
| A10 long kx | 1249 | −0.42 | 0.70 | 0.005 | 28 | longer T, still cond=0 |
| A12 slit | 62 | −0.27 | 0.80 | 0.016 | 0.5 | directed beam structure |
| A13 spin+kx | 548 | −0.29 | 0.49 | 0.08 | 0 | combined progressive |
| S0b noise2 | 2745 | −0.03 | 1.51 | 0.01 | 370 | energy≠structure |

**Kill-pair confirmation (equal Ee):**

- A6 vs A7 at Ee=1645: Jx −1.04 vs ~0, R2f 0.001 vs 0.73  
- S0 vs A2 at Ee=548: same pattern  

**Avoid:** ranking by R2f or by qcnvF alone (noise/high amp raise credits without genesis).

**Apparatus failures / weak arms:**

- A11_vac_kx: `bath=0` pulse → z_live=0 (no substrate) — **invalid**, do not reuse  
- B6_ring_mw: runs but tiny system; use for winding, not genesis  

---

## 3. Progressive dense / Bohmian-method matrix

| arm | fd | fl | Em_tot | Em_max | PR | F_tag | Jd |
|---|---:|---:|---:|---:|---:|---:|---:|
| B1 unitary linear | **+284** | 0 | 93 | 2.79 | 146 | 0.70 | 0.30 |
| B2 additive | 0 | **−109** | 78 | 0.76 | 336 | 0.83 | 0.23 |
| B3 unitary gated | +5.5 | 0 | 93 | 0.74 | 413 | 0.89 | 0.25 |
| B4 uni melt qd=0 | 0 | 0 | 93 | 0.86 | 390 | 0.73 | 0.16 |
| B5 uni deep qd=3.6 | 0 | 0 | 186 | 2.39 | 306 | 0.86 | 0.55 |
| B6 ring seed_mw | 0 | 0 | 4.9 | 0.81 | 6 | 1.0 | **1.42** |

**Works:** linear unitary hop current as guidance meter; deep trap raises Em_max/Jd; ring holds high dense |J| on cycle.  
**Avoid:** gated τ for transport (fd collapses); comparing Em_tag across additive vs unitary; treating melt residual R2 as success (prior Face C lesson).

---

## 3b. Door-side matrix (probe_20260808_102702)

Hypothesis: workfn presence-gate blocks pure-field birth. Test wf_on=0 and
matter+field co-presence.

| arm | setup | cond | Em_tot | Jx | R2f | fd |
|---|---|---:|---:|---:|---:|---:|
| **C1_nowf_kx** | pulse amp=3 kx=1.2 **wf_on=0** | **488** | 474 | **−0.16** | 0.03 | 0 |
| **C2_nowf_k0** | pulse amp=3 kx=0 **wf_on=0** | **472** | 452 | ~0 | **0.36** | 0 |
| **C3_blob_field** | blob + noise_amp=2 wf_on=1 | **498** | 586 | −0.02 | 0.02 | 0 |
| **C4_blob_quiet** | blob no field | **0** | 84 | 0 | 0 | 0 |
| **C5_uni_noise** | unitary blob + noise | **220** | 325 | ~0 | 0.02 | **+440** |

### Door-side reading

1. **wf_on=0 unlocks massive condensation** from field alone (C1/C2). The
   V3a workfn was the genesis wall for pure-field arms in the full matrix.
2. **Shear is not the main birth switch:** C1 vs C2 at equal amp — cond
   488 vs 472 (~3%), Em 474 vs 452 (~5%). Progressive Jx still distinguishes
   them (structure yes; genesis boost marginal).
3. **Matter presence + field** (C3) births under wf_on=1; quiet blob (C4)
   does not — door needs field energy when presence-gated.
4. **C5** combines progressive transport (fd=+440) with genesis (cond=220)
   — most “complete” single arm for metabolic follow-ups.

**Shear-genesis thesis (strong form) is not supported:** organization does
not dominate birth once the door is open. **Weak form** (slight assist)
remains possible; needs tighter energy-matched ensembles.

---

## 4. Parallel approaches — keep / drop / next

| approach | keep? | reason |
|---|---|---|
| Phase-wind pulses (vary kx, amp) | **keep** as structure probe | Jx bar solid |
| Standing equal-Ee control | **keep** | mandatory kill pair |
| High-amp wind | **keep** | max \|Jx\| |
| Noise / noise+wind | **keep** | energy≠structure |
| Vacuum bath=0 pulse | **drop** | empty geometry |
| Unitary linear e3b | **keep** Bohmian-method core | fd real |
| Unitary gated | demote for transport | fd~0 |
| **wf_on=0 field birth** | **keep** as genesis apparatus | door wall identified |
| **blob + field (C3/C5)** | **keep** | genesis under V3a workfn |
| Shear as sole birth mechanism | **demote** | C1≈C2 birth |
| Coherent phase-matched driver | **next** | Face B open; may beat random noise for *ordered* refill |
| Energy-matched scramble of same pulse | **next** | cleaner than noise S0 for structure bars |
| Metabolic I=O on C5-like objects | **next** | after genesis, stability |
---

## 5. How to re-run (Go only)

```bash
cd /home/d/code/scp/v93
make -C theory probe          # build
make -C theory probe-run J=6  # full matrix
make -C theory probe-shear    # field arms only
make -C theory probe-prog     # dense progressive only
# or:
./theory/bin/probe -j 8 -only shear,ctrl
./theory/bin/probe -name hiamp
```

Outputs land in `theory/runs/probe_<timestamp>/`.

---

## 6. Relation to earlier theory notes

| note | update |
|---|---|
| `SHEAR_GENESIS_DESIGN.md` | design stands; Phase-0 result = genesis **negative**, structure meters **positive** |
| `PROGRESSIVE_CELL.md` | confirmed by Go remeters (no Python) |
| `BOHMIAN_HYBRID.md` §10 | method morph validated via fd / Jx — not tracers |
| `METABOLIC_MATTER.md` | I=O still unblocked; shear is not T and not yet I |

---

*Genesis wall: organized field ≠ condensation under present door.
Progressive wall: mostly cleared for meters; next physics is door-side.*

---

## 7. Long-horizon sustain T=3000 — COMPLETE

**Dir:** `theory/runs/sustain_T3000_20260808_104548` — **29/29 DONE** (`.done` 17:23 CDT).  
**Final structural report:** run-dir `SUSTAIN_REPORT.md` · `MECHANISM.md` (Em-only auto-HOLD retired).

| claim | status |
|---|---|
| Quiet unitary compact HOLD | **YES** S00/S00b/S04/S10/S16 |
| Progressive kx | **HOLD_ext** S14 fd=+650 |
| Ring topology | **HOLD_tiny** R01 Em=4.88 nHot=6; R02 bath dephases nHot→1 |
| Additive/gated/qd=0 | structure fail |
| Radiance quiet | SLOW_DRAIN S05≡S18; deep S06 partial resist |
| Noise fuel | box-fill Em~2200–2370 nHot~1300+ — not particle |
| Pure-field genesis long | G01/G02 Em~640 **diffuse**; G03 empty; G04/G05 taxed ~−33% |
| Shear as genesis switch | **weak** G01≈G02 / G04≈G05 (~2%) |
| Seed robustness | S00b HOLD; S11b box-fill |

### Compact HOLD scoreboard (last frame)

| arm | Em | max | nHot | class |
|---|---:|---:|---:|---|
| S00 | 148.3 | 1.93 | 50 | HOLD |
| S00b s111 | 146.9 | 1.82 | 56 | HOLD |
| S04 deep | 147.4 | 2.51 | 51 | HOLD |
| S14 kx | 147.7 | 1.84 | 56 | HOLD_ext |
| S16 mw | 148.1 | 1.90 | 51 | HOLD |

**R1 endurance FINAL:** `R1_T10000_20260808_172409` — **8/8 done** (see `R1_SCORE.md`):

| arm | Em_end | max | nHot | verdict |
|---|---:|---:|---:|---|
| R1_S00 | 148.30 | 1.71 | 52 | **PASS HOLD** |
| R1_S04 | 147.29 | 1.99 | 55 | **PASS HOLD** |
| R1_S14 | 147.71 | 1.56 | 50 | **PASS HOLD_ext** |
| R1_S16 | 148.07 | 1.66 | 50 | **PASS HOLD** |
| R1_S00_s111 | 146.85 | 1.68 | 53 | **PASS HOLD** |
| R1_S00_amp12 | 206.75 | 1.80 | 106 | **PASS HOLD** |
| R1_S00_rad | 105.41 | 0.75 | 18 | **SLOW_DRAIN** no Em* |
| R1_S04_rad | 113.90 | 0.66 | 22 | **SLOW_DRAIN** deep +8 |

**P2 quiet park T=2000 FINAL** (`p2quiet_T2000_*`): qd12 d0=113 drain / d1=134 **recover** (ratio 1.18); qd36 d0=122 / d1=129 (ratio 1.05). No box-fill; not Em\* park.

**Quiet dual T=5000 FULL FINAL** (`p2quiet_T5000_20260809_010524`, `.done` 04:02):

| arm | Em_end | slope 3k→5k | class |
|---|---:|---:|---|
| qd12 d0 | **108.1** | −0.0009 | SLOW_DRAIN |
| qd12 d1 | **140.4** | **+0.0023** | **SLOW_REFILL** (ratio 1.30) |
| qd36 d0 | **117.2** | −0.0012 | SLOW_DRAIN (deep +9 vs law) |
| qd36 d1 | **139.1** | **+0.0030** | **SLOW_REFILL** (ratio 1.19) |

→ Dual = real quiet I; **no Em\* park** at either detune; deep trap does not stop refill. See `P2QUIET_T5K.md`.

**Mono scan FINAL:** `mono_T2000_*` — **8/8**. k_rad>0 → SLOW_DRAIN only. See `MONO_FINAL.md`.

**P2 bath suite FINAL** (`p2_T600_*`): dual raises Em (1.08→2.74) but **box-fill**. **Not adopted.** See `P2_RESULTS.md`.

**Battery p2fix2:** **ALL GREEN 87/87** (exit 0). abx PASS after Em_tot mirror + banner/credit mask. Dual default-off.

**qgraph P3/P5:** p_rad=2 drains 3× vs p=4; door0≈door1 box-fill. See `P3_P5_RESULTS.md`.

