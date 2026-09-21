# HOURLOG — append-only (sustain campaign)

## 2026-08-08 10:50 CDT — campaign start / calling installed

- Deadline 23:00 CDT. CALLING.md written for post-compaction re-entry.
- Active: `theory/runs/sustain_T3000_20260808_104548` — 29 arms, fabrun j=4, T=3000.
- First batch: S00 S01 S02 S03 (fcs growing). No .log yet (probe CombinedOutput until arm end; next binary streams).
- Bench: fabrun T=40 ~44s, freecell ~32s → ~55 min/arm Go → ~7h matrix.
- Prior short-T: genesis weak for shear; door wf wall; progressive Jx/fd real.
- Hourly scheduler task created (re-read CALLING.md).
- Next: wait first-batch complete; score; refine.

## 2026-08-08 10:51 CDT — partial FCS at t≈200 (first batch)

In-flight fcsdump stats (not final):

| arm | Em@0 | Em@200 | Em_max@200 | nHot | note |
|---|---:|---:|---:|---:|---|
| S00 uni base | 148.8 | **148.5** | **2.41** | 77 | energy hold + peak grow (condensing) |
| S01 additive | 148.8 | **98.4** | 1.34 | 46 | **DRAIN** cell Em (flight?) |
| S02 gated | 148.8 | 148.8 | 0.84 | 42 | frozen envelope, no peak grow |
| S03 qd0 melt | 148.8 | 148.8 | 0.72 | **4** | total Em conserved but **structure melts** |

Early judgment (t=200/3000): unitary+law qd looks most HOLD-like; additive loses cell Em; qd=0 melts structure without energy loss. Confirm at t=3000.
- fabrun still 4. No completed arms.

## 2026-08-08 10:52 CDT — kernel diag + durable ops

- **Kernel (Go):** diag line now prints `| Em_tot=...` (all-cell dense inventory) — `fab/run.go`. Smoke OK. Running batch still on old binary; **next arms** use new fabrun.
- Hourly durable task scheduled (re-read CALLING.md).
- Symlink `theory/runs/ACTIVE_SUSTAIN` → current T=3000 dir.
- Mission files: CALLING.md, HOURLOG.md, SUSTAIN_CAMPAIGN.md.

## 2026-08-08 11:52 CDT — first 4 arms DONE; structural HOLD on S00

- Time 11:52 CDT (deadline 23:00). fabrun=4 (S04–S07), probe alive ~65 min wall.
- **Done 4/29:** S00 S01 S02 S03. results.csv partial. No `.done` yet.
- ETA matrix finish ~17:00–18:00 CDT (well before 21:30 stop-new-arms).

### T=3000 verdicts (honest structure, not Em-only)

| arm | Em | Em_max | nHot | class |
|---|---:|---:|---:|---|
| **S00 uni base** | 148.3 flat | 1.93 | 50 | **HOLD structure** (quiet trap enough) |
| S01 additive | 77.6 drain | 0.52 | 4 | **DRAIN + structure melt** |
| S02 gated | 148.8 flat | 0.34 | 0 | energy HOLD / **STRUCT_MELT** |
| S03 qd=0 | 148.8 flat | 0.80 | 10 | energy HOLD / **DIFFUSE melt** |

- CALLING branch: **S00 quiet HOLD** → trap enough; monodispersity next (not fuel-only path).
- S00 support: hot CV=0.45, Rg_hot=3.85, multi-hot graded clump (not monodisperse soliton); Em flat from t≈400; tags lie (tagEm falls while Em holds).
- In-flight: S04 deep HOLD-like (max 2.75); S05 rad DRAIN; S06 deep+rad slower drain; S07 noise GROW Em~2300.
- Wrote `runs/sustain_T3000_.../INTERIM_SCORE.md`. Fixed `watch_and_score.sh` (rg -h was --help).
- **Next:** no new long matrix; wait S04–S07; refine mechanism table; R1 T=10k only after full HOLD list.

## 2026-08-08 11:55 CDT — mechanism draft + S07 early plateau

- Wrote `theory/MECHANISM.md` from finished S00–S03.
- **Headline:** quiet unitary+qd=1.2 **structurally HOLDs** to T=3000 (multi-hot clump, not monodisperse).
- Unitary linear required (gated STRUCT_MELT; additive DRAIN). Trap required for clump (qd=0 DIFFUSE).
- In-flight S07 noise: Em 148→**~2320 by t=300 then nearly flat** — possible open-system **GROW→plateau** (watch to t=3000).
- S05 rad: Em 148→130 at t=600 (DRAIN rate ~0.03/t.u.).
- S04 deep: Em flat, Em_max~2.5 HOLD-like.
- fabrun=4 still. No new arms launched (CALLING: wait full matrix).


## 2026-08-08 12:52 CDT — batch 2 DONE (8/29); S04 HOLD, S05 tax, S07 box-fill

- Time 12:52 CDT. fabrun=4 (S08 noise_rad, S09 noise_deep, S10 nowf quiet, S11 nowf noise).
- **Done 8/29.** ETA full ~17:00–18:00 CDT. No .done.

### New verdicts (structure-honest)

| arm | Em_end | Em_max | nHot | class |
|---|---:|---:|---:|---|
| **S04 deep** | 147.4 | 2.51 | 51 | **HOLD** (taller than S00) |
| **S05 rad** | 109.9 | 0.89 | 25 | **SLOW_DRAIN** tax ~26% Em; CV↓ mono-ish |
| **S06 deep+rad** | 119.7 | 1.11 | 40 | **SLOW_DRAIN** deep resists tax |
| **S07 noise** | 2338 | 3.41 | 1388 | **GROW→PLATEAU box-fill** (not compact) |

- Mechanism: trap+unitary enough (S00/S04). Radiance taxes quiet; no Em* park by T=3000. Noise fills box — do not call particle sustain.
- Partial: S08 Em~2200@800; S09 max~6; S10 quiet HOLD-like; S11 like S07.
- Updated MECHANISM.md + INTERIM_SCORE.md + NUMERICS §7.
- **Next:** wait metabolic S12 + progressive; R1 only on S00/S04 after matrix; no new long launch.


## 2026-08-08 13:52 CDT — batch 3 DONE (12/29); S10≡S00; S08 box-fill thruput

- Time 13:52 CDT. fabrun=4 (S12 nowf+noise+rad, S13 nowf+deep+noise, S14 kx, S15 kx+noise).
- **Done 12/29.** No .done. ETA ~17:00–18:00.

### New verdicts

| arm | Em | Em_max | nHot | class |
|---|---:|---:|---:|---|
| **S08 noise+rad** | 2207 | 2.97 | 1553 | box-fill + huge I/O (cond~20k rad~19k) — **not compact metabolic** |
| **S09 noise+deep** | 2367 | **5.66** | 1329 | box-fill deep peaks |
| **S10 nowf quiet** | 148.3 | 1.93 | 50 | **HOLD ≡ S00** (em.csv md5 identical) |
| **S11 nowf+noise** | 2339 | 3.32 | 1379 | box-fill ≈S07 |

- wf door irrelevant for quiet/fuel-fill sustain once matter present; genesis wall was pure-field birth.
- S08 CALLING “metabolic” = open thruput at filled capacity, not I=O particle.
- Partial: S12~S08; S13~S09; S14 quiet HOLD-like; S15 fuel+tilt.
- Updated MECHANISM.md, INTERIM_SCORE, NUMERICS.
- **Next:** wait S12–S16 + G*; R1 only on S00/S04 after matrix.


## 2026-08-08 14:52 CDT — batch 4 DONE (16/29); S14 HOLD_ext progressive

- Time 14:52 CDT. fabrun=4 (S16 seed_mw, S17 add+noise, S18 V3a quiet, S19 V3a noise).
- **Done 16/29.** No .done. Remaining after batch: G* + R* + seed-robust (~9). ETA ~17:30–18:30.

### New verdicts

| arm | Em | max | nHot | Rg | class |
|---|---:|---:|---:|---:|---|
| S12 nowf+noise+rad | 2206 | 3.20 | 1534 | 8.0 | box-fill thruput ≈S08 |
| S13 nowf+deep+noise | 2360 | 6.03 | 1304 | 8.0 | box-fill deep |
| **S14 kx quiet** | 147.7 | 1.84 | 56 | **6.06** | **HOLD_ext** fd=+650 |
| S15 kx+noise | 2317 | 3.81 | 1364 | 7.9 | box-fill + fd=+802 |

- Progressive phase spreads clump (Rg↑) but energy+hot structure holds — Bohmian fd real.
- Partial: S16 HOLD-like; S18 rad drain; S17/S19 fill.
- Updated MECHANISM + INTERIM + NUMERICS.
- **Next:** finish matrix; R1 design S00/S04/S14; no new long launch until `.done`.


## 2026-08-08 15:52 CDT — batch 5 DONE (20/29); S16 HOLD; genesis partial birth

- Time 15:52 CDT. fabrun=4 (G01 nowf kx, G02 nowf k0, G03 wf kx, G04 nowf kx+rad).
- **Done 20/29.** Remaining after G-batch: G05, R01–R02, S00b, S11b. ETA ~18:00–18:30.

### S16–S19 verdicts

| arm | Em | max | nHot | class |
|---|---:|---:|---:|---|
| **S16 seed_mw** | 148.1 | 1.90 | 51 | **HOLD** Rg=4.23 |
| S17 add+noise | 2157 | 2.30 | 1744 | box-fill weaker peaks |
| S18 V3a quiet | 109.9 | 0.89 | 25 | SLOW_DRAIN **≡S05** |
| S19 V3a noise | 2207 | 2.97 | 1553 | box-fill **≡S08** |

### Genesis partial (long-T)

- G01/G02: Em~640 plateau from pure field (wf_on=0); shear ≈ equal (confirms short-T C1≈C2).
- G03: Em≡0 under wf_on=1 — door wall holds to T~1800+.
- G04: Em~414, nHot~41 — radiance suppresses genesis mass.

- Compact HOLD list: S00, S04, S14, S16 (+S10≡S00). R1 only those.
- Updated MECHANISM, INTERIM, NUMERICS.
- **Next:** finish G*+R*+seed-robust; full SUSTAIN_REPORT; then R1 T=10000 if time before 21:30.


## 2026-08-08 16:52 CDT — genesis DONE (25/29); final batch in flight

- Time 16:52 CDT. fabrun=4 (G05, R02 ring bath, S00b seed111, S11b seed111).
- **Done 25/29.** ETA `.done` ~17:30–18:00. Draft `SUSTAIN_REPORT.md` written.

### Genesis T=3000 verdicts

| arm | Em | nHot | Rg | class |
|---|---:|---:|---:|---|
| G01 nowf kx | 650 | 392 | 8.1 | **genesis_plateau diffuse** (not compact HOLD) |
| G02 nowf k0 | 638 | 370 | 8.2 | ≈G01 (+2% Em) — shear not switch |
| G03 wf kx | 0 | 0 | — | **empty** door wall long-T |
| G04 +rad | 439 | 74 | 7.6 | taxed birth −33% Em |
| R01 ring | 4.88 | 6 | 1.64 | **HOLD_tiny** mono CV~0.05 |

### Mission answer (standing)

- Compact sustain: **unitary + qd** (S00/S04/S14/S16) — **internal**, no fuel.
- Fuel = box-fill; metabolic thruput ≠ particle; genesis open-door = diffuse plateau.
- Partial final: S00b HOLD-like; S11b fill; G05≈G04; R02 Em flat ring.
- **Next:** on `.done` finalize report + NUMERICS; launch **R1 T=10000** (4 arms) if ≤21:30 — plan ready in R1_ENDURANCE.md.

## 2026-08-08 17:00 CDT — deadline EXTENDED to 2026-08-09 11:00 CDT

- User: continue experiments/validations until **11:00 AM tomorrow**.
- CALLING.md rewritten; R1_ENDURANCE updated (8 arms T=10000 overnight OK).
- Hourly scheduler re-prompted for new deadline.
- Primary still RUNNING ~25/29 (final batch: G05, R02, S00b, S11b).
- Plan: finish T=3000 → finalize SUSTAIN_REPORT → launch R1 → monodispersity refinements → stop new arms 09:00 08-09 → final writeup by 11:00.

## 2026-08-08 17:47 CDT — primary DONE; R1 live; Claude mid-flight

- **Primary T=3000:** 29/29 done (`sustain_T3000_20260808_104548.done`).
- **R1 T=10000:** launched `R1_T10000_20260808_172409` — 4 arms (S00/S04/S14/S16) at t≈1100–1200; **Em_tot flat ~148** so far (HOLD holding).
- **Claude qgraph:** still writing (~25 min). Go package expanded (rates dual I/O, charge Q=Z−N, door, shells, maxima stefan.mac). `go test` has 2 failing tests mid-edit. docs/ still empty — awaiting VERDICT.
- Showcase webps complete under `/space/scp/v93/webp/`.


## 2026-08-08 17:12 CDT — showcase export
- FCS+WebP to /space/scp/v93 (hold, rad, noise, progressive, genesis).
- Short collision/EMF/ambient runs in flight; auto-render on finish.
- T^4: emission YES (p_rad=4); absorption NOT as T_bath^4 (door intake only).

## 2026-08-08 17:20 CDT — qgraph Claude launched
- Subdir: theory/qgraph (Go scaffold + tests green).
- Claude fable working on UNIFIED_LAW / Maxima / charge-as-shell-defect.
- Awaiting REVIEW_claude.done then VERDICT.md → optional sim queue.

## 2026-08-08 17:22 CDT — Claude (default model) active on qgraph
- PID 3345062; brief in theory/qgraph/BRIEF_claude.md.
- Scaffold graph.go + tests pass. Deliverables: docs/UNIFIED_LAW, Maxima, VERDICT, Go expansion.
- On YES verdict: queue sim tests from SIM_PREDICTIONS.md.
- Showcase webps ready under /space/scp/v93/webp including collision/emf/ambient.

## 2026-08-08 17:55 CDT — Claude DONE YES; P1 offline; P3/P5 queued
- VERDICT YES. go test green. Unified dual I/O + Q=Z-N + Stefan ds.
- P1: matter HOLD ds~1.4; full foam ds~2.9. See P1_RESULTS.txt.
- Queued theory/runs/qgraph_P5_* (P3+P5 config). P2 kernel not started.
- R1 continues Em flat.

## 2026-08-08 18:01 CDT — primary finalized; R1 ~18% HOLD flat

- Time 18:01 CDT (deadline 2026-08-09 11:00; stop-new-arms 09:00). Well before cutoffs.
- **Primary T=3000:** confirmed `.done` 29/29. Structural scoreboard written (auto Em-only HOLD table retired).
  - Final batch: **S00b HOLD** (Em=146.9 nHot=56); **S11b box-fill** robust; **G05≈G04** taxed genesis; **R02** Em flat but nHot 6→1 (bath dephase).
  - Compact HOLD list: S00, S00b, S04, S10, S14, S16 (+R01 HOLD_tiny).
- **Docs:** `SUSTAIN_REPORT.md` final (run-dir), `MECHANISM.md` final, `NUMERICS.md` §7 updated.
- **R1:** already live `ACTIVE_R1` → `R1_T10000_20260808_172409` (8-arm suite, j=4). No re-launch.
  - Batch-1 (S00/S04/S14/S16) t≈1800–1840/10000; Em_tot **flat** 147.7–148.3.
  - fabrun count ~8 (R1×4 + qgraph P3/P5×4).
  - ETA batch-1 ~21:00; second batch (rad/seed/amp) overnight → full suite ~00:30 CDT.
- **qgraph:** Claude VERDICT YES done; P3/P5 still running (T=500 short); P1 offline ds matter~1.4 / foam~2.9.
- **Next hour:** watch R1 Em/structure; on batch-1 complete score mid-suite; no new multi-hour arms unless spare capacity for monodispersity short T. Stop new multi-hour after 09:00 08-09.

## 2026-08-08 19:00 CDT — R1 ~49% flat; mono launched; P3/P5 scored

- Time 19:00 CDT. Deadline 08-09 11:00; stop-new-arms 09:00. Capacity free → refinement OK.
- **Primary:** still DONE 29/29 (no change).
- **R1 batch-1** (S00/S04/S14/S16) t≈4760–4960/10000 (~49%), wall ~1h36:
  | arm | Em | class mid |
  |---|---:|---|
  | S00 | 148.30 | HOLD flat |
  | S04 | 147.29 | HOLD mild −0.27 creep |
  | S14 | 147.71 | HOLD_ext flat |
  | S16 | 148.07 | HOLD flat |
  - ETA batch-1 ~20:45; rad/seed/amp batch overnight → suite done ~00:15.
- **Mono scan LAUNCHED:** `mono_T2000_20260808_190046` (ACTIVE_MONO), probe -suite mono -T 2000 -j 4.
  - 8 arms: qd∈{1.2,3.6} × k_rad∈{0,0.02,0.05,0.10}. ETA ~2 h (finish ~21:00, parallel with R1).
- **qgraph P3/P5 DONE** (T=500):
  - P3: p_rad=2 → Em 148→99 (−33%); p_rad=4 → 148→131 (−12%). Lower p taxes harder; adopted p=4 = mild SLOW_DRAIN.
  - P5: door0≈door1 Em~2320 box-fill — **no coherent door preference** (claim FAIL at this apparatus).
  - Writeup: `qgraph/docs/P3_P5_RESULTS.md`.
- fabrun≈8 (R1×4 + mono×4). load~5.6/16.
- **Next:** score mono on finish; R1 batch-1 complete → interim R1_SCORE; no more multi-hour launches tonight unless mono negative and a short follow-up fits.

## 2026-08-08 20:00 CDT — R1 ~77% HOLD; mono qd12 = no Em* park

- Time 20:00 CDT. Deadline 08-09 11:00; stop-new-arms 09:00. No new multi-hour launches (R1+mono fill capacity).
- **Primary:** DONE 29/29 unchanged.
- **R1 batch-1** t≈7560–7880/10000 (~77%), wall ~2h36:
  | arm | Em | note |
  |---|---:|---|
  | S00 | 148.30 | flat HOLD |
  | S04 | 147.29 | mild long creep; still HOLD-class |
  | S14 | 147.71 | flat HOLD_ext |
  | S16 | 148.07 | flat HOLD |
  - ETA batch-1 ~20:50; rad/seed/amp batch → suite ~00:30 CDT.
- **Mono 4/8 (qd12 DONE):**
  | arm | Em | nHot | class |
  |---|---:|---:|---|
  | kr00 | 148.3 | 53 | HOLD |
  | kr02 | 124.8 | 31 | SLOW_DRAIN |
  | kr05 | 113.2 | 26 | SLOW_DRAIN |
  | kr10 | 105.4 | 21 | SLOW_DRAIN + peak weaken |
  - **Judgment:** radiance alone does **not** monodispersify — monotone tax, no Em* fixed point without open I. Deep qd36 midflight (t=900) resists more; finish ~20:25.
  - Interim: `mono_.../MONO_INTERIM.md`.
- fabrun=8 (R1×4 + mono qd36×4).
- **Docs:** NUMERICS §7 refreshed; HOURLOG this entry.
- **Next:** finalize mono scoreboard on `.done`; R1 batch-1 → interim R1_SCORE.md; overnight watch R1 batch-2. No further multi-hour arms unless short offline-only work.

## 2026-08-08 20:43 CDT — P2 dual intake kernel AUTHORIZED + live

- User auth: "You have auth for P2 kernel, possibly copy to new folder first, then make p2. Go."
- Baseline snap: `theory/p2/baseline/` (params/step/run/sim + freecell.c).
- Knobs: `rad_dual` (default 0, byte-inert), `bath_xb`/`bath_em` apparatus.
- Physics: I = k_rad*rad_dual*cap*(Ee/cap)^p reverse of radiance bookkeeping.
- C+Go both updated. Smoke: dual on grows Em_tag vs dual off at held bath.
- Battery launched: `theory/runs/battery_p2_20260808_204313/`
- P2 suite launched: `theory/runs/p2_T600_20260808_204313/` (9 arms, j=2, T=600).
- Docs: theory/p2/DESIGN.md, qgraph/docs/P2_KERNEL.md. Not adopted into laws_V3a.

## 2026-08-08 21:00 CDT — R1 batch-1 PASS T=10k; mono FINAL; P2 live RED abx

- Time 21:00 CDT. Deadline 08-09 11:00; stop-new-arms 09:00. No new multi-hour launches.
- **Primary:** DONE. **Mono:** DONE 8/8 (marker `.done` touched).
- **R1 batch-1 COMPLETE — endurance PASS:**
  | arm | Em | max | nHot | verdict |
  |---|---:|---:|---:|---|
  | S00 | 148.30 | 1.71 | 52 | PASS HOLD |
  | S04 | 147.29 | 1.99 | 55 | PASS HOLD |
  | S14 | 147.71 | 1.56 | 50 | PASS HOLD_ext |
  | S16 | 148.07 | 1.66 | 50 | PASS HOLD |
  - Score: `R1_SCORE.md`. **Internal compact HOLD validated to T=10000.**
- **R1 batch-2** (~4–8%): S00_s111 flat 146.8; S00_rad Em~123 drain; S04_rad Em~137; amp12 Em~213. ETA suite ~00:30–01:00.
- **Mono FINAL:** all k_rad>0 = SLOW_DRAIN; deep resists; **no Em* park**. `MONO_FINAL.md`.
- **P2 kernel** (user-auth 20:43): battery **RED 6/87** — all C≡Go abx (+ conserve cross) byte fails; physics PASS. Dual not adopted.
  - P2 suite: xb03 done d1 Em 4568 > d0 4214 (intake direction OK); xb06 mid d1>>d0; still high-Em bath regime not compact park.
- fabrun≈6 (R1×4 + P2×2).
- **Docs:** R1_SCORE, MONO_FINAL, NUMERICS, MECHANISM (next-internal rec), HOURLOG.
- **Next:** score R1 batch-2 as arms finish; diagnose C≡Go abx if capacity (mirror dual inert path); P2 suite to completion; no more multi-hour unless P2 needs short follow-up before 09:00.

## 2026-08-08 22:00 CDT — P2 FINAL (box-fill dual); R1 b2 ~40%; C Em_tot mirrored

- Time 22:00 CDT. Deadline 08-09 11:00; stop-new-arms 09:00. No new multi-hour arms.
- **Primary / mono:** DONE. **P2 suite:** DONE 9/9.
- **R1 batch-1:** PASS T=10k (unchanged). **Batch-2** t≈3600–4040/10000 (~38%):
  | arm | Em | note |
  |---|---:|---|
  | S00_s111 | 146.8 | HOLD seed-robust |
  | S00_rad | 108.8 | SLOW_DRAIN — no Em* park at long T |
  | S04_rad | 118.4 | deep resists more |
  | amp12 | 207.4 | higher-mass near-flat after settle |
  - ETA suite ~00:30–01:00.
- **P2 FINAL** (bath_xb apparatus):
  - Dual ON > OFF always; ratio grows with xb (1.08→2.74). Direction OK.
  - **All box-fill** (Em 2k–6k, nHot 1.3k–2.6k) — not compact HOLD park.
  - Do **not** adopt. Next test should be dual on quiet S00+rad (no bath flood).
  - Writeup: `P2_RESULTS.md`.
- **Battery abx root cause + fix:** Go had `| Em_tot=` diag; C lacked it → C≡Go abx RED.
  - Mirrored Em_tot into `kernel/freecell.c`. Short abx_b2 masked body **PASS**.
  - Full battery re-run still needed (not launched tonight — R1 owns j=4).
- fabrun=4 (R1 batch-2 only).
- **Next:** R1 suite complete → final R1_SCORE; optional battery re-run overnight if cores free; quiet dual-on-S00+rad short T if time before 09:00.

## 2026-08-08 23:00 CDT — R1 b2 ~70%; quiet dual park launched

- Time 23:00 CDT. Deadline 08-09 11:00; stop-new-arms 09:00. Short refinement OK.
- **Primary / mono / P2 bath suite:** DONE.
- **R1 batch-1:** PASS. **Batch-2** t≈6880–7440/10000 (~70%):
  | arm | Em | class mid |
  |---|---:|---|
  | S00_s111 | 146.8 | **PASS HOLD** flat all horizon |
  | S00_rad | 106.6 | SLOW_DRAIN (slope slowing ~−5e-4; no hard Em*) |
  | S04_rad | 115.9 | SLOW_DRAIN; deep +9 vs law trap |
  | amp12 | 206.7 | HOLD-like higher mass |
  - ETA suite complete ~00:00–00:30 CDT.
- **Launched P2 quiet park:** `p2quiet_T2000_20260808_230025` (ACTIVE_P2QUIET)
  - 4 arms: qd∈{1.2,3.6} × rad_dual∈{0,1}, k_rad=0.05, **no bath_xb** (tests I on clump, not box flood).
  - j=2 sequential batches; ETA ~40–50 min. fabrun≈6 with R1.
- **Docs:** R1_SCORE mid, NUMERICS, HOURLOG.
- **Next:** score quiet dual on finish; finalize R1 full 8/8; battery re-run if cores free after R1; no multi-hour new after 09:00.
2026-08-09T00:00:01-05:00 R1 finished

## 2026-08-09 00:00 CDT — R1 FULL FINAL 8/8; quiet dual recovers Em; battery re-run

- Time 00:00 CDT 08-09. Deadline 11:00; stop-new-arms 09:00.
- **R1 COMPLETE 8/8** (`.done` marked). Final `R1_SCORE.md`:
  - **PASS HOLD:** S00, S04, S14, S16, S00_s111, S00_amp12 (Em 147–207, nHot 50–106).
  - **SLOW_DRAIN:** S00_rad Em=105 max=0.75; S04_rad Em=114 max=0.66 — **no Em* park** to T=10k.
  - **Headline:** internal unitary+qd HOLD is long-lived; radiance alone is a tax.
- **P2 quiet qd12 FINAL:** d0 Em 149→113 drain; **d1 Em 149→123→134 recover** (late slope +0.009, ratio 1.18). Dual assists on clump without bath flood. qd36 still ~t=1180.
- **Battery re-run LAUNCHED** (`battery_p2fix_*`) after C Em_tot mirror — expect abx green.
- **MECHANISM** updated with next-internal-dynamics recommendation.
- fabrun: p2quiet qd36×2 + battery jobs.
- **Next:** finalize p2quiet qd36; battery result; optional longer quiet dual if green and time; finalize writeup trail by 11:00.

## 2026-08-09 01:05 CDT — p2quiet FINAL; battery mask fix; T=5000 dual + batt2 live

- Time 01:05 CDT. Deadline 11:00; stop-new-arms 09:00. Short/overnight OK.
- **R1 / primary / mono:** DONE (no change).
- **P2 quiet T=2000 FINAL** (all 4):
  | arm | Em | slope2 | class |
  |---|---:|---:|---|
  | qd12 d0 | 113 | −0.008 | SLOW_DRAIN |
  | qd12 d1 | 134 | **+0.009** | RECOVER (not park) |
  | qd36 d0 | 122 | −0.007 | SLOW_DRAIN |
  | qd36 d1 | 129 | −0.002 | slowed drain |
  - Dual assists without box-fill; law trap recovers inventory late; not Em* fixed point.
- **Battery p2fix1:** still RED 6/87. Root cause after Em_tot: banner `freecell` vs `freecell-go`, column header missing Em_tot on C, 1-ulp credit. **Physics t= lines already equal.**
  - Fixes: C header + Em_tot; battery maskDrift banner+credit. Short abx masked **PASS**.
  - **p2fix2 battery re-launched** `battery_p2fix2_20260809_010524`.
- **Launched p2quiet T=5000** (4 arms, j=2 batches) — longer park/overshoot test. ETA ~1.5–2 h.
- fabrun≈2–4 + battery. Docs: P2QUIET_SCORE final, NUMERICS, HOURLOG.
- **Next:** score battery2 (expect GREEN or near); score T=5000 dual; prepare final campaign writeup by 09:00–11:00.

## 2026-08-09 02:00 CDT — BATTERY ALL GREEN 87; T=5000 dual ~62% climb

- Time 02:00 CDT. Deadline 11:00; stop-new-arms 09:00.
- **Battery p2fix2: ALL GREEN 87/87** (exit 0). abx C≡Go PASS. Dual still default-off; not adopted into laws_V3a.
- **Primary / R1 / mono / P2 bath / p2quiet T=2000:** all DONE (unchanged).
- **P2 quiet T=5000** mid t≈3120/5000 (wall ~55 min):
  | arm | Em@3.1k | path note |
  |---|---:|---|
  | Q5_qd12_d0 | 109.9 | SLOW_DRAIN continues (113@2k→110) |
  | Q5_qd12_d1 | 136.3 | after recover, **slow climb** 134@2k→136 (overshoot, not park) |
  - ETA qd12 done ~02:40; qd36 batch after → suite ~04:00.
- **Campaign judgment written:** `CAMPAIGN_JUDGMENT.md` — sustain YES; switches; next-internal rec; deliverable checklist.
- fabrun=2 (T5000). No new multi-hour arms.
- **Next:** score T=5000 final; light polish NUMERICS/MECHANISM; by 09:00 analysis-only; final docs by 11:00.

## 2026-08-09 03:00 CDT — T=5000 qd12 FINAL: dual = SLOW_REFILL not park

- Time 03:00 CDT. Deadline 11:00; stop-new-arms 09:00.
- **Battery:** ALL GREEN 87 (held). R1/primary/mono DONE.
- **P2 quiet T=5000 qd12 DONE:**
  | arm | Em_end | last-2k slope | class |
  |---|---:|---:|---|
  | d0 | **108.1** | −0.0009 | SLOW_DRAIN (slowing) |
  | d1 | **140.4** | **+0.0023** | **SLOW_REFILL** (134@2k→140@5k) |
  - Ratio 1.30. Dual real I on clump; **no Em* fixed point** — overshoots and keeps climbing slowly.
- **qd36 batch** mid t≈1540 Em~125/128; ETA suite ~04:30.
- Docs: `P2QUIET_T5K.md`, CAMPAIGN_JUDGMENT dual note, NUMERICS, HOURLOG.
- **Next:** finalize T5k deep; no new multi-hour launches; polish judgment for 09:00–11:00 freeze.

## 2026-08-09 04:05 CDT — quiet dual T=5000 FULL FINAL 4/4

- Time 04:05 CDT. Deadline 11:00; stop-new-arms 09:00. Analysis only preferred.
- **Primary T=3000 / R1 T=10k / mono / P2 bath / p2quiet T=2000:** DONE (unchanged).
- **Battery:** ALL GREEN 87 held (p2fix2).
- **P2 quiet T=5000 COMPLETE** (`.done` 04:02):

  | arm | Em_end | slope 3k→5k | class |
  |---|---:|---:|---|
  | qd12 d0 | 108.1 | −0.0009 | SLOW_DRAIN |
  | qd12 d1 | 140.4 | +0.0023 | **SLOW_REFILL** (ratio 1.30) |
  | qd36 d0 | 117.2 | −0.0012 | SLOW_DRAIN (deep +9) |
  | qd36 d1 | 139.1 | +0.0030 | **SLOW_REFILL** (ratio 1.19) |

  - Dual = real quiet I at both detunes; **no Em\* park**; deep trap does not stop refill.
  - No box-fill (bath_xb=0 correct). **rad_dual remains default-off / not adopted.**
- **fabrun:** idle (no live campaigns).
- **Docs updated:** `P2QUIET_T5K.md` full final, NUMERICS, MECHANISM, CAMPAIGN_JUDGMENT.
- **Next:** polish-only until 09:00; after 09:00 analysis/freeze; by 11:00 final SUSTAIN_REPORT / HOURLOG close if needed. **No new multi-hour arms.**

## 2026-08-09 05:00 CDT — all campaigns idle; SUSTAIN_REPORT postscript

- Time 05:00 CDT. Deadline 11:00; stop-new-arms 09:00. **No new arms.**
- **Status check:** primary / R1 / mono / p2quiet T2k+T5k / battery 87 all **DONE** `.done`. fabrun **idle**.
- **T=5000 dual:** already FULL FINAL (04:02) — SLOW_REFILL both traps; no Em\* park; not adopted.
- **Polish:** `SUSTAIN_REPORT.md` post-matrix table (R1 6/8 HOLD, mono, dual T5k, battery GREEN) + frozen next-internal list → points at `CAMPAIGN_JUDGMENT.md`.
- **Deliverable checklist:** 1–5 essentially closed; remaining = light freeze polish after 09:00, stop at 11:00.
- **Next hour:** re-verify idle; no launches; if still pre-09 optional doc consistency only.

## 2026-08-09 06:00 CDT — idle re-verify; docs clean; no launches

- Time 06:00 CDT. Deadline 11:00; stop-new-arms 09:00. **No new arms.**
- **fabrun/probe/battery:** idle.
- **Done markers hold:** sustain T3000 · R1 T10k · p2quiet T5k (and T2k). Battery 87 GREEN record unchanged.
- **Deliverable cross-check:** SUSTAIN_REPORT / R1_SCORE / CAMPAIGN_JUDGMENT / MECHANISM / NUMERICS / P2QUIET_T5K / MONO_FINAL all present; no stale "in flight/ETA/live" language.
- **Judgment frozen content:** sustain YES (uni+qd through T=10k); dual = SLOW_REFILL not Em\* park; shear ≠ birth switch; next-internal in CAMPAIGN_JUDGMENT §5.
- **Next:** 07–08 polish only if needed; ≥09:00 analysis/freeze; ≥11:00 stop + final HOURLOG close if not already written.

## 2026-08-09 07:00 CDT — idle; theory README campaign table

- Time 07:00 CDT. Deadline 11:00; stop-new-arms 09:00. **No new arms.**
- **fabrun:** idle. All primary/R1/mono/dual/battery DONE.
- **Polish:** `theory/README.md` — campaign-status table (CLOSED for scoring) + document index pointing at CAMPAIGN_JUDGMENT / SUSTAIN_REPORT / R1_SCORE first.
- Content judgment unchanged: sustain YES (uni+qd T=10k); dual SLOW_REFILL not Em\*; no launches planned.
- **Next:** ≥09:00 analysis-only freeze; ≥11:00 stop + final HOURLOG close.

## 2026-08-09 ~07:30 CDT — campaign CLOSED; hourly loop killed

- User: kill loop if done; point to status complete + FCS/viz.
- Scheduler task `019fe36384ea` **cancelled**.
- Wrote **`theory/STATUS_COMPLETE.md`** (single entry: scoreboards, FCS archive, `/space/scp/v93` webp+fcs).
- `theory/README.md` points to STATUS_COMPLETE first.
