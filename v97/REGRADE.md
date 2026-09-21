# REGRADE — exact-ledger meters + leads re-grade (2026-08-13, PLAN step 1′)

**Verdict: H1 still EMPTY, 0/10 leads arms — now on honest books.**
All three pre-registered predictions (PLAN §1) hit. The Φ story is
corrected: the door barely fires and fires one-sided (I ≪ O, a slow
drain); the old ⟨I⟩ ≈ ⟨O⟩ near-balance was hop noise. The promo
orientation discovery *survives and strengthens* where the door
actually fired both ways (polar 0.36 → 0.88 on p2_space).

## 1. What changed (offline `hunt/` only; no kernel, no bars moved)

1. **Exact ledger is schema-level, not opportunistic.** When the FCS
   carries `promo_i/promo_o`, `MeasureConvert` uses them on every
   frame, zeros included (`convert.go`). The old per-frame switch fell
   back to the hop-shaped heuristic on every quiet snap. The heuristic
   remains the fallback for streams without the columns, and is still
   computed as a side-channel (`Ih/Oh`) so its inflation stays
   measurable.
2. **`sep` of an empty book is undefined, not zero:** averaged only
   over frames where both I and O fired (matches the existing `polar`
   convention of skipping den≈0 frames). `nBook/nLive` is printed so
   sparsity is visible. Direction acknowledged: this can only raise
   sepBar; the Φ clause remains the liveness gate.
3. **Reporting (bars unchanged):** Φ/Em also as a per-t.u. rate with
   the frame Δt (the 0.015 bar is snap-dependent — charter debt,
   REVIEW_grok_next Q5c); `min(⟨I⟩,⟨O⟩)` and the I/O ratio alongside
   ⟨min(I,O)⟩, separating mean imbalance from asynchrony; the I/O
   source (`exact`/`heuristic`) named in every clause line.
4. **New H0 class `exact`** (`runs/h0.sh`, now a script): ledger
   present on every frame; graded I/O ⊆ box ledger; heuristic
   inflation reported. **Regression tests** pin the graded totals to
   the kernel's own RESULT lines (`hunt_test.go`:
   p2_space cond = 0.524917, evap+rad = 14.184963; p0_ctl must have
   no ledger). `go test ./hunt/` green.

## 2. H0 recheck — boring, as required

`runs/regrade/H0_recheck.out`: all 8 standing classes PASS unchanged
(archive streams have no ledger → heuristic path identical). New:

| stream | inflation measured by the meter itself |
|---|---|
| p2_space | heuristic ΣI = 20.44 vs ledger 0.5249 — **I × 38.9**, O × 1.33 |
| p2_hot | I × 13.9, O × 1.71 |

Grok's hand check (~40× on p2_space) is now a standing instrument.

## 3. Leads re-grade — `runs/regrade/H1_exact.out`

| arm | clauses old→new | Φ/Em exact | I/O | book frames | sep clause |
|---|---|---|---|---|---|
| p0_ctl | 2/5 → 2/5 | (0.0036 heuristic — no ledger at promo=0) | 1.73 | 26/26 | FAIL |
| p1_slot | 2/5 → **1/5** | 0 (⟨I⟩ = 0 in window) | 0 | 0/26 | FAIL |
| p2_space | 3/5 → 3/5 | **7.3e-5** (was 1.8e-3) | **0.070** | 1/26 | PASS 2.36 / polar **0.88** |
| p3_close | 2/5 → **1/5** | 0 | 0.073 | 0/26 | FAIL |
| p2_hot | 3/5 → 3/5 | 7.9e-4 | 0.456 | 17/26 | PASS 1.91 / 0.60 |
| p2_noise | 3/5 → 3/5 | 1.1e-3 | 2.81 | 23/26 | PASS 3.10 / 0.26 |
| p2_ring | 1/5 → 1/5 | 0 (door never fired) | — | 0/26 | FAIL |
| p2_dual | 3/5 → **2/5** | 0 | 0 | 0/26 | FAIL (old sep PASS was hop fiction) |
| p2_jm | 2/5 → **1/5** | 0 | 0 | 0/26 | FAIL |
| p2_tight | 0/5 → 1/5 | 7.7e-4 | 0.44 | 1/26 | PASS (gained by the defined-book convention) |

**Predictions (PLAN §1) scoreboard:** exact Φ/Em drops — **HIT**
(24× on p2_space); I/O imbalance appears — **HIT** (I = 0 exactly in
four arms' best windows); sep/polar survives on promo>0 where the
door fired — **HIT, strengthened** (hop noise had been *diluting*
the polarization: 0.358 → 0.880).

## 4. Corrections to the LEADS.md-level reading (dated here)

- "Every promo>0 blob with door traffic passed I/O spatial split"
  **narrows to**: where the door fired both ways (p2_space, p2_hot,
  p2_noise, p2_tight), the split is real and strong; in p1_slot /
  p2_dual / p2_jm / p3_close the best window had **zero intake
  events** — their old sep/polar PASSes were heuristic fiction.
- "Φ/Em ~10× below bar" **worsens to ~20–200×** on exact books, with
  the named cause: the door's intake is nearly dead (cond = 0.52 over
  the whole p2_space run) and one-sided. The Φ failure is an
  **apparatus/intake question** (field supply, workfn, atom grain),
  not a pairing-noise question — exactly what PLAN step 2 probes.
- p2_dual (rad_dual=1) has ⟨I⟩ = 0 in its best window: the spontaneous
  dual did not deliver intake *where the motif lives* at this bath.

## 4b. Battery events (2026-08-13, authorized kernel edits)

With the user's sign-off ("kernel edit accepted") two `fab/run.go`
changes went in, and the full battery ran per the ratchet:

1. **Composed helix seed:** the ring `kx` tilt now applies AFTER
   `seed_mw` and composes (`th2 = m·φ − kx·x`). kx-only and
   seed_mw-only seeds are byte-unchanged; kx+seed_mw was previously a
   silent no-op (rg2d ≡ rg2c byte-identical — `A5_RING.md`).
2. **Pre-existing parity break found and fixed:** the first battery
   run came back RED 6/87 — all six C≡Go byte-parity bars. Proven
   pre-existing by stash-and-rerun at HEAD
   (`runs/regrade/BATTERY_pre_edit.log`, same 6 RED): the leads
   commit added an unconditional `# v97 promo leads` header that the
   C kernel never prints, and the full battery was not re-run after
   (the green `BATTERY.log` is 12:29, the leads runs 13:17+). Fix:
   the header prints only when `promo != 0`, restoring byte-equality
   at defaults.

**Battery after both edits: ALL GREEN 87**
(`runs/regrade/BATTERY_helixseed.log`).

## 5. Named caveats (kept visible, not fixed silently)

- p2_space's sep PASS rests on **1/26** book frames (polar 0.88 rests
  on ~17 frames with any traffic). Whether the sep bar should require
  a minimum book-frame count is a **bar sharpening → user decision**;
  nBook is printed so nobody mistakes thin evidence for thick.
- With I = 0, `polar` degenerates to outtake-side asymmetry (it reads
  the oriented outtake deposit, not an I/O book). The clause is
  guarded by nBook > 0, so no false pass; noted for readers of the
  polar column.
- promo=0 controls can never carry the ledger (fab writes the columns
  only at promo>0), so control rows stay heuristic-graded and are not
  directly comparable. Making the columns unconditional is a
  physics-inert `fab/` edit that would change default-output schema —
  **needs sign-off**, listed in PLAN §4.
- Archive comparability: the sep defined-book convention could shift
  archive numbers slightly (heuristic books are dense, so almost
  none expected). Full archive re-check **done**:
  `runs/regrade/H1_archive_recheck.out` — **0 inhabited / 259 files**
  (135 qualifying, 124 short, 0 errors). H1.md's EMPTY verdict
  survives the meter update. (The original 203-file list is not
  byte-recoverable — dirs gained files since; this sweep is a
  superset.)

## 6. Post-consult queue results (PLAN §5b, 2026-08-13 evening)

**Item 0 — dest-only re-grade: the pre-registered KILL band hit.**
`-dest` (replayed promoPick fore/aft picks — `hunt/dest.go`; unit
test pins scope) on the old arms: rg2c Φ/Em = **0.0012** (inside the
0.000–0.005 "there was no Φ pass" band), rg2f 0.0057, rg2e late 0.
`runs/regrade/DEST_regrade.out`. **Confound A confirmed in full: the
"first Φ-clause PASS" was carried by the dilate halo.** Even beyond
the true door destinations, ~68% of the dilated intake was bath
contamination. A5_RING's headline is corrected accordingly; `-dest`
is now the default grading scope for promo streams in this
programme's write-ups (`-dilate` remains a diagnostic).

**Item 3 — `pJ` integrated hop-current meter: built, H0 GREEN,
battery GREEN.** `pj_meter=1` (Go, default 0, byte-inert) books the
per-cell fd (half per endpoint; Σcells ≡ p1fd exactly). H0 on three
fresh e3b-style seeds at L=24 dt=0.01: identity ratio
1.000126–1.000207 (f32 floor), fd seed-robust **+x** (316–360) while
the snapshot bilinear reads **−78 to −92 ± ~250 — the wrong sign,
not just noise**. The H0 current-direction SKIP is retired for
pj-instrumented streams. `runs/pj/e3b_s*.{fcs,log}`. Battery after
the kernel edit: **ALL GREEN 87** (`runs/regrade/BATTERY_pjmeter.log`).
The H1 J clause now grades on the integrated channel when the
columns exist (same bars 0.25/0.80; snapshot bilinear remains the
labeled fallback).

Items 1/2/4 (3-seed replication, quench park probe, promo=3 ring)
run at the user-directed scale L=24 dt=0.01 snap_every=400 →
`runs/rep/`, plus `rf_pj` (pj-instrumented helix) for the first
honest J-clause read. Results: `REPLICATE.md`.
