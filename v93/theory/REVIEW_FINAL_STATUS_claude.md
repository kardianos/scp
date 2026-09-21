# External review — `theory/FINAL_STATUS_ADOPTION.md`

**Reviewer:** independent (claude), read-only
**Date:** 2026-08-09
**Scope:** primary doc + SUSTAIN_REPORT / R1_SCORE / MONO_FINAL / P2QUIET_T5K /
P2_RESULTS / qgraph UNIFIED_LAW §3.3 / P1_RESULTS.txt / P3_P5_RESULTS /
`fab/step.go` / `fab/params.go` / `kernel/freecell.c` / `laws/laws_V3a.cfg` /
`cmd/battery/main.go` / `theory/STATUS_COMPLETE.md`. Also read (not in the
brief, but load-bearing): `v93/RECKONING.md`, `v93/FACE_C.md`,
`v93/RING_DNLS.md`.

No simulations re-run. Two cheap re-reads of already-written campaign data
(`results.json`, `*.em.csv` in the campaign run dirs) were used to check the
document's own claims against its own recorded meters.

---

## Overall verdict: **MOSTLY SOUND**

The document is unusually careful *as an adoption ledger*. Every numeric I
spot-checked against the five cited reports matches. The code→theory mappings
are accurate down to the line number. The reject list is honest, the "do not
adopt / do not silently flip" posture is ratchet-correct, and all 21 inventory
items really are covered in §2.

It is **not** sound in one place, and that place is the headline: the positive
claim ("**Sustain is real** … compact HOLD … channel of record") is stated more
strongly than the campaign's own recorded meters support, and it is not
reconciled with `v93/RECKONING.md`, the committed same-week v93 verdict that
this exact channel produces a **transient**, not matter. Fix those two things
and the document is SOUND.

---

## What checks out (verified, not assumed)

- **Numerics.** Items 1, 2, 3, 4, 7, 9, 12, 13, 15, 16, 17, 18, 21 all match
  their sources: quiet dual 108.1 / 140.4 / 117.2 / 139.1 and ratios 1.30 /
  1.19 and late slopes +0.0023 / +0.0030 (`P2QUIET_T5K.md`); bath arms Em
  2229–6258, nHot 1316–2616, dual ratios 1.08→2.74 (`P2_RESULTS.md`); mono
  ladder 148.30 / 124.81 / 113.15 / 105.36 (`MONO_FINAL.md`); R1 rad 105.41 /
  113.90 with Em_max 0.75 / 0.66 (`R1_SCORE.md`); S01 77.6/0.52/4, S02
  148.8/0.34/0, S03 148.8/0.80/10, S05 109.9, S06 119.7, G01 650 vs G02 638,
  G03 empty (`SUSTAIN_REPORT.md`); P1 matter d_s 1.14–1.59, foam 2.88–2.94
  (`P1_RESULTS.txt`); P3 −33% vs −12% and P5 2327 vs 2321 ≈1.00
  (`P3_P5_RESULTS.md`). No fabricated or drifted numbers found.
- **Code grounding.** Line citations land exactly: `fab/step.go:1405` radiance,
  `:1432` dual, `:1509` bath, `:202` `det = 1+QDetune*x`, `:414`/`:482` pass U.
  The rad_dual description (gate `KRad>0 && RadDual>0 && Em>0 && Ee>0`, formula,
  atom fire/clamp, reverse bookkeeping — `fa` shrink, Es reverse-backs,
  `Em += di`, `radInTotal`) is faithful. The workfn mapping
  (`econdI = Em>=WfFloor ? 0 : WfFar`) is faithful. `laws_V3a.cfg` values
  (k_rad 0.05, p_rad 4, rad_clock 0, wf_on 1, wf_floor 0.01, wf_far 99,
  q_detune 1.2) are quoted correctly. `bath_xb` is correctly described as
  injecting from outside the ledger — the code does set `Ee = bath_xb*cap` with
  no matching debit.
- **Theory quotes.** Theorems 3.4 (`x* = xb`) and 3.5 (`u⁴ − g u³ − 1 = 0`,
  `u*>1`, restoring, `Φ*>0`) match `UNIFIED_LAW.md` §3.3 in substance.
  VERDICT's two admitted gaps ("no force law for Q", "graph is given, not
  selected") are quoted fairly.
- **Battery.** `theory/runs/battery_p2fix2_20260809_010524/BATTERY.log` ends
  `ALL GREEN (87 bars)`. Confirmed.
- **Pointers.** `/space/scp/v93/{fcs,webp}` exist (14 FCS, 38 WebP);
  `hold_S00_em.webp` and `hold_S00_uni_quiet.fcs` resolve. P1's
  `hold_S00_uni_quiet.fcs` (559 cells) is the same object as `S00_uni_base.fcs`
  (559 cells) — internally consistent with STATUS_COMPLETE.
- **Coverage.** Items 1–21 each have a §2 section. Nothing in the inventory is
  orphaned.

---

## Findings

### HIGH

**H1 — The document does not reconcile with `v93/RECKONING.md`, which reaches
the opposite verdict on the same channel.**
`RECKONING.md` (2026-08-08, HEAD commit 723068d, titled *"Face C is a TRANSIENT,
not matter"*) ran three gates on the **same** unitary `amp_nat` + `q_detune`
substrate and concluded the condensate is (a) not stable (peaks ~t=300, then
monotone Em_max decay + PR/rms spread through T=3000), (b) phase-blind /
phase-random, (c) not particulate (one connected condensate). Its NOT-to-do list
names, explicitly, *"adoption of amp_nat>0 into a law table … gated on a
confirmed stable bound state, which does not yet exist."*
`FINAL_STATUS_ADOPTION.md` cites neither RECKONING nor FACE_C nor RING_DNLS, and
its §4 freeze opens with "**Sustain is real**". A reader of the v93 tree sees a
direct contradiction with no argument offered. The two may well be reconcilable
(RECKONING's R1 scored a *spontaneously condensing broad packet* on Em_max/PR/rms;
this campaign scores a *seeded blob* on Em_tot/nHot/Em_max) — but that
reconciliation has to be written, not left implicit.
*Severity: high — it is the document's central claim.*

**H2 — "Compact HOLD" is contradicted by the campaign's own recorded PR meter.**
The probe harness records participation number `PR` in `results.json`
(`prog.PR`) but neither `R1_SCORE.md` nor this document reports it. Reading it
out of the campaign's own files, for the identical arm/trajectory at t=3000 vs
t=10000:

| arm | PR @T=3000 | PR @T=10000 | Δ |
|---|---:|---:|---:|
| S00 | 297.2 | 371.5 | **+25%** |
| S04 | 304.6 | 395.6 | **+30%** |
| S16 | 298.0 | 405.5 | **+36%** |
| S14 | 371.0 | 430.2 | +16% |

Em_max on the same arms, windowed means from `R1_*.em.csv`: S00 2.03 (t 0–1k) →
1.71 (t 9–10k); late/early ratios across the PASS arms are **0.67–0.82**; nHot
falls ~66→51 on S00. That is a slowly de-peaking, monotonically spreading
object — the same signature RECKONING used to declare the physics transient,
now extended to T=10000. "Em flat 147–207 … PASS HOLD" is true only of `Em_tot`.
*Severity: high — the word "compact" is not supported by the meter that measures
compactness, and that meter was recorded.*

**H3 — The load-bearing pass bar is near-tautological in the quiet arms.**
R1 bars are `|slope|≲0.02` and `Em_end≳100` **on `Em_tot`**. In the quiet arms
(`k_rad=0`, `noise_amp=0`, presence-gated door) `Em_tot` is a conserved ledger
quantity: in `R1_S00.em.csv` it is *literally constant* at 148.298 (and `Ee` at
0.4059) from t≈2000 through t=10000, and `results.json` reports `Cond=0`,
`Rad=0`, `Evap=0.4668` total over 10 000 t.u. So the flat-Em bar is measuring
conservation, not sustain. The genuinely structural bars are `nHot`, `Em_max`,
and (unreported) `PR`/`rms` — and those move.
*Severity: high — affects how much weight §1/§4 can carry.*

### MEDIUM

**M1 — "Channel of record" over-reaches the standing v93 position.**
§1(2) and §2 item 5 call unitary+linearize+Strang the "channel of record" for
"particle-like HOLD". Given RECKONING's gate, the defensible phrasing is
"the best measured transport algebra for long-lived energy concentration"
— and the qualifier "not particulate / phase-blind" belongs next to it. The
document's own hedges ("not adopted", "needs sign-off") are right; the noun
phrase is not.

**M2 — Zero battery coverage for every new/campaign knob, and no measured
C≡Go at nonzero.**
`grep` over `cmd/battery/main.go` returns **no occurrence** of `amp_nat`,
`rad_dual`, `bath_xb`, `amp_door`, or `hop_order`. Consequences the document
should state plainly: (a) the campaign channel (`amp_nat>0`) has **no ratchet
protection at all** — a future kernel edit could break it silently; (b) the abx
C≡Go bars run at defaults, so the **C-side pass U and C-side dual code paths are
exercised by no cited evidence**, while §2 items 1/5/12 say "C: same placement"
in a way a reader will take as verified parity. All campaign arms are Go
(`fabrun`) only. This is a registered debt, not an error — but it should be
registered.

**M3 — The coded dual is not the exact time-dual of radiance (clock
normalization).**
Both blocks use `comp = P.W2 / s.w2e[i]` — the **dense** clock — even though the
dual's source is the field. The atom grain *is* correctly swapped
(`eFi = A0eff·w1e/2π`, `eDi = A0eff·w2e/2π`), so the asymmetry is specifically in
the per-time compensation. The document's "same (k, p)" is true; "spontaneous
dual of radiance" is stronger than the code delivers. This is exactly the kind
of factor that decides whether a balance root exists, so it belongs in item 1's
redesign note.

**M4 — `rad_dual` is a strength multiplier, so Theorem 3.4's root shifts.**
As coded, `I = k·rad_dual·cap·xb^p`, hence the g=0 balance root is
`x* = rad_dual^{1/p}·xb`, coinciding with `xb` only at `rad_dual=1` (which is
what was run). §2 item 1 quotes `x* = xb` without the caveat, and §1 treats
`rad_dual` as an on/off flag. Harmless for the runs; misleading for the redesign.

**M5 — "Real I channel" in the quiet apparatus is (at least partly)
re-absorption of the object's own emission.**
With `bath_xb=0` the box is closed and total energy conserved (`E0=2829`, of
which `Em~148`); the quiet-dual refill from 108→140 necessarily draws on the
`Ee`/`Es` inventory that radiance itself deposited. "I exceeds O" is correct
bookkeeping, but it is not evidence of a *sustainable* intake, and the late
positive slope is equally consistent with the closed box relaxing toward a
global `x* = xb`. `P2_RESULTS.md` §5 already flags that `Em_tot` is the wrong
observable and asks for core-cell `x` vs `xb`; that caveat did not survive into
FINAL_STATUS. Recommend reporting `Ee_tot`/`Es_tot` next to `Em_tot`, and
stating the boundedness explicitly.

**M6 — Two levers missing from the 21-item inventory.**
(a) The **topological ring** (`RING_DNLS.md` A.1: W=+2 exact, which RECKONING
calls *"the ONE genuinely stable object"*) appears only as SUSTAIN's `R01/R02`
HOLD_tiny rows and is not an inventory item, despite being the only retained
structure in the programme. (b) The **J-dependent nonlinearity / dynamical
Hodge star** — RECKONING's *recommended next face* — is absent, while §4's
"Next" recommends dual-balance + `g` instead. Recommending a different next face
than the standing recommendation is legitimate; doing so without naming it is not.

**M7 — Item 10 ("`g` is the highest-upside unfinished theory") omits its
negative priors.**
`g` is uncoded *for the radiance dual*, correctly stated. But the tree already
carries two coded coherence-bias mechanisms — `amp_drv` (phase-current bias on
the dense want, `fab/step.go` ~400–412) and `amp_door` — and both have negative
records (v92 closed the feedback hypothesis in three forms; P5 ratio ≈1.00).
"Highest-upside" is defensible; stating it without those priors is optimistic.

### LOW

- **L1** — Item 4: "deep retains ~+9–15 Em". `MONO_FINAL.md`'s own table gives
  **+7.3 / +9.3 / +10.9** at k_rad = 0.02 / 0.05 / 0.10. Should read "+7–11".
  (`MONO_FINAL.md`'s own "+10–15" is likewise off and should be fixed at source.)
- **L2** — Item 13: "S07–S15" as a range sweeps in **S10** (quiet HOLD) and
  **S14** (HOLD_ext), which are not box-fill. Enumerate:
  S07, S08, S09, S11, S11b, S12, S13, S15, S17, S19.
- **L3** — §3 index drops the `theory/` prefix on every run dir
  (`runs/sustain_T3000_…`, `runs/R1_T10000_…`, `runs/p2quiet_T5000_*`,
  `runs/battery_p2fix2_…/BATTERY.log`). Actual paths are under
  `v93/theory/runs/`. `STATUS_COMPLETE.md` has them right — the two disagree.
- **L4** — Item 5: "C default `hop_order=1`". Both kernels default to 1
  (`fab/params.go:250`, `kernel/freecell.c:408`). Say "both kernels".
- **L5** — Name collision: "R1" is the T=10000 endurance suite here and the
  endurance *gate* in `RECKONING.md`. Disambiguate (e.g. "R1-endurance (theory)"
  vs "reckoning-R1").
- **L6** — Item 5: "Em flat 147–207" — the minimum is 146.85 (`R1_S00_s111`).
- **L7** — Radiance and dual are described as "pass-6" blocks; both sit inside
  the `beatFire` branch, so I and O are *both* beat-rate gated. Worth a
  half-sentence, since it is part of why the two share a rate scale.
- **L8** — Ratchet bookkeeping: this campaign added code to both kernels and
  produced claims, but contributed **no new battery bar**. That is acceptable for
  default-0 knobs, but the ratchet's "experiments that pass and defend a claim
  are added to the suite" is unaddressed. Register it as debt.

---

## Adoption hygiene (brief item 5)

**Consistent with the ratchet, with one gap.** "Do not adopt `rad_dual` into
V3a", "keep default 0", "no silent `amp_nat` flip, needs user sign-off + a
dedicated battery panel" are all correct, and they agree with `RECKONING.md`'s
NOT-to-do list. The gap is **M2**: the ratchet protects only what the battery
gates, and the battery currently gates *none* of the campaign's physics — the
green 87 proves inertness at defaults, nothing more. The document should say
so in one sentence rather than letting "87 GREEN" sit adjacent to the sustain
claim, where it reads as support for it.

Note also that `hop_order=1` (Strang) is *already* the shipped default in both
kernels. It is byte-inert at `amp_nat=0` and was battery-gated at adoption, so
this is fine — but a document that freezes adoption judgments should say that a
schedule default was flipped, rather than mentioning `hop_order` only as a knob.

---

## Corrections list (shortest path to SOUND)

1. Add a subsection reconciling with `RECKONING.md` / `FACE_C.md` — state the
   apparatus difference (seeded compact blob vs spontaneously condensing broad
   packet), and say which verdict governs. **(H1)**
2. Report `PR` (and `rms`/`FTag`, all already in `results.json`) for the R1 PASS
   arms, and either drop "compact" or defend it against PR 297→371 / 305→396 /
   298→406. **(H2)**
3. Restate the R1 pass bar: note that `Em_tot` is conserved in quiet arms
   (constant to 3 decimals from t≈2000) and promote `Em_max`/`nHot`/`PR` to the
   load-bearing bars. Report the Em_max late/early ratios (0.67–0.82). **(H3)**
4. Downgrade "channel of record" → "best measured transport algebra for
   long-lived concentration"; carry "phase-blind, not particulate" alongside. **(M1)**
5. Add one sentence: no battery bar exercises `amp_nat>0`, `rad_dual>0`,
   `bath_xb>0`, or `amp_door>0`; C≡Go is verified only on the `amp_nat=0`
   surface; all campaign arms are Go-only. **(M2, L8)**
6. In item 1, note the dense-clock `comp` asymmetry and that
   `x* = rad_dual^{1/p}·xb` as coded. **(M3, M4)**
7. In item 1/§4, note that quiet-dual refill is bounded by the closed box's own
   `Ee+Es`, and carry `P2_RESULTS` §5's core-`x`-vs-`xb` caveat forward. **(M5)**
8. Add the topological ring and the J-dependent/Hodge-star face to the
   inventory (or say why they are out of scope). **(M6)**
9. Cite `amp_drv`/v92-feedback and P5 as priors next to the `g` recommendation. **(M7)**
10. Mechanical: "+7–11 Em" (L1); enumerate box-fill arms (L2); add `theory/`
    prefixes (L3); "both kernels" for `hop_order=1` (L4); "146.85–206.75" (L6);
    note the beat gating (L7).

---

## One caution about this review

Findings H2/H3 rest on re-reading meters already written by the campaign
(`results.json` `prog.PR`, `*.em.csv`) — no new physics was run. They show the
*positive* claim is weaker than stated; they do **not** show it is false. A
slowly de-peaking, slowly spreading object that survives 10 000 t.u. with a
conserved ledger and no fuel is a real and reportable result — it is just not
"compact HOLD", and it is the same object RECKONING already classified. The
comparative claims of the campaign (unitary ≫ additive: PR 297 vs 504,
Em_max 1.93 vs 0.53; gated τ melts: PR 1721, Em_max 0.34; `q_detune=0` diffuses:
PR 900) are unaffected by any of this and stand as measured.

---

*End of review.*
