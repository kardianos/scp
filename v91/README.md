# v91 — RADIANCE: stable mass from steep emission (+ the coherent channel)

**THIS IS THE SINGLE ENTRY POINT. An agent starting v91 work reads this
file top to bottom, then follows "NEXT STEPS" at the bottom. Everything
needed is in this directory or explicitly linked.**

Opened 2026-08-04, user-authorized ("Yes. Write up the v91 version...").
Prepared by the closing v90 session; the baseline is BUILT and VERIFIED
(kernels byte-identical C↔Go, battery green at the inert defaults —
see `runs/BATTERY.log`), and the first law-candidate knob is wired but
INERT. **No new-concept experiment has been run. That is the v91
programme, and it starts here.**

## The one-sentence thesis

Give the dense sector a steep graded radiance (emission rising with
fullness below cap, on a clock that does not slow) and — second stage —
a coherent amplitude channel, and this substrate should select interior
flux fixed points: **stable masses with derived sizes, forged in its
own density wells** (the FORGE conclusion, `carried/FORGE_EVIDENCE.md`).

## What v91 carries (all COPIED into this tree; originals frozen)

| path | what | provenance |
|---|---|---|
| `kernel/freecell.c` | **C kernel of record** — full v90 apparatus surface (slit/rings/blob2/xsec/p1/sect/convtag/grad/tag_r) **+ the v91 radiance candidate, inert at k_rad=0** | copy of `v90/kernel/freecell.c` + one guarded term |
| `fab/` + `cmd/fabrun` | the Go kernel A/B experiment (same surface incl. radiance; byte-identical to C on the bath — re-verified at v91 open) | copy of `v90/fab`, imports rewritten to `scp/v91` |
| `cmd/battery` | the battery harness; law-purity check now ALSO pins `k_rad=0 p_rad=4 rad_clock=0` (the inert defaults are part of purity until a table change is adopted) | copy of `v90/cmd/battery` + 1 header line |
| `cmd/volview` | the volumetric viewer (`-i`, `-avg`, `-info`, `-follow`) | copy |
| `laws/laws_V2g.cfg` | the standing table, verbatim copy (canonical stays `v89/battery/laws_V2g.cfg`) | copy |
| `THEORY.md` | the standing law + the two modifications specified + back-references | written at open |
| `carried/FORGE_EVIDENCE.md` | WHY v91 exists: the balance-curve kill (cap-wall step), dark lensing, forging-follows-density | copy of `v90/FORGE.md` |
| `carried/S2_CHANNEL.md` | the coherent-channel entry measurement + design memo (stage 2) | copy of `v90/S2.md` |
| `carried/REALITY_AUDIT.md` | the correlation-to-reality audit v91 exists to improve | copy of `v90/REALITY.md` |
| `carried/VERIFY_PROTOCOL.md` | the C↔Go A/B rules and byte claims | copy of `v90/VERIFY.md` |
| `carried/P2_SCHEDULER.md` | the local-clock prototype record + production checklist (needed when boxes grow) | copy of `v90/P2.md` |
| `Makefile`, `go.mod` | build; `go.work` at repo root already includes `./v91` | adapted |

**Build:** `cd v91 && make all` (gcc + go ≥ 1.27; `make freecell` for C
only). **Battery:** `./battery` from `v91/` — must end `ALL GREEN`.

## The law change, precisely (implemented, inert)

In pass 6, on each cell's beat, AFTER the evaporation block:

    if k_rad > 0 and Em > 0:
        x      = Em/cap                      # holdings (decision R-D1)
        comp   = rad_clock ? 1 : w2/w2e      # = det: beat-stretch compensation
        demand = k_rad · cap · x^p_rad · comp
        → atoms_fire/clamp (grains ε, shared F-credit with evap)
        → routed as evaporation: Em−, field += dr−bs, space += bs

Keys `k_rad` (0), `p_rad` (4), `rad_clock` (0) print on their own
header line; the battery pins them at these defaults (purity). With
`k_rad=0` the kernel is **byte-identical** to v90's (verified: C≡Go
bath diff empty; full battery green). Wiring verified live once:
`exp=pair bath=0 T=60 pair_x0=1.0 k_rad=0.5` → `rad=4.584987`,
identical in both kernels, drift ~1e-14. Decision points R-D1..R-D4
and the full rationale: `THEORY.md` §2.1.

## Rules (carried, non-negotiable)

1. **The ratchet:** every kernel/table change runs the FULL battery
   before commit; passing claim-defending experiments join the suite;
   bars sharpen by measurement, never soften to pass; a green test
   leaves only by explicit user decision. A LAW ADOPTION (k_rad > 0
   becoming the table) is a conscious re-gauging event: every bar that
   moves gets an explicit entry in the campaign doc.
2. **Working rules:** derive the criterion before the sweep; print the
   config beside every result; controls verifiably run; magnitudes not
   booleans; pre-register predictions BEFORE first run; state where
   every model's state lives (no background).
3. **Pre-v89 ban unchanged and permanent** (root `CLAUDE.md`).
4. **C is the kernel of record; Go is the verified A/B experiment.**
   New apparatus lands in C first; the abx battery experiment gates
   the pair.
5. **Instrument lessons that will bite you** (from the FORGE/E1 runs):
   the pair-picker moves with d*(x) — compensate `pair_doff` to pin
   the probed cell across an x sweep; sub-grain demands accrue credit
   silently (a zero can mean "below one atom over this T", so size T
   accordingly); single-seed ANGULAR ratios are foam-speckle-dominated
   (λ/dmin 3–7) — ledger claims are the seed-robust ones; the sector
   meter is cylindrical in 3D.

## NEXT STEPS — the R campaign (radiance), in order

Create `v91/RADIANCE.md` as the campaign doc (pre-register before
running; FORGE.md is the format precedent). Then:

**R0 — baseline re-verify (5 min).** `make all && ./battery` → ALL
GREEN (the log at open is committed; re-verify on your machine).

**R1 — the balance curve with the knob on.** Rebuild FORGE E1
(fixed-cell sweep: `exp=pair bath=1 freeze_geo=1 pair_x1=0`,
compensated `pair_doff`, `convtag=1`, noisy ambient) over a
(k_rad, p_rad, rad_clock) grid × fullness x ∈ [0.1, 1.05]. Deliverable:
intake(x) and outtake(x) curves per parameter point; find the interior
crossing x* and its slope. Pre-register: candidate A with rad_clock=0
produces a crossing for k_rad in some band; rad_clock=1 (ablation)
reproduces the anti-Stefan pathology. Select (k_rad, p_rad) by the
crossing quality, not by taste. Also record the comb-resonance
structure of outtake(x) — the emission-line probe.

**R2 — stability of real objects.** At the selected point: single
blob, ring6, UUD, composite rings (nv=6) in vacuum AND in bath, long T.
Bars: lifetime ≥ 10× the v90 ceilings (ring death t≈1900, composite
dissolution t≈480); perturbation return to x*; intake=outtake in the
ledger (the flux-machine criterion). The v90 meters are all present
(convtag, sect, p1, Em_tag diags).

**R3 — size selection (the WHY answer).** Measure the surviving
objects' bond lengths against the prediction d* = πC/w(x*). If sizes
land on the derived rung: particle sizes are wavelengths of
flux-balanced clocks — write it as a gated bar. Enumerate the
(interval, x*) zoo the comb admits.

**R4 — forging closure.** Re-run FORGE E3 (density well + beam,
`grad_r0/grad_frac/tag_r` apparatus is in the kernel) at the selected
point: condensation in a well must now SATURATE into a finite stable
object instead of a sink. This closes the star-formation loop.

**R5 — the ratchet reckoning.** Full battery at the selected
(k_rad, p_rad): list every bar that moved; for each, decide (with the
user) sharpen/re-gauge/reject. Only then does laws_V3r become a table
(`laws/laws_V3r.cfg`) and the purity line change.

**Then stage 2 — the coherent channel** (`carried/S2_CHANNEL.md` §2):
design decisions C0, implementation C1, acceptance surface C2 (e3b
coherent, the fifth survives, composite EMF radiates, exclusion
precursor). Radiance-stabilized objects are its test bodies.

**Parallel/opportunistic:** the no-law quantitative targets from
`carried/REALITY_AUDIT.md` C7 (DS visibility-vs-geometry curve;
Born-rule KS distance) — they upgrade the reality audit regardless of
how R goes. The fission lane (THEORY.md §2.3) stays design-only until
R lands and the user re-authorizes.

## Status at open (2026-08-04)

| item | state |
|---|---|
| tree carried + imports rewritten | done |
| radiance candidate A implemented, inert | done (both kernels) |
| C≡Go byte identity at defaults | **verified** (bath diff empty) |
| wiring check k_rad>0 fires identically | **verified** (rad=4.584987 both) |
| baseline battery | **see `runs/BATTERY.log` — must read ALL GREEN (93 bars)** |
| RADIANCE.md campaign doc | not started — the first v91 agent writes it (pre-registration first) |
| laws_V3r table | does not exist yet — R5 creates it |
