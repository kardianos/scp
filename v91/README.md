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

## Status after the R campaign (2026-08-04, same day — R0–R5 EXECUTED)

**`RADIANCE.md` is the complete record** (pre-registrations, ~470
runs, results, the reckoning). Headlines: the interior fixed point
EXISTS and is gated (x̂* = 0.62 ± 0.02 at the selected point
k_rad=0.05 p_rad=4 rad_clock=0; quarter-power law in k and p;
ablation and Go-port verified; conservation ≤ 1.35e-15 throughout);
the bath GLOWS (nucleated proto-matter, supply-throttled); radiance
TAXES v90-era hoard-objects (t_half 80–140 vs k0's 260–510) whose
ruins then track the ambient's fixed point; and **structures cannot
hold the balance without the coherent channel** — measured from both
sides (bath pile-pump-out; vacuum-at-x* bond-gate break under
thermal atom jitter). The flux-machine interior (S2) is thereby
measured load-bearing for stable mass. **laws_V3r was NOT created**
— agent recommendation in RADIANCE.md §6: adopt radiance + coherent
channel TOGETHER when the coupled acceptance surface passes; the
per-bar reckoning (11/93 movers, all classified) is in §5 awaiting
user sign-off. The battery at defaults stays ALL GREEN (93) —
re-verified after the harness gained the `-extra` override flag.
Next agent: stage 2 (S2 coherent channel), with radiance as its
test harness, per THEORY.md §2.2 and `carried/S2_CHANNEL.md`.

## Status after the CANTUS campaign (2026-08-05, user-directed: "atoms are NOT cells")

**`CANTUS.md` is the complete record** (pre-registrations §1–§3,
three design rounds §3.2–§3.4, results §4). The user-directed
candidate B — a superimposed field fitted to the cells, holding its
own gauges, simulating the overall harmonic lock — is IMPLEMENTED in
both kernels and INERT at defaults (`k_cant=0 k_tune=0` = byte-
identical to the pre-cantus kernels, verified against the archived
binary; battery ALL GREEN 93 after every change; C≡Go byte-equal
with the candidate firing). Final form v1.1i: per-BOND gauge memory
`sgg` (the chord chart is link-borne), Kuramoto lock on matter
clocks correcting only the differential closure error, within-mode
retune current on holdings-memory deviations (COE-metered, p1-clean).

Headlines: **(1)** the wall at the fixed point is GEOMETRIC — the
unison rung exceeds contact at x̂\* (the R-campaign's x62 vacuum arm
was dead at birth, slots=0); only p+q ≥ 5 intervals fit — the comb
admits only CHORDS at the balance. **(2)** The V2g churn bath is a
persistent-closure medium (66% of cells hold time-averaged gate
> 0.5 at zero force — the null-meter) ⇒ every self-growing coherence
rule ignites a bath-wide Kuramoto transition and suppresses the glow
50× (v1) / 10× (v1.1): **closure statistics cannot distinguish bond
from churn; only IDENTITY can — the exchange registry is measured
load-bearing for coherence itself**, not just exclusion. **(3)** The
coupled lock+tune mechanism holds pairs (90–94% vs 13% control under
the confound; ret > 1 at frozen-gauge bound in the honest medium)
and each half alone dies — but in the HONEST medium (v1.1i
instrument, bath exactly V2g+radiance) self-tracking gauges buy only
1.4–1.8×, and the frozen-full-strength upper bound splits by
topology: unison rings stay at 1.4× (winding conservation is the
measured killer inside the coherence sector) while **the UUD triad
reaches 8.3× with its retention curve turning back up — the chord
the geometry demands is the object that nearly holds.**
**(4)** H5 reckoning: 83/93 with both candidates live; movers = the
R5 classes + one cantus signature (coherent matter ABSORBS above the
band — the lock raises the cross-section). **(5)** OUTLOOK E-A:
the integrated nv=48 object holds ret 0.95±0.02 to **t=50,000**
(≥100× v90 ceilings, no decay trend). laws_V3 remains NOT adopted;
all knobs inert at defaults; every decision awaits the user per the
ratchet. Next agent: the exchange registry (S2 §2.3) is now the
single measured door — design identity-carrying transfers with the
cantus as carrier apparatus and the v1.1i instrument as the
honest-medium harness. **[2026-08-05: the user opened this lane —
"Use the CANTUS experiment as a mechanism of discovery... Do open
exchange registry lane and do go for a UUD chord." The campaign doc
is `REGISTRY.md`; standing concepts were adjusted inline the same
day (RADIANCE §4.8 correction, OUTLOOK §2/§4, THEORY §2.2,
carried/S2_CHANNEL.md addendum).]**

## Status after the REGISTRY campaign (2026-08-05, user-directed: "do open exchange registry lane and do go for a UUD chord")

**`REGISTRY.md` is the complete record** (pre-registration §0–§3,
three registered design rounds §3.4/§3.5, results §4, verdict §5).
The apparatus is IN both kernels, inert at defaults (`reg_tau=0
reg_gate=0 reg_f0=0` purity-pinned; -O0 byte-equal to the
pre-registry source including drift; battery ALL GREEN 93 twice;
C≡Go byte-equal with the meter and both gate forms firing).

Headlines: **(1)** the churn bath is a persistent-EXCHANGE medium —
per-line balance EMERGES with the ledger window (bath ρ median
0.49/0.68/0.79 at τ=10/30/100; 92% of slots flowing): the v90
"detailed-rate balance without per-line balance" gains its scale,
and balance is a unison instrument (pair ρ 0.95) while chord bonds
run a standing directed current (UUD net ≈ 0.75·gross — the
flux-machine circulation observed per link for the first time).
**(2)** The lock manufactures its own gate evidence: all four
registry-gated self-growing combos ignite the bath (economy
−60/−83%, f0-insensitive) — the third measured ignition of the
statistical-gate family — and the same statistics are impotent as
maintenance (honest-instrument UUD 152 vs control 140). **(3)** The
UUD chord with ASSERTED identity does not die: ret 0.50–0.63 to
t=5000, bonds alive, books balanced, mean load at the 0.466 sweet
spot — superseding the CANTUS 8.3× bound. **(4)** The door in final
measured form: **identity must be born with the matter and carried
by its exchanges — ontological, not statistical** (three independent
closures). Next lane recommended in §5: parcel-carried birth
identity (gid at the conversion door, no background), acceptance
target = the frozen-bound chord. Nothing adopted; decisions are the
user's.

## Status after the COMBINE campaign (2026-08-06, user-directed: "do combined experiments to see which effects survive in the presence of others")

**`COMBINE.md` is the complete record** (registration §1–§4 committed
before any arm ran; per-experiment results §5.1–§5.7 written at each
harvest; synthesis §6). Config-only on the standing binary; battery
re-verified after (`runs/BATTERY_combine.log`); nothing adopted.

Headlines — the attribution map (§6): **(1)** chord stability is
carried by the coupled lock+TUNER (both required: no-lock dies t=36,
frozen-lock-without-tuner t=48 — gauges die with their bonds) and by
nothing else — it survives radiance removal, the cool medium, and
coexistence. **(2)** The doublet law belongs to laws_V2g itself — it
survives lock removal at frozen-grade tightness (1.6%), coexistence,
and even the starving loser (3.7%); the per-event residual is the
spectroscopy's invariant. **(3)** METABOLISM AND LIGHT are
radiance-era phenomena: without the radiance drive the door
collapses ~200×, the chord's tagged books stop at evap = 0.000
exactly, and the object holds as an economically static statue —
stability≠metabolism, cleanly split. **(4)** Two-body physics
measured for the first time: isolated coexistence preserves every
effect with pooled luminosity EXACTLY 2.00× (per-object conserved),
zero mutual force at range 7 (sep slope −6.8e-5±4.2e-5/t.u. ≈ 0, as
π-flat demands), a SHARED box balance locus (three chords ±0.001 —
an operating point, not a constant); contact coexistence adds
**competitive starvation** — probabilistic (1 of 2 contact pairs),
killed by feeding-connection loss (1 vs 18–33 intake events), read
live by the spectrometer (loser's line walks blueward), with no
door-to-door coupling (ILAG stays at the 1.4× bath common-mode).
**(5)** CO-RL closes the no-kernel-work path: the lock on the
embedded nv48 is the registered NULL (1.67×, fold and starvation
schedule untouched) while 26/48 gauges hold rock-stable through the
whole death — assertion preserves SKELETONS, not living contents;
the nv6 winding wall does not generalize (fold regime ≠ winding
regime). **(6)** One interference found and it was a RECORD error:
ASTRO 6.5.4's UUD brightness baseline 0.0230 unreproducible →
dated correction (UDD is 1.1–1.5× brighter absolute, sub-linear per
emitter; ratios 2.1–2.6× stand); zero apparatus-on-apparatus
interference anywhere. Instrument notes: kernel diag `sep=` is
wrap-broken for straddling objects (display-only; fcs true-sep = the
standing method); door-net drain warns, intake disconnection kills.

## Status after the ASTRO campaign (2026-08-05, user-directed: "make it so" on the far-field + spectroscopy probes)

**`ASTRO.md` is the complete record** (pre-registration §0–§3 with
dated shakedown corrections, results §4, verdict §5). A MEASUREMENT
campaign: config-only arms on the standing binary, no kernel or law
change; battery re-verified at defaults after (`runs/BATTERY_astro.log`).

Headlines: **(1)** pad 11 (g4 reopened) is CLOSED with the mechanism:
around every stable body the medium's transport potential
π = Es + s_disp·(Em+Ee) is FLAT at ±1×10⁻⁵–2×10⁻⁴ (blob, chord,
floors) because **a stable object's space cycle closes at its own
conversion door** — footprints are contact-local (< 1 cell
screening), the corpse control carries none, and the only resolvable
medium disturbance anywhere is the radiative dM/dt transient of a
DYING object (r^−0.18 — luminosity leaves no space far field). B6
stays ABSENT with a measured prerequisite: a far field needs books
that run THROUGH the medium while the object lives (the S2/identity
currency again). **(2)** B4 moves — the first species spectra: the
UUD chord's parameter-free emission doublet (mean residual 2×10⁻⁴,
seed-robust, load-tracking), a D line 30–60× over a spectrally-dark
bath (species detectable by light alone), the anti-Stokes gap at the
law constant to 0.3%, metabolism visible as eat-at-U/shine-at-D, and
a redward death spectrogram on the dying control. **(3)** New
measured limitation: the nv=48 flagship's stability is
VACUUM-CONDITIONAL — first-ever real-medium embedding folds it at
settle and starves it (ret 0.038 @ t=1500) while the frozen chord
lives its full 5000 embedded: embedded stability exists only below
the fold scale. NOTHING ADOPTED; the §5 recommendation (unchanged:
the parcel-identity lane, now also far-field-motivated) awaits the
user.

## Status after the pad + integration campaigns (2026-08-04, same day)

User-directed sequels, all recorded in-tree: **`pad/`** (34 crazy
ideas, up to 3 believer-rounds each — `pad/RESULTS.md` is the index;
apparatus added under the ratchet: QATOM stream with cell/Em tags +
`qatom_every`, `RESULT credit` line, `cmd/fcsdump`, battery `-extra`;
ALL GREEN 93 re-verified after each change). **`pad/INTEGRATION.md`**:
the six working pads COMPOSE — an nv=48 low-curvature ring at the
ambient's fixed point holds Em, bonds, balanced books, and a
two-sided-selected size (d = 1.62 ± 0.01) for ≥ 5000 t.u. at every
ambient (≥10–60× the v90 ceilings), while curved controls die and the
k0 control merely sleeps — acceptance items 2/3/4/7 now MEASURED for
large straight objects; "structures cannot hold the balance" is
CONDITIONAL (small/curved matter still needs S2). **`OUTLOOK.md`**:
ability / limitations / next experiments and theory.
