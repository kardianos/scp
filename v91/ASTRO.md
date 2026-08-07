# ASTRO — the far field and the first species spectra (pre-registered)

Opened 2026-08-05, user-directed ("Okay, make it so. Proceed." — on the
recommendation: re-open g4 at the newly stable bodies + species
spectroscopy of the chord). This campaign is a MEASUREMENT campaign:
no law candidate, no kernel change, config-only arms on the standing
apparatus. It continues pad `11-g4-reopened.md` (R1: 2-cell null at the
±0.002 floor; R2 blob left pending with `runs/pad/w2_blob.fcs` already
on disk) and targets two rows of `carried/REALITY_AUDIT.md`:

* **B6 — no far field (ABSENT).** g4 closed honestly in v89: a leaking
  blob's space flux is all mass-rate bookkeeping; the far field "awaits
  a stable particle's internal space cycle." That prerequisite now
  exists: the nv=48 ring (ret 0.95 @ t=50,000, books 0.87–0.96 at
  100–200× sealed-door metabolism) and the frozen-identity UUD chord
  (ret 0.50 @ t=5000, books 29/31, standing per-link current
  net ≈ 0.75·gross, REGISTRY §4). The confound of the old null (huge
  net dM/dt) is gone: these bodies run throughput at net ≈ 0.
* **B4 — no spectra of species (ABSENT), partially moved.** The pads
  measured the CAVITY line (peak w=1.90, anti-Wien, two-family
  anti-Stokes); INTEGRATION measured a bonded object's emission line.
  What has never been measured: a species FINGERPRINT — an object
  whose spectrum has more than one line because the object has more
  than one voice class. The UUD chord is that object by construction
  (two pitch classes, wD/wU = 2/3 exactly at seed).

State declaration (no background): all state lives in the cells and
slots of the running kernel; the FCS stream and QATOM prints are pure
output (consume no RNG, alter nothing — verified property of the
apparatus, FCS.md). Analysis is offline over `fcsdump` output.

## 1. Bodies and meters — the mass ladder

The one design axis the pad-11 round 1 result forces: footprint vs
MASS. The 2-cell pair (tag Em ≈ 3.1) showed no radial Es structure
above the ±0.002 meter floor. So the campaign runs a ladder:

| body | tag Em (measured) | medium | stream |
|---|---|---|---|
| UUD chord, frozen instrument | ~3.5 → 1.7 | bath L=16, amp 0.5 | new run (A-F1) |
| nv=48 ring, EMBEDDED | ~54 at seed | bath L=36, amp 0.15 | new run (A-F2) |
| blob (bondless condensate) | ~180–190 | bath L=16, amp 0.5 | EXISTS (`runs/pad/w2_blob.fcs`, T=2000) |
| 2-cell pair (anchor, done) | 3.1 | bath L=16 | pad-11 R1: NULL at ±0.002 |
| no object (floor) | 0 | bath L=16 / L=36 | `runs/pad/m1_glow_k005.fcs` + new L=36 arm |

Meters, all existing: FCS cell frames (`t i x y z r es em ee xload
tag`) → radial profiles Es(r), Em(r) about the tag centroid (angular
average is the seed-noise reducer; time average over frames is the
churn reducer); the CONVTAG books (object-attributed cond/evap/rough/
backs); `# QATOM t dir w e i Em` events (dir=DF: dense→field
emission — rough/evap/radiance atoms at the source's dense pitch;
dir=FD: field→dense condensation at the field pitch) with per-event
source cell id — object cells are the last ids (tri: NC−3..NC−1,
verified against the FCS tag column).

Note on what F-arms can and cannot see: the FCS link record carries no
`swl`, so the direct per-slot space current is not dumped. The F-arm
measurands are therefore the PROFILE (the g1-analog footprint and its
tail shape) and its STEADINESS, not the pointwise flux; continuity
(shell dEs/dt vs object books) is the cross-check where a profile
resolves.

## 2. Closed-form expectations (written before any new run)

* **Pressure balance.** Pass S equalizes π = Es + s_disp·(Em+Ee)
  across live slots (s_k=0.06, s_disp=0.3). A standing mass is a
  standing π-excess; if equalization were complete and unscreened, the
  surrounding medium would hold an Es DEPRESSION whose integral tracks
  s_disp·(Em+Ee)_object. [Instrument correction, 2026-08-05, from the
  sanctioned shakedown BEFORE any physics run: the bath-class medium
  is THREE-dimensional (2700 cells in 16³, dbar 1.41; the tri object
  sits in the z=8 mid-plane) — the earlier 2D wording was wrong.] In
  a 3D medium the steady unscreened monopole tail is 1/r — reality's
  own far-field form: ΔEs(r) = a + c/r in the far zone. A
  churn-screened footprint is instead short-ranged: ΔEs(r) ~
  exp(−r/λ_s) with λ_s of order dbar. The discriminator between
  "footprint" (g1-analog, near field) and "far field" (B6) is the
  tail SHAPE, not the depth. Shells: spherical about the tag centroid
  (compact bodies); cylindrical about the ring axis for the nv=48
  near zone.
* **Sign ambiguity registered honestly:** the glow's space-return
  (pad 10: +0.8% uniform Es lift at k005) is a GLOBAL offset, and the
  object is also a local GLOW ANOMALY (its books run hotter than
  bath). The near-object sign could be negative (pressure push-out) or
  positive (local cond starvation of space refilled from around). The
  bars below score magnitude and shape; the sign is a finding either
  way.
* **Emission lines are pitch atoms.** Every DF event's atom pitch is
  the firing cell's dense pitch w2e = w2/(1 + q_detune·x_cell). The
  chord's two voice classes therefore predict a DOUBLET with zero fit
  parameters, per time window, from the diag's own xUDD column:
  w_U(t) = 2.9/(1+1.2·x_U(t)), w_D(t) = 2.9/(1+1.2·x_D(t)). At seed:
  2.171 / 1.447 (ratio 3:2, the fifth). At the frozen chord's t=5000
  state (x = 0.06/0.05/0.58): ≈ 2.70 / 1.71. The bath's own emission
  band sits at ≈1.90 (pad 06) — the U line rides ABOVE it, the D line
  BELOW it, both drifting with load: the line is a load meter.
* **The absorption family** (FD events) sits at field-pitch atoms
  w1e = w1/(1+q_detune·x); emit/absorb pitch ratio AT THE SAME CELL is
  w2/w1 = 1.7576 exactly (det cancels) — the anti-Stokes gap is a law
  INVARIANT and doubles as an instrument check.

## 3. Arms, predictions, bars (registered BEFORE first run)

Configs print in every log header; seed 20260802 primary; the frozen
chord F1 arm gets one replicate seed (314159). All arms
k_rad=0.05 p_rad=4 rad_clock=0 (the selected radiance point — the
medium that maintains the bodies); registry/cantus knobs only where
the body needs its instrument frame (the frozen chord), exactly the
REGISTRY goal-arm configs. Logs + streams → `runs/astro/`.

* **A-F0 (no run — analysis of existing streams).** Blob Es(r,t) from
  `w2_blob.fcs` (+ `w2_blob2.fcs` as the two-body variant), floor from
  `m1_glow_k005.fcs`/`m1_glow_k0.fcs`. This is pad-11 R2, finally.
* **A-F1 (chord).** The REGISTRY frozen arm + streams:
  `exp=tri tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 k_cant=1
  k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10
  qatom_every=1 snap_every=500` (10 t.u. cadence — shakedown-sized,
  500 frames); ctl twin (no cantus/registry keys)
  identical otherwise. The ctl chord dies by t≈150 — its stream is the
  CORPSE arm (spectrum and footprint of a dying species).
* **A-F2 (embedded nv=48 — first contact with a medium).** `exp=ring
  ring_n=48 ring_x=0.47 bath=1 noise_amp=0.15 convtag=1 L=36 T=1500
  qatom_every=200 snap_every=500`. The carve rule leaves bath OUTSIDE
  clear=14.15 in an L=36 box (far zone r ≈ 14.2→17). This arm is
  double-duty: (a) the mid-ladder far-field body; (b) the first test
  of the flagship object against a real cell medium (INTEGRATION only
  ever ran it with all bath carved; its "ambient" was the noise
  drive). FLOOR twin: same box, `exp=bath`, no object, T=800.
  [Shakedown sizing, before physics: 3D at L=36 is ~30,700 cells and
  ~6.6 µs/cell-step single-threaded — T=1500/800 replaces the 2D-era
  T=3000 guess; the settled-window frames are what the bars read.]
* **A-L0 (bath spectrum baseline).** `exp=bath bath=1 noise_amp=0.5
  qatom_every=1 T=2000` at k005, L=16 — the medium's own two-family
  spectrum at full statistics. (Chord spectra come from A-F1's own
  QATOM streams; no separate runs.)
* **Shakedown protocol:** T≤100 smoke runs are permitted to verify
  snap/QATOM mechanics and wall-time ONLY; no bar may be read from a
  smoke run.
* **Log convention (one log path, one writer, two grades):** the
  committed `runs/astro/*.log` are QATOM-stripped (diag + RESULT
  rows); the full event-grade originals stay on disk as
  `*.log.raw` (gitignored, like the `.fcs` streams) — the derived
  event tables `spec_*.txt` and profiles `prof_*.txt` are committed
  and reproduce from the raws via `spec.awk`/`prof.awk` verbatim.

**Pre-registered predictions.**

1. Chord footprint: BELOW the floor (same mass class as the pad-11
   pair null) — registered as an expected null anchoring the ladder.
2. Ring and blob footprints: resolve above 3× floor, depth monotone
   in tag Em.
3. Tail shape where resolved: screened-exponential (λ_s ~ 1–2 dbar)
   is the honest default guess; a measured ln-r tail is the upset
   that moves B6 — either outcome is scored, neither is assumed.
4. The chord doublet: YES at the closed-form positions (±2%, the
   pad-23 same-beat systematic allowed), U above / D below the bath
   band, positions tracking xUDD per window.
5. Embedded nv=48: survives contact (ret ≥ 0.5 at T=3000) — 60%
   confidence; bath-crush death is the informative alternative
   (would make the flagship's stability vacuum-conditional, a
   sharpened OUTLOOK §2 limitation).

**Bars (magnitudes, not booleans).**

* **P-A1 (ladder):** report ΔEs(contact..r) depth per body vs tag Em;
  "resolved" = |ΔEs| > 0.006 (3× the ±0.002 pad-11 floor) sustained
  over ≥ 20 time-averaged frames and surviving the angular average.
* **P-A2 (shape):** where resolved, fit the far zone both ways
  (a + c/r vs A·exp(−r/λ_s), 3D forms per the §2 correction); report
  c (or λ_s) with residuals; claim "monopole tail" ONLY if the 1/r
  fit wins on residuals AND c's sign matches the measured object
  books. This is the B6 bar.
* **P-A3 (steadiness):** the resolved profile at matched windows
  early/late at net-dM/dt ≈ 0 — steady means late depth within 30% of
  early depth while books stay balanced (the "not bookkeeping" test).
* **P-A4 (attribution):** depth tracks throughput (books gross) or
  net (books in−out)? Report both correlations across the ladder —
  this is the Bondi-duality question pad 11 asked.
* **P-L1 (doublet):** chord-cell DF pitch histogram shows two
  resolved lines at the closed-form positions (±2% per window);
  resolved = each line's peak ≥ 5× the bath continuum in the same
  pitch bin, from ≥ 30 object events per window.
* **P-L2 (anti-Stokes invariant):** per cell class, (DF pitch)/(FD
  pitch) = 1.7576 ± 1% — the instrument-check bar.
* **P-L3 (species distinguishability):** line separation vs measured
  linewidth: |w_U − w_bath|, |w_D − w_bath| > 2σ_line — the "you can
  tell WHAT is shining" bar, the B4 mover.
* **P-L4 (death spectrogram):** the ctl chord's lines drift along
  w2e(x(t)) as it dies — report tracking residual; the line as a
  load/flight meter on a dying species.

**Reality-audit mapping (registered):** P-A2 pass at a stable body ⇒
B6 moves ABSENT → measured (first far-field candidate; g2/lensing
re-runs become motivated). P-A2 fail everywhere ⇒ B6 stays ABSENT
with a SHARPER statement: even a stable metabolizing mass with
standing internal current does not source a far field — the space
sector needs the coherent/identity construction too (S2-currency,
same door as everything else). P-L1+P-L3 pass ⇒ B4 gains its first
species row (STRUCTURAL → the first quantitative species statement:
line positions parameter-free). P-L1 fail ⇒ the emission machinery is
bath-collective, not species-borne — a new measured limitation.

## 4. RESULTS

### 4.1 A-F0 — the blob (pad-11 R2, closed from the standing stream)

Analysis: `runs/astro/prof.awk` over `runs/pad/w2_blob.fcs`
(151 frames, window t ∈ [500, 2000], tag-centroid shells, medium
cells only) vs the object-free floor `runs/pad/m1_glow_k005.fcs`.
Numbers: `runs/astro/prof_blob.txt`, `prof_floor16.txt`.

* **The medium around the 110-Em blob is in COMPLETE pressure
  equilibrium.** π(r) = Es + 0.3·(Em+Ee) is flat at 1.0563 ± 0.0002
  across the whole resolvable zone (r 4.75 → 8.75; the blob's tagged
  cells occupy r < ~4.5). The visible Es structure (0.966–0.972,
  anti-tracking Em at the s_disp factor) is pressure balance, not a
  footprint. The floor bath is equally flat at its own level
  (π = 1.0430 ± 0.0002). Meter resolution ±0.0002 — 10× sharper than
  pad-11 R1's ±0.002 Es floor, because π is the transport potential
  (pass S moves space on π differences ONLY): a flat π is a direct
  null on sustained space transport, not just on density structure.
* **Continuity bound.** |∇π| < ~1e-4/unit caps the per-slot sustained
  current at s_k·dt·O(1)·1e-4 ≈ 4e-6/t.u.; summed over the ~10³
  slots crossing a mid-zone shell: ≲ 5e-3 units/t.u. of
  medium-carried space current. The blob's conversion door cycles
  cond 2873 / evap 3336 / backs 436 over 2000 t.u. (≈ 3 units/t.u.
  gross, net −36 = 1.2% of gross), and the whole box's Es is steady
  to −2.6e-5/t.u. in the window. **The object's space cycle closes at
  its own door (backsplash refills what cond pulls, pad-10's 0.3%
  book balance), so the medium is never asked to carry a net current
  — the Bondi in/out duality cancels pointwise at the conversion
  door, not as separated accretion/wind zones through the medium.**
* Scored: P-A1 blob — NO radial footprint above the 3×-floor bar in
  its own far zone (a −0.005 GLOBAL Es offset vs the object-free
  floor is inventory, not profile); P-A2 — no tail of either shape to
  fit; P-A3/P-A4 — moot at null. Pad-11 R2 verdict, now measured:
  the old "no steady monopole" survives at 10× sharper resolution,
  and the MECHANISM is identified (door-closed cycles). The far
  field, if this substrate ever has one, must come from an object
  whose books do NOT close at the door — no such object exists yet.

### 4.2 A-F1 — the chord's footprint (primary + replicate + corpse)

Runs: `runs/astro/uud_frozen{,_s314159,_ctl}.log/.fcs`, T=5000,
451 analysis frames each, window [500, 5000]. The seed-20260802
frozen arm reproduces the REGISTRY run's t=5000 diag row EXACTLY with
the QATOM/FCS meters on — the ASTRO apparatus is physics-silent
(determinism cross-check, free). Replicate ret 0.61 (registered band
0.50–0.63), with 20× the primary's standing circulation (circ +10.2
vs +0.52) — the chord's internal current is seed-variable; its
stability is not.

* **The living chord's footprint is CONTACT-LOCAL.** First-medium-bin
  Es depression: −0.059 ± 0.006 (primary, r≈0.75, ~9σ) / −0.0097 ±
  0.0008 (replicate, r≈1.25, 12σ) vs each run's own bulk. By the
  next half-cell bin the medium is bulk-flat: the g1-analog
  footprint EXISTS at chord mass (pad-11 R1's pair could not see it)
  and its screening length is BELOW ONE CELL SPACING.
* **No far field.** π(r) flat from r=1.25 to 8.75 in every arm;
  far-zone aggregate π scatter ±2×10⁻⁵ (the sharpest transport null
  the programme has produced). No 1/r and no exponential component
  above it — there is nothing to fit.
* **The corpse control separates mass from structure:** the dead
  chord (ctl, ret 0.20, husk) shows NO contact dip above 3.5σ — the
  footprint tracks LIVING Em (pressure balance against held
  holdings), not the object's history.
* Books: frozen cond 29.2 / evap 31.1 / net +2.1 over 5000 t.u. —
  balanced at door-gross ~0.012/t.u.; the medium-current bound from
  flat π (~5×10⁻³/t.u., §4.1 method) is NOT even exercised: the
  chord's door turnover is itself of that order. Scored: P-A1 chord
  = contact bin only (registered prediction 1 confirmed for the far
  zone, sharpened by the resolved contact dip); P-A2 = no tail;
  P-A3 = the contact dip is steady (present in both window halves);
  P-A4 = across blob/chord: footprint depth tracks LIVING Em at
  contact; NOTHING tracks throughput — the medium carries no
  current at any mass on the ladder measured so far.

### 4.2b Registered BEFORE the A-F2 read (2026-08-05, ring still running)

The two measured ladder points imply a mechanism: the footprint is
contact-local pressure balance against HELD Em (chord D cell Em 1.45
→ dip −0.059; blob surface cells Em ≈ bath 0.18 → no dip). The ring's
cells hold Em ≈ 1.1–1.2 at the sweet spot. Prediction for A-F2,
written before any ring frame is read: contact-bin Es dip −0.03 to
−0.06 along the hoop (torus-ρ first bin), far zone π-flat at the
meter floor, no tail. If instead a tail resolves at the ring (48
doors sharing one medium shell) the door-closure picture is wrong at
extended objects — that is what the arm decides.

### 4.3 A-L — the first species spectra (P-L1..P-L4)

Analysis: `runs/astro/spec.awk` (per-event zero-parameter
prediction w = base/(1+1.2·x_cell(t)), base 2.9 emission / 1.65
absorption, x from the run's own diag xUDD); numbers in
`runs/astro/spec_*.txt`.

* **P-L1 PASS (the doublet is real, parameter-free, seed-robust).**
  Object emission (DF) residuals vs the closed form: primary
  n=80, mean −1.75×10⁻⁴, rms 0.017; replicate n=72, mean +1.3×10⁻³,
  rms 0.024. Per-cell, per-window line positions match to 0.15–0.5%
  and DRIFT WITH LOAD exactly as predicted (D line 1.658→1.684
  measured vs 1.660→1.687 predicted across the run; replicate
  1.664→1.706 vs 1.661→1.704). Event-count caveat, recorded not
  softened: the D line carries 29–32 events/window (≥30 bar met);
  the U lines carry 1–10/window (below the registered 30 — the U
  voices sit nearly empty and are nearly dark; position information
  is consistent but statistics-poor).
* **The bath is spectrally dark at the D line — the astronomer's
  condition.** The object-free bath's emission band peaks at
  w ≈ 2.35 (measured here; pad-06's 1.90 was a walled-cavity value,
  not an open-glow value — context corrected by measurement). In the
  1.66–1.70 bins the bath emits 1–3 events TOTAL (of 22,339) over a
  whole run; the chord puts 59 there. **A UUD chord in a glowing
  bath is detectable by its D emission line alone, no cell
  attribution needed** (line/continuum ≈ 30–60×). The U line sits ON
  the bath band (camouflaged); the D line is the species marker.
  P-L3 PASS for D (separation ≈ 0.66 from band peak ≈ 15 linewidths);
  U-line separation from the band is < 1 bandwidth — FAIL for U,
  recorded.
* **P-L2 PASS (anti-Stokes invariant ≈ 0.3%, bar ±1%):** both
  families match their matched-x predictions (emission ~0.2%,
  absorption FD means: U cells 1.52–1.54 ↔ x 0.065–0.09, D cell
  0.9715 ↔ x = 0.582, each within 0.3% of 1.65/(1+1.2x)); the
  emit/absorb pitch ratio at matched load is the law constant
  2.9/1.65 = 1.7576 to that precision.
* **The flux machine has a spectral signature (unregistered
  finding).** The U voices are net ABSORBERS (FD 34/31 events vs DF
  6–10) at the near-empty field pitch 1.52–1.54; the D voice is the
  net EMITTER (DF 59 vs FD 44). Cross-attributed with the convtag
  books (cond 29 in ≈ U-side absorption, evap 31 out ≈ D-side
  emission) and the REGISTRY per-link standing current: **the
  chord eats at U, circulates U→D through the bonds, and shines at
  D — metabolism is visible as a two-line spectrum with an
  absorption line above the emission line's load.**
* **P-L4 (death spectrogram):** the dying ctl chord's tracking
  degrades 3× (rms 0.045 vs 0.017 living; mean −0.006, emitted
  atoms redder than the load prediction) — the fade outruns the
  clock's retune; the line residual is a death meter.

### 4.4 A-F2 — the embedded nv=48: the medium kills the flagship

Run: `runs/astro/ring48_L36.log/.fcs` (T=1500, NC=22421, carve
leaves bath outside clear=14.15); floor `bath36` (T=800, π = 1.01414
± 1.1×10⁻⁵ flat — the sharpest floor of the campaign).

* **Registered prediction 5's ALTERNATIVE is what happened.** The
  medium folded the ring AT SETTLE: gyration shape (83.9, 83.9) →
  (12.6, 14.1) already by t=200 while bond lengths stayed near-rung
  (edge_dev mean 0.61) — a dropped necklace, not snapped links. The
  fold destroys the load-bearing straightness (the ~170° angle gate),
  bonds die (gg 0.15 @ t=200 → 0.02), and the object then STARVES:
  books cond 15.4 in / evap 58.9 out / net −36.1 over 1500 t.u. —
  intake collapsed in the sub-nucleation medium while evaporation
  continued. ret 0.038 at t=1500. **The E-A/INTEGRATION stability
  (ret 0.95 @ 50k) is hereby measured as CONDITIONAL on the
  carved-vacuum + noise-drive configuration: "at every ambient" was
  never "in a medium."** (New OUTLOOK §2 limitation.)
* **§4.2b's contact-dip prediction: UNTESTABLE AS REGISTERED** — its
  premise (a living hoop holding Em against the medium) never
  existed; the fold scrambles the contact geometry. Not scored as
  pass or fail; the dip-tracks-living-Em mechanism stands on the
  chord/blob/corpse triple of §4.1–4.2.
* **The death era shows the only resolved medium disturbance of the
  campaign — and it is NOT a space monopole.** Against the run's own
  late-era level, the early window [100,400] carries a box-wide π
  excess +0.0066 (contact) → +0.0024 (r=16), best power-law
  r^(−0.18) — nowhere near a 1/r monopole, and flat again (±1×10⁻⁴)
  by the corpse era. Composition: box-wide Em roughly doubles
  (0.02 → 0.04/cell) — the bleed leaves through the FIELD channel as
  radiation, crosses the box, and recondenses everywhere; the
  released Es of the drained cells shows at contact as an EXCESS,
  not a dip. **A dying object has luminosity, and luminosity leaves
  no space far field** — the transient is dM/dt bookkeeping in
  profile form, the v89 g4 verdict caught in the act.
* **P-A4 closes across the ladder:** the medium shows a gradient
  ONLY where the books do not close (dying ring, net 4× outflow);
  where they close (living chord net ≈ 0, blob net 1.2% of gross)
  the medium is π-flat to its floor. Depth-tracks-throughput is
  REFUTED; depth tracks NET, and net ≈ 0 is what stable matter IS
  here.

## 5. VERDICT (agent's readings; decisions are the user's)

1. **B6 (far field) stays ABSENT — with the sharpest and most
   mechanistic statement the programme has produced.** Around every
   stable body the medium's transport potential is flat at the meter
   floor (chord far zone ±2×10⁻⁵; blob ±2×10⁻⁴; floors ±1×10⁻⁵),
   because a stable object's space cycle closes at its own
   conversion door (backsplash refills cond's pull to 0.3%; the
   Bondi in/out duality cancels pointwise, not through the medium).
   The only medium disturbance ever resolved is the transient
   radiative bookkeeping of a DYING object (r^−0.18, field-channel,
   gone with the corpse). Consequence, stated for the record: in
   this law the far field cannot come from stable matter as it now
   exists — it requires an object whose books run THROUGH the
   medium while it lives, i.e. medium-carried persistent currents —
   the same S2-coherent/identity construction every other absence
   already points at. g2/lensing re-runs stay unmotivated until
   then.
2. **B4 (spectra of species) MOVES.** The campaign delivers the
   programme's first species spectroscopy: a parameter-free doublet
   whose positions track load to 0.02% (mean) / 1.7% (per-event
   rms), seed-robust; an emission line standing 30–60× above a
   spectrally-dark bath (the astronomer's condition — species
   identifiable at range by light alone); the anti-Stokes gap at
   the law constant to 0.3%; metabolism visible as
   eat-at-U/shine-at-D; a death spectrogram (3× tracking loss,
   redward). Proposed audit re-grade: B4 ABSENT → STRUCTURAL with
   the first QUANTITATIVE species row (positions parameter-free; no
   real-atom numeric correspondence claimed). Honest edges: U lines
   statistics-poor (1–10 events/window) and camouflaged on the bath
   band; only the D line is astronomer-detectable.
3. **A NEW measured limitation joins OUTLOOK §2: flagship stability
   is vacuum-conditional.** Real-medium contact kills the nv=48 by
   crush-then-starve (fold at settle → angle-gate death →
   evaporative drain in a cold medium). Small chords survive
   embedded (the frozen UUD lives its full 5000 in the same medium
   class that killed the ring in 600) — at this law, EMBEDDED
   stability exists only below the fold scale. Untested and cheap
   if wanted: a WARM medium (amp 0.5 glow) embedding arm.
4. **Nothing was adopted.** No kernel, law, or table change anywhere
   in the campaign (config-only arms on the standing binary; the
   seed-20260802 frozen arm byte-reproduces the REGISTRY diag row
   with the meters on). Battery re-verified at defaults after the
   campaign: see `runs/BATTERY_astro.log` (must read ALL GREEN 93).
5. **Recommendation (for the user to accept or reject):** the next
   lane is unchanged by this campaign — parcel-carried ontological
   identity (REGISTRY §5.5) — and is now motivated from one more
   side: medium-carried books are the measured prerequisite for a
   far field (this campaign), identity is the measured prerequisite
   for medium-carried coherent structure (REGISTRY), and the
   spectroscopy built here is the instrument that will fingerprint
   whatever that lane builds, the day it lives. Cheap parallel
   probes if wanted: warm-medium embedding; a D-heavy chord for
   U-line statistics; the C7 no-law audit upgrades (DS visibility
   curve, Born KS) remain open from the R-campaign queue.

---

## 6. DECIDE — the decision arms (registered 2026-08-06, BEFORE any §6 run)

User directive (2026-08-06): *"Do initial runs to clarify the
decision. Need more data to decide."* The decisions on the table are
§5's: **(a)** the next-lane choice (parcel-carried ontological
identity, REGISTRY §5.5) and **(b)** the B4 audit re-grade. The arms
below are config-only measurements on the standing binary
(`v91/freecell`, 2026-08-05, battery-green; kernel and laws
untouched). The identity lane itself is **not** opened here — it
needs kernel work and explicit user authorization. Everything in
§6.1–§6.4 was written and committed before the first §6 run.

### 6.1 The confound these arms remove

The two decision-bearing stability results sit at opposite corners of
a 2×2 in (scale, medium temperature), with apparatus varying too:

| body | apparatus | medium | box | result |
|---|---|---|---|---|
| chord nv=3 | frozen cant+tune, reg meter | amp 0.5 (warm) | L=16 | LIVES — ret 0.5101@1500, 0.4975@5000 |
| ring nv=48 | none (V2g+radiance) | amp 0.15 (cool) | L=36 | FOLDS at settle (rms 12.95→5.32 by t=200), starves — ret 0.0376@1500, evap/cond 3.83 |

"Vacuum-conditional" (§5.3) attributes the ring's death to embedding
*as such* — but medium temperature, scale, and apparatus are all
confounded between the two corners. The lane decision needs to know
**which knob kills**: if the cool medium kills everything (even
chords), the frozen-bound reference itself is warm-conditional and
the lane's acceptance target needs a medium spec; if warmth saves the
ring, the limitation softens to cool-conditional and embedded large
matter is viable at this law; if neither — fold and starve persist
warm, chord lives cool — the limitation is scale/identity-shaped and
the lane's urgency is confirmed. §6's E-arms fill the two missing
matrix cells with single-knob differentials, plus one scale rung
(frozen nv=6). The S-arms clear B4's two registered edges (U-line
statistics; a second species). D-A1 measures the identity lane's
design envelope from standing raws — no new run.

### 6.2 Arms (exact commands, standing binary, from `v91/`)

```
# D-E1  ring48 WARM — §4.4 command, ONLY noise_amp 0.15→0.5
./freecell exp=ring ring_n=48 ring_x=0.47 bath=1 noise_amp=0.5 convtag=1 L=36 T=1500 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 qatom_every=200 snap_every=500 snap_file=runs/astro/de1_ring48_warm.fcs > runs/astro/de1_ring48_warm.log 2>&1

# D-E2  chord COOL — the frozen reference command, ONLY noise_amp 0.5→0.15
#       (box stays L=16: the chord's footprint is contact-local (§4.2),
#        so box size is not a medium knob; L=36 would cost ~5h for nothing)
./freecell exp=tri tri_xU=0.28 bath=1 noise_amp=0.15 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 snap_every=500 snap_file=runs/astro/de2_uud_cool.fcs > runs/astro/de2_uud_cool.log 2>&1

# D-E3  frozen nv=6 — REGISTRY i5 geometry, apparatus switched
#       maintained→frozen (cant_tau 50→1e18, reg_gate 2→0, reg_f0 2e-4→0;
#        = the exact frozen-chord apparatus on the i5 ring)
./freecell exp=ring ring_n=6 ring_x=0.47 noise_amp=0.15 bath=1 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 snap_every=500 snap_file=runs/astro/de3_ring6_frozen.fcs > runs/astro/de3_ring6_frozen.log 2>&1

# D-S1  UDD chord — second species, config-reachable via tri_kind=1
#       (smoke T=100 first to verify diag/QATOM formats for kind=1, then:)
./freecell exp=tri tri_kind=1 tri_xU=0.28 bath=1 noise_amp=0.5 convtag=1 T=5000 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 k_cant=1 k_tune=0.2 cant_grow=0 cant_seed=1 cant_tau=1e18 reg_tau=10 qatom_every=1 snap_every=500 snap_file=runs/astro/de4_udd_frozen.fcs > runs/astro/de4_udd_frozen.log 2>&1

# D-S2  two more seeds of the UUD frozen reference (U-line statistics)
#       seed=271828 → de5_uud_s271828.{log,fcs}; seed=141421 → de5_uud_s141421.{log,fcs}
```

UDD cell roles (from the seed code, read before registering): i0 = U
apex, i0+1 and i0+2 = D base; the diag `xUDD` tuple is positional
(per-cell (Em+flload)/cap in cell order), so `spec.awk`'s per-event
law needs no change. The kernel's own seed comment carries a
registered geometry prediction for UDD: the D–D edge d\* = π/w_D
exceeds contact ⇒ **predicted open chain, not triangle**.

### 6.3 Registered predictions, bars, discriminators

**P-D1 (E1, warm ring).** Prediction: the fold is mechanical
(contact/settle), not nutritional — it persists warm: rms ≤ 6 by
t=400 (cool: 5.32@200). Feeding improves the books: evap/cond@1500
< 3.83 and ret@1500 > 0.038. Discriminators for the decision:
  - **No fold** (rms ≥ 10 at t=400) and ret@1500 ≥ 0.3 → the
    limitation was COOL-specific → OUTLOOK "vacuum-conditional"
    softens to "cool-conditional"; embedded large matter is viable
    at this law; identity-lane urgency drops a notch.
  - **Fold, but fed persistence** (ret@1500 ≥ 0.3 as a folded body)
    → EXISTENCE is embeddable, SHAPE is not → the lane's target
    sharpens to shape-holding identity.
  - **Fold and starve** (ret@1500 < 0.15) → the limitation is
    general (not nutritional) → identity-lane urgency confirmed.
  - 0.15 ≤ ret < 0.3 → gray; report as measured, no claim.

**P-D2 (E2, cool chord).** Prediction: the frozen bound is not
warm-conditional — the chord survives the ring's medium class:
ret@1500 ∈ [0.36, 0.66] (±0.15 of the warm 0.5101), ret@5000 ∈
[0.35, 0.65], conn = 1.0 (all three voices bonded) at 5000. Bonus
spectral registration: its doublet must track its OWN xUDD
parameter-free at the new (cooler) operating point — per-event resid
rms ≤ 5%. If instead the chord dies cool → §4.4's scale attribution
weakens (cool medium suspect for everything), the frozen-bound
reference is warm-conditional, and the lane's acceptance target must
carry a medium spec.

**P-D3 (E3, frozen nv=6).** Prediction: the winding wall is not an
identity problem — freezing does not rescue the 6-ring:
t_{ret<0.25} < 1200 (i5 references: ctl 600, maintained 436).
If instead ret@5000 ≥ 0.25 with the chain connected → asserted
identity DOES scale past the chord → the lane's acceptance target
moves up from "chord bound" to "nv6 bound" before it opens.

**P-D4 (S1, UDD).** Geometry: open chain (kernel-comment prediction
on record) — read edge/conn at settle. Spectra: the same per-event
law with NO new constants (resid rms ≤ 5% on ≥ 30 events); line
roles INVERT with multiplicity — UDD = 1 absorber (U) + 2 emitters
(D), so the D line stands ~2× the UUD chord's D line per object and
the U line is weaker still; species distinguishable by light alone
(loci track each species' own xUDD; multiplicity ratio flips). If
UDD won't hold embedded at 5000, that is itself the finding
(species selection — which chords live — recorded, spectra taken on
whatever window it survives).

**P-D5 (S2, pooled U-line statistics).** Pooling 4 seeds (20260802,
314159 standing + 271828, 141421 new) over the settled window:
≥ 30 U-cell DF events pooled, pooled U-line resid |mean| ≤ 0.02 and
rms ≤ 5%. Passing clears §4.3's registered U-row statistics edge →
the B4 re-grade is unblocked on that edge (the camouflage edge — U
loci sit on the bath band — is physics, not statistics, and stays).

**D-A1 (analysis-only, standing raws; descriptive, no bar).** From
`uud_frozen.log.raw` (+ REGISTRY's published per-link ledger): the
identity lane's design envelope, measured before the lane exists —
per-voice door event rate (events/t.u.), gross and net door energy
rates, turnover time τ_turn = Em_voice/(gross door rate) = the
lifetime a parcel gid must survive at the door to be load-bearing,
and the same-cell DF→FD lag distribution (does the door re-eat its
own emissions, on what timescale — the §4.1 door-closure mechanism
seen per-event).

### 6.4 Registered as NOT run, and why

- **Chord-patterned nv=6/12 rings** (the direct scale ladder): not
  config-reachable — `exp=ring` is unison-only; per-cell tune
  patterns on rings are kernel work = the lane itself.
- **Warm-L36 chord** (4th matrix cell at the ring's exact box):
  inferable from E1+E2 unless they contradict; the chord's footprint
  is contact-local (§4.2) so box size is not a medium knob; ~5h/run.
- **C7 audit upgrades** (DS visibility curve, Born KS): not
  decision-critical; stays on the R-queue.

### 6.5 RESULTS (2026-08-06; run in the order registered)

Smoke first: `de4_udd_smoke` (T=100) verified the kind=1 formats —
UDD seeds with wD/wU = 2/3 exactly, the D–D contact check prints
`DD edge OPEN (no contact)` at seed, the diag xUDD tuple is
positional, QATOM fires. NC=2703 at seed 20260802 → i0=2700.

**6.5.1 D-A1 — the door envelope (standing raws; the lane's design
parameters, measured before the lane exists).** Both frozen-chord
raws, settled window [1000, 5000], per voice (v0,v1 = U; v2 = D):

| voice | rate DF/FD (ev/t.u.) | net (energy/t.u.) | Emfl | τ_turn |
|---|---|---|---|---|
| U (seed 20260802) | 0.0020 / 0.0073 | +0.00124 | 0.146 | 46 |
| U′ | 0.0010 / 0.0057 | +0.00118 | 0.145 | 69 |
| D | 0.0118 / 0.0085 | −0.00196 | 1.463 | 278 |
| U (s314159) | 0.0010 / 0.0037 | +0.00065 | 0.142 | 93 |
| U′ | 0.0005 / 0.0080 | +0.00204 | 0.145 | 58 |
| D | 0.0123 / 0.0090 | −0.00219 | 1.398 | 259 |

Seed-robust structure: **U voices are net door-absorbers, the D
voice a net door-emitter** (the flux machine seen per-event at the
door, matching §4.3's eat-at-U/shine-at-D), and the whole-object
turnover time is Σ Emfl / Σ gross = **167–179 t.u.** Lag structure
(same-voice DF→FD vs a memoryless null at the measured rates): the D
door re-eats its own shine at median 29–43 t.u. against a ~59 t.u.
null (~1.4–2× faster than chance), and U absorptions trail the most
recent emission by another object voice at median 34–50 t.u. against
a ~59–85 null (~1.4–1.7×) — mild but consistent across both seeds
and all voices. And a fact worth a line: the chord's door runs at
**bath-cell loudness** (per-voice event rates 0.003–0.020/t.u. vs
bath 0.0095–0.0104/t.u./cell) — what distinguishes living matter at
the door is the standing *direction* of its flow, not its rate.
Envelope for the lane: a door-carried gid is load-bearing if it
survives ~30–500 t.u. (lag quartiles) at a traffic of ~0.04
tags/t.u. per chord — computationally trivial to carry.

**6.5.2 P-D2 (E2, cool chord) — PASS. The frozen bound is not
warm-conditional.** ret@1500 = 0.4258 (bar [0.36, 0.66]), ret@5000
= **0.4114** (bar [0.35, 0.65]), conn = 1.000 with all edges alive
at 5000. The chord survives the ring-killing medium class at ~83% of
its warm retention. The cool medium does redistribute the object:
the U voices run near-empty (final Em 0.0097/0.0084 vs warm ~0.15)
while the D voice holds 1.40 — a starved-lean operating mode, still
bonded, still metabolizing. Bonus spectral registration PASS: the
doublet tracks the cool operating point parameter-free at resid rms
0.34% (DF, n=32) / 0.58% (FD, n=31) — *sharper* than the warm arm.

**6.5.3 P-D3 (E3, frozen nv=6) — as predicted. Freezing does not
rescue the winding wall.** t_{ret<0.25} = **364** (bar < 1200;
references: ctl 600, maintained 436). The full ranking is the
finding: **ctl 600 > maintained 436 > frozen 364** — on a wound
unison ring the lock apparatus doesn't just fail to help, each
escalation of identity assertion dies *earlier*. The winding wall is
geometric, not an identity deficit. Ember state at 5000: ret 0.0326
with the chain still connected (conn 1.0) — same ember class as
i5_ctl (ret 0.0464).

**6.5.4 P-D4 (S1, UDD) — the second species lives; the multiplicity
prediction is REFUTED and the refutation is the discovery.**
Geometry: open chain as predicted by the kernel's own seed comment —
the D–D edge is born beyond contact (`DD edge OPEN`), channel-dead
from t=0.4; the object runs as a bent U–D–D chain. Survival: alive
at 5000, ret 0.3675, with the U apex EMPTY from t≈2000 (x→0.000)
while its bond persists and one D voice drains late — matter
persists, full topology doesn't. Spectra: parameter-free PASS — the
D emission line carries n=61 at resid rms 1.5%, mean −0.003 (≥30
events, ≤5% bar); role structure inverts with the species (U is
almost pure absorber: 31 FD vs 7 DF). Two registered expectations
broke, informatively:
  - **Luminosity is intake-limited, not emitter-limited.** Windowed
    [1000,2000] with both D voices healthy, UDD's D line emits
    0.0170 ev/t.u. against the UUD reference's 0.0230 — *dimmer with
    twice the emitters*, because it has half the intake structure
    (U-absorption 0.0077 vs 0.0130 ev/t.u.). The registered "~2×
    brighter" bar FAILS; what the measurement says instead is that a
    species' brightness is set by its metabolism (intake), not its
    emitter count.
    **[CORRECTION, dated 2026-08-06 (COMBINE §5.2 uniform
    recompute):** the UUD baseline "0.0230" above is NOT
    reproducible — the uniform metric (D-cell DF events/t.u.,
    window [1000,2000], per-run i0) measures the UUD singles at
    **0.0110** (seed 20260802) and **0.0150** (s314159), UDD at
    0.0170 (reproduced), U-absorption UUD 0.0150/0.0170 vs UDD
    0.0090. Corrected statement: UDD is 1.1–1.5× BRIGHTER than a
    single UUD in absolute, but **sub-linear per emitter**
    (0.55–0.77× per D voice) because its intake structure is
    halved — the intake-limitation thesis survives at per-emitter
    grain; the "dimmer in absolute" clause is RETRACTED. The
    species line-ratio discrimination below survives with updated
    numbers: UUD absorption:emission 1.13–1.36 vs UDD 0.53 =
    2.1–2.6× apart. Independent confirmation of the corrected
    baseline: two coexisting UUD chords emit pooled 0.0220 = 2.00×
    the 0.0110 single (COMBINE CO-T2F).]**
  - **The healthy-D locus is universal, not a species tag.** UDD's
    healthy D voice emits at w = 1.6774 ± 0.0071; UUD's at 1.6710 ±
    0.0052 — indistinguishable (Δ = 0.006 ± 0.009), because both D
    voices sit at the same balance attractor. Species identification
    by light alone therefore runs on the **line pattern** — the
    absorption:emission rate ratio differs ~2.3× between species
    (UUD 0.0130/0.0230 vs UDD 0.0077/0.0170 in the same windows) —
    the way real spectroscopy reads abundances from line ratios, not
    from shifted transition frequencies. (The drained D voice's
    locus slides blueward as it empties — w 2.263 ± 0.088 — the
    line *tracks* the voice's load; loci are operating points, not
    species constants.)

**6.5.5 P-D5 (pooled U emission) — n and mean PASS, rms FAILS as
registered, and the decomposition identifies the physics.** Pooled
over 4 seeds: n = 55 (bar ≥ 30 ✓), resid mean = −0.0035 (bar |mean|
≤ 0.02 ✓), resid rms = **0.063 vs the 0.05 bar ✗**. Per-seed: the
one seed whose triangle stayed fully intact to 5000 (20260802) has a
clean U line — rms 0.018, n=21 ✓ — while all three seeds with
degraded topology broaden with redward tails (rms 0.037 / 0.054 /
0.103, outliers to −0.40): §4.3's death-spectrogram signature
operating on *partial* deaths. The U absorption line, by contrast,
is tight in all four seeds (rms 0.010–0.022, pooled n≈198). Stated
conditionally, no bar softening: U emission from intact chords
conforms to the law at ~2%; topology degradation is spectrally
visible in the U line long before the object dies.

**6.5.6 Seed-robustness of the frozen bound, reframed.** Four seeds
at t=5000: ret 0.4975 (conn 1.0, full triangle), 0.6177 (conn 1.0,
full triangle), 0.5430 (conn 0.667, edges opened, all voices
energized), 0.2268 (conn 0.667, one U voice dead). **Matter alive
4/4; full chord topology 2/4.** REGISTRY's "NO DEATH" claim holds
for matter and should be read as scoped: the asserted-identity chord
keeps its *stuff* alive at 5000 in every seed; it keeps its *shape*
in half of them. The topology failures are spectrally visible as
they happen (6.5.5).

**6.5.7 P-D1 (E1, warm ring) — all three registered sub-bars PASS;
the discriminator lands on "fold and starve"; and the late window
adds a measured refinement.** Single knob vs §4.4's A-F2: noise_amp
0.15→0.5. Scored as registered:
  - **Fold is mechanical, not thermal:** rms = 5.52 @ t=400 (bar ≤ 6;
    cool 5.16) — the warm medium folds the ring identically, on the
    same schedule (warm crosses ret<0.25 at t=376, cool at t=460).
  - **Feeding improves the books:** evap/cond@1500 = 165.72/113.20 =
    **1.46** (bar < 3.83) — the warm medium runs both ledger sides
    ~3–7× the cool arm's (cond 113 vs 15, evap 166 vs 59) and closes
    the deficit from 3.83 toward parity.
  - **Feeding improves retention:** ret@1500 = **0.1069** (bar >
    0.038; 2.8× the cool arm), conn 0.438 at close.
  - **Discriminator: ret@1500 = 0.107 < 0.15 → fold and starve → the
    limitation is general (not nutritional) → identity-lane urgency
    confirmed.** Scored exactly as registered.
The refinement the trajectory adds (measured, not registered): the
two arms are indistinguishable to t≈1000 (ret 0.119 warm vs 0.115
cool), then SEPARATE — cool drains monotonically to 0.035 while warm
**plateaus** at ret ≈ 0.12–0.15 for the whole last third (1000→1400:
0.119, 0.140, 0.135, 0.146, 0.123; close 0.107). A warm medium
arrests the terminal drain at a folded remnant of ~¼ the matter — a
FED EMBER steady state — but neither prevents the fold nor recovers
the object. Nutrition sets the remnant's asymptote; it does not
touch the fold or the ~75–89% loss. (Bond census inverts the matter
story: warm holds 25/48 edges alive vs cool 33/48 — the warm medium
churns bonds while retaining more matter; edge count and matter
retention are separate meters.)

---

## 7. DECISION BRIEF (all §6 evidence mapped to §5's two pending decisions; the decisions are the user's)

**(a) Next lane = parcel-carried ontological identity (REGISTRY
§5.5)?** What §6 measured, for and against:

- **FOR — the config-level assertion path is now exhausted by
  measurement.** The three stability arms close the 2×2 confound:
  the frozen chord survives BOTH medium classes (P-D2: cool ret
  0.4114@5000, the only object class that lives embedded anywhere);
  the nv=48 folds and starves in BOTH medium classes (P-D1: warm
  ret@1500 0.107 < 0.15, fold on the same schedule); and MORE
  assertion actively hurts wound topologies (P-D3: frozen 364 <
  maintained 436 < ctl 600). Medium temperature is not the knob;
  scale/topology is. Identity assertion at the config level helps
  exactly one object class (chords) and cannot be escalated — the
  next carrier has to be built into the matter, which is the lane.
- **FOR — the lane's design envelope is measured and cheap (D-A1):**
  a door-carried gid is load-bearing if it survives ~30–500 t.u. at
  ~0.04 tags/t.u. per chord; the door's directed sign pattern
  (U+,U+,D−) is the identity-bearing observable, running at
  bath-loudness — identity by direction, not by rate, exactly the
  REGISTRY §5 thesis at the event grain.
- **SCOPE the acceptance target — matter vs shape (6.5.6):** the
  frozen bound keeps matter alive 4/4 seeds but full topology only
  2/4. The lane's acceptance bars should be split: matter+books
  (the chord meets them everywhere) vs shape (fragile even frozen).
- **NEW — a nutritional floor exists but is not a path (6.5.7):**
  warm feeding buys a fed-ember remnant (~0.11–0.15 plateau), not an
  object. Any future "keep large matter alive by feeding" proposal
  is bounded by this measurement.

**(b) B4 audit re-grade ABSENT → STRUCTURAL?** The two registered
edges from §4.3, as they landed:

- **Statistics edge: cleared conditionally (P-D5).** Pooled U
  emission n=55 ≥ 30 ✓, mean −0.0035 ✓; rms 6.3% misses the 5% bar
  pooled but passes on intact topology (1.8%, seed 20260802) — the
  excess is measured physics (partial-death redward broadening), not
  instrument noise. U absorption is tight in all seeds (n≈198).
- **Species row: STRENGTHENED and corrected (P-D4).** The second
  species lives (UDD, open chain exactly as the kernel's seed
  comment predicted, ret 0.37@5000) — but species identification
  works by **line pattern/ratios** (absorption:emission 2.1–2.6×
  apart between species, numbers per the dated 6.5.4 correction),
  NOT by shifted loci: the healthy-D locus is a
  universal balance attractor (UUD 1.6710±0.0052, UDD 1.6774±0.0071,
  Δ = 0.006±0.009). Brightness is intake-limited at per-emitter
  grain — sub-linear per emitter when intake halves; the "dimmer in
  absolute" clause was a baseline error, RETRACTED in the dated
  6.5.4 correction (UDD is 1.1–1.5× brighter absolute). If the
  re-grade is accepted, the row should read: species by ratios and
  pattern on a universal line locus — which is how real spectroscopy
  reads abundances.

Nothing in §6 was adopted; every arm was config-only on the standing
binary, and the §6.2 registration was committed before the first run.

**[DECIDED 2026-08-06 — user: "first regrade B4, then open
parcel-carried ontological-identity lane." Both executed in that
order: (b) v90/REALITY.md §B4 re-graded ABSENT → STRUCTURAL as
proposed above (the dated row carries the corrected brightness
numbers and the radiance-conditionality edge from COMBINE CO-AR);
(a) the identity lane is OPEN — `v91/IDENTITY.md` is the lane's
registration and record.]**

---

## 8. FINITE-SIZE addendum (registered 2026-08-06, user task #23, BEFORE the runs)

The fold's geometric suspicion, stated: the nv=48 ring (R ≈ 12.9)
sits ~1.2 cells from uncarved bath at seed in the L=36 box (carve
clear 14.15) — the fold could be a BOX/CARVE artifact, not medium
physics. Two single-knob discriminators vs §4.4's A-F2:

* **SZ-1: nv=24 ring, L=36 cool** (double clearance margin, same
  box). No fold (rms stays > 10, ret@1500 ≥ 0.3) → the fold is
  size-relative-to-box; fold-and-starve anyway → the embedding
  limitation is size-independent.
* **SZ-2: nv=48 ring, L=48 cool** (same ring, bigger box). No fold
  → THE FOLD IS A BOX ARTIFACT and the "vacuum-conditional /
  limitation-is-general" verdicts take a dated correction; folds
  anyway → the fold is contact-with-medium physics, box-independent,
  and the standing verdicts hold.

P2 note, honest: the P2 batching economy is a Go-kernel benchmark
experiment (`exp=p2lc`), not a general run mode — wiring it into the
general step loop is future engineering, recorded, not done here.

### 8.1 Results

(registered empty)
