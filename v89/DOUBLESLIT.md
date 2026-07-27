# DOUBLESLIT — interference and measurement in the cell fabric

**Proponent postulate + instrument construction.** Subordinate to
`PRINCIPLE.md`; companion to `CELLFAB.md` (the kernel) and `CONSONANCE.md`
(joining as consonance). Direction (2026-07-27): *postulate the double-slit
experiment theoretically, then simulate it, including detectors.*

Tags: **[D]** derived from existing kernel law, **[P]** postulated,
**[M]** measured below.

The important structural claim up front: **nothing new is postulated about
interference or measurement back-action.** Both are consequences of
mechanisms already in the kernel — the tail-phase gate and clock
entrainment. What is postulated is only the *instruments* (a wall, a
screen, a which-path recorder), each built as an energy structure, each an
exactly paired ledger move. Energy is conserved through every event of the
experiment, including detection.

---

## 1. The experiment, restated without a traveling particle

Nothing traverses the apparatus (PRINCIPLE §4.4). A source region converts;
the conversion advances as a phase-tilted re-alignment of cell clocks (the
field-sector wave of CELLFAB §2). A **wall** interrupts it: a region whose
field plane runs at a strongly shifted rate — a *spectrally different
medium*. Resonant joining fails against it (G_res ≈ Γ²/Δω² ≈ 10⁻³), so the
wall neither accepts nor transmits: it is a **mute**. Two **slits** are
windows where the fabric keeps vacuum pitch — the only strings left
sounding. Downstream of the wall, every cell is reachable by conversion
paths through slit A and through slit B, and by nothing else.

So the double slit, in this ontology: **one advancing conversion, offered
two families of routes, rejoining in a shared region.** The question "which
slit did it go through" presumes a traveler; there is none. There are two
route families of one process — exactly the structure of the Bell pair in
CONSONANCE §6 (one process, two ends), here as one process, two paths.

## 2. Interference is joining consonance of two retarded tails [D]

A downstream cell at distances r_A, r_B from the slits is brushed by two
tails carrying retarded phases θ_A − ω·r_A/C and θ_B − ω·r_B/C. The slits
are driven by one incident wave, so θ_A, θ_B hold a fixed relation. The
cell's clock entrains to the flux-weighted mean, and the tail gates then
score each arrival against that clock:

* **Δr = r_A − r_B = m·λ** (λ = 2πC/ωᵉ): the tails agree; the clock can
  satisfy both; gates open; energy joins and flows onward. A bright fringe
  is a locus of **consonance of the two routes**.
* **Δr = (m+½)·λ**: the tails quarrel; whatever the clock does, both gates
  score poorly (each ~cos²ᵖ(π/4) at the compromise); joining is refused.
  A dark fringe is not cancellation of stuff — nothing is destroyed — it is
  **refusal of joining**, and the energy takes the consonant routes
  instead. Conservation does not permit "destructive interference" as
  annihilation; it *requires* redistribution.

Musically: the fabric behind the wall is an Aeolian harp sounded by two
voices. Strings placed where the voices arrive in tune ring; strings where
they arrive opposed cannot ring, and the air's energy goes to the strings
that can. The fringe pattern is the room's consonance map of the two
sources.

The wavelength is not a free parameter: **λ = 2πC/ωᵉ** from the kernel's
own clock rate and conversion rate. Fringe loci therefore follow the exact
two-source relation cos²(ωᵉ·(r_A−r_B)/2C) with everything known in
advance — a parameter-free spatial prediction.

## 3. A screen is a condensation medium; detection is conversion [P]

A detector is not an observer; it is a region of fabric that **converts
arriving field into locked record** — dense-mode grains that stay
(CONSTRUCTION §5: a record is a lock; conversion proceeds through complete
cycles, and a partial cycle records nothing). The screen here: cells that
condense incoming field energy at a fixed rate and hold it (deep record
medium: no re-emission, no capacity throttle, grains exempt from
evaporation). Exposure = accumulated dense energy per cell — the
photograph. Grain completions ("clicks") are logged each time a cell's
record crosses a whole grain quantum: **discrete events whose spatial rate
follows the joining flux**, i.e. the gated intensity — a mechanical,
Born-flavored statement [P: rate ∝ flux is the kernel's absorb law; the
claim that nothing more is needed at the single-quantum level is theory,
not yet simulation].

## 4. Which-path recording destroys fringes by the mechanism that carries the wave [D]

To record passage **at a slit** is to convert there: a which-path detector
taps a fraction of the transiting field into a local dense record. But in
this kernel every delivery **entrains the receiver's clock, and every
re-emission carries the sender's clock** — recording and back-action are
not two effects, they are **one event**. A recorder's clock is its own
(thermal, uncorrelated with the source: modeled as a strong random walk).
So the wavelet that continues past a monitored slit leaves re-timed to the
detector's pitch: *the detector cannot listen without humming its own pitch
into the line.*

Consequence, derived: the two routes lose their fixed phase relation; the
downstream consonance map shifts randomly in time; the accumulated exposure
is the time-average over a wandering relative phase — **fringes dissolve
into the sum of two single-slit humps.** One monitored slit suffices (the
relative phase is randomized by either voice going off-book).
Complementarity is here a theorem of the imprint mechanism: you cannot
entrain a record out of a passage and leave the passage's phase entrained
to its source, because both are the same entrainment.

Partial coupling (weak tap, slow detector walk) randomizes the relative
phase only partially per transit → intermediate visibility. [P → D4]

## 5. Instrument audit (constraints 2 and 3)

| instrument | construction | audit |
|---|---|---|
| wall | field plane detuned by `wall_detune` (seeded at init, like any energy structure) | a property of the region's energy structure, not forbidden connectivity; transfer laws untouched |
| slits | absence of wall (vacuum pitch windows) | geometry of the seed |
| screen / edge sinks | absorb: paired move Ee→Em per step; head=1 (deep medium); grains locked (no evap, no dense transfer) | conversion, conserved; idealized record medium [P] |
| which-path recorder | tap: paired move Ee→Em (the record) + clock random walk (its own pitch); transmits otherwise | back-action = the kernel's own entrainment; no new coupling invented |

No imported field; no background: every instrument is cells with seeded
energy structure. Conservation is reported end-to-end including all
instrument conversions.

## 6. Pre-registered predictions

* **P1 (fringes):** two slits open → screen exposure shows maxima/minima at
  the parameter-free loci cos²(ωᵉΔr/2C); visibility against
  predicted-position bins well above the which-path level. (Tempered
  honestly: free-flight clock dephasing ~30 t.u. is comparable to the
  slit→screen transit, so V < 1 is expected; the *positions* are the sharp
  claim.)
* **P2 (control):** one slit blocked → single hump centered off-axis toward
  the open slit; V ≈ 0 on the fringe template.
* **P3 (which-path):** recorders at both slits → fringes collapse to two
  humps; V collapses; per-slit records tally ≈ symmetric.
* **P3b:** recorder at slit A only → V collapses likewise.
* **P4 (partial):** weak recorder → intermediate V between P1 and P3.
* **P5 (the law):** total energy conserved at the floating-point floor in
  every run, screen and detector records included; every joule on the
  screen is a joule that left the source slab.

## 7. Results (2026-07-27, seed 20260727; logs `cellfab_runs/d*.log`)

Eight runs; conservation between 0 and 1.9e−16 in every one, screen and
recorder ledgers included (**P5 confirmed** — every joule on the screen and
in the meters left the source slab).

### What the instruments did correctly [M]

* **Opacity and transmission.** The detuned wall passes nothing measurable;
  the slits pass beams; the screen condenses arriving field into held
  grains; the shutter freezes the record after the first transit (needed:
  the reverberant pool in front of the wall re-feeds the slits incoherently
  — unshuttered visibility-against-template fell from 0.46 to 0.18).
* **Additivity, exactly.** Fitting the two-slit exposure to the *measured*
  single-slit envelopes gives I_AB = 1.01·I_A + 1.00·I_B. The apparatus is
  linear in its sources with no cross-talk — a sharp conservation
  statement, and the correct classical-transport baseline.
* **Which-path recorders record.** Both-slit meters tallied A = 21.1,
  B = 21.3 (symmetric, as pre-registered in P3's record clause); the
  one-slit meter run measured its own insertion loss cleanly (envelope fit
  a = 0.15, b = 0.85 — the metered slit transmits 15%). Strong meters
  attenuate what they dephase — measurement costs transmission, paid in
  energy, conserved.

### The load-bearing negative: no fringes [M — P1 refuted]

After envelope normalization, the two-slit record shows **no interference
modulation**: ratio r(y) = I_AB/(fit) is flat at 1.00 ± 0.1 across the
predicted maxima and minima (V_norm = −0.008). The instantaneous arriving
flux just before the screen is likewise pure envelope — the modulation
never forms even transiently. The earlier fixed-template "V = 0.46" was
envelope shape leaking through the template, and is retracted.

**Why, mechanically.** Two structural facts of the present kernel:

1. **A cell holds one clock per sector, not an amplitude.** Where two
   wavelets overlap, entrainment settles each cell to one compromise
   phase. Anti-phased arrivals are then partially *refused* — but refusal
   is not cancellation. A driven oscillator with signed superposition
   responds weakly to opposed drives; a gated energy-transport cell merely
   queues the energy it declined.
2. **Conservation turns refusal into delay.** Refused energy is not
   destroyed (the law); it waits upstream and passes when foam disorder
   drifts the gates. Integrated over the exposure, every route delivers
   what was headed down it: the record relaxes to the envelope sum.
   A dark fringe requires energy to be *redirected at amplitude level*,
   not declined at transport level.

So the double slit **falsifies the single-phase field sector**: gated
positive-energy transport with one phase per cell reproduces geometric
optics, opacity, records, additivity, and complementarity's cost side —
but cannot produce two-path interference, even in principle, because it
has no signed superposition to cancel with. This is the first experiment
the fabric fails, and the failure is exact (additivity to 1%).

### What the failure points at [P — next construction]

The repair is already latent in the cell's own hardware: each cell has
**two planes**. Let the field sector carry a genuine two-component
amplitude — the pair (a₁, a₂) on the two planes, deposits adding *signed*
in the plane pair, energy being the squared norm (conserved as ever;
transport moves the norm, phases live in the component ratio). That is
precisely the chiral pair u± = (u₁ ∓ i·u₂)/√2 of CONSTRUCTION §2.2, built
from existing structure rather than imported. With signed two-component
superposition, anti-phased arrivals cancel *in the components* while
energy redistributes to the constructive loci — interference as
consonance geometry, conservation intact. Bell (CONSONANCE §6) already
pointed the same direction: the pair-level physics wanted a joint
two-plane amplitude; the double slit now demands it at the single-wave
level. *(The repair was authorized and made — §8.)*

### Scorecard against §6

| prediction | outcome |
|---|---|
| P1 fringes at computed loci | **refuted** — additivity instead (the finding) |
| P2 single-slit envelope | confirmed (and became the analysis tool) |
| P3 which-path records + hump sum | records confirmed; "V collapse" moot — the baseline itself has V ≈ 0 at record level |
| P3b one meter suffices | moot for V; insertion loss measured |
| P4 partial coupling | moot for V; weak meters transmit (a = 0.60, b = 0.56) and are re-entrained by the beam — a pointer dragged along by what it measures [M, worth keeping] |
| P5 conservation through detection | **confirmed at the floor, all runs** |

*(The scorecard above is round 1 — the transport field sector. The repair
changes P1's verdict below.)*

---

## 8. The repair, made — and fringes (2026-07-27, same day)

Authorized: *make the repair and try again.*

### 8.1 The repaired field sector

The field sector is now a **two-component signed amplitude** (a₁, a₂) per
cell — the chiral plane pair, ψ = a₁ + i·a₂ — evolved by exact unitary
steps: an onsite rotation at ωᵉ (load and wall detuning enter here), then
**exact pairwise hop rotations** over live channels; each 2×2 rotation
conserves the pair norm to roundoff, so Σ|ψ|² is conserved without a
ledger clause. Field energy is E = |ψ|², non-negative by construction;
every exchange with the other modes (screen condensation, recorder taps,
beat conversion, roughness injection) scales |ψ| by the exact energy moved
— paired ledger moves as ever. Superposition is native: anti-phased
arrivals cancel in components while the norm redistributes to constructive
loci. The dense sector — every consonance result — is untouched: harmony
stays gated transport, melody became amplitude mechanics.

Instruments under the repair: the wall is still a detuned mute (now via
onsite frequency — off-band, evanescent, reflective); the screen still
condenses |ψ|² into grains; the which-path recorder becomes **unitary
dephasing** (a random phase kick — back-action without imposed
absorption) plus a small record tap.

### 8.2 Three walls hit on the way, kept for the record [M]

1. **Sign–seed pairing.** The seed convention θ = −k·x (inherited from the
   retarded transport era) means the packet carries κ = −k; the hop sign
   must pair with it or everything propagates backward. The front-radius
   diagnostic is direction-blind and "confirmed" the wrong sign; the
   energy-weighted centroid settled it. (Instrument lesson: never measure
   propagation with a radial metric.)
2. **Sub-wavelength channels cut off.** With a broad hopping band, a
   two-cell-wide slit walled by detuned cells is below cutoff and
   transmits ~0.2% by tunneling. Wider apertures (≈ 0.7λ) and a thin wall
   restored %-level transmission.
3. **The foam speckles waves.** Random contact areas and coordination are
   strong diagonal disorder for a linear wave: measured patch coherence
   0.02–0.07 at the slit exits — dead before arrival (and the packet
   centroid decelerated: diffusion, not flight). The framework-honest cure:
   **a cell's total joining bandwidth is a property of the cell, not of
   its accidental contact geometry** — symmetric normalization
   ŵᵢⱼ = wᵢⱼ/√(sᵢsⱼ) of the hop weights. Coherence rose to 0.76 and the
   fringes appeared. [P, load-bearing: vacuum is spectrally uniform even
   though the foam is geometrically irregular.]

Linear regime for the optics runs: q_detune = 0 (the weak-field limit made
exact), else the standing wave before the wall self-phase-modulates.

### 8.3 Results — the experiment, repaired [M]

Seven runs; conservation 5.6e−16 … 1.1e−15 through all detection.
λ = 2π/k = 5.71 seeded; screen exposure 61 units (223 click grains);
envelope fit I_AB against measured single-slit envelopes: a = 0.99,
b = 0.98.

**P1, now confirmed.** The raw record shows the two-slit structure at the
parameter-free loci: central maximum 4.04 at y = 15.4 (pred max), minima
0.5–1.1 at y ≈ 10.9 and 19.9 (pred minima), secondary maxima at y ≈ 4.9
and 25.9 (pred maxima). Envelope-normalized ratio r(y): 1.61–1.68 across
the central maximum, 0.59–0.73 in both predicted minima, recovering to
0.9–1.0 at both secondary maxima. **V_norm = 0.316** (fixed-template
V_central = 0.466). Visibility is ceilinged by residual off-diagonal foam
scattering and source bandwidth, not by the mechanism.

**Complementarity, measured as a dial** (envelope-normalized V at
predicted loci; recorder kick strength σ, tap 0.02):

| run | σ | V_norm | screen exposure |
|---|---|---|---|
| D1 no meters | 0 | **0.316** | 61.0 |
| D4b weak | 0.3 | 0.312 | 43.8 |
| D4 moderate | 0.8 | 0.282 | 17.0 |
| D3b one slit, strong | 2.5 (A only) | **0.025** | 30.8 |
| D3 both, strong | 2.5 | starved (0.9 exposure) — strong dephasing also back-scatters | — |
| D2 single-slit controls | — | ≈ 0 by construction | ~31 |

The punchline row is D3b: fully metering **one** slit collapses the
fringes to the single-slit floor (V 0.316 → 0.025) while the apparatus
still transmits half its light — and the envelope fit shows the metered
slit's *coherent* contribution gone (a = 0.02, b = 0.99). Knowing one
path suffices. Weak meters (σ = 0.3) record (A ≈ 23, B ≈ 24 tallied) and
barely disturb (V 0.312): the weak-measurement corner, mechanically.

**Regressions under the repair:** E5 Bell identical (2.826/1.414); dense
sector untouched — E3a containment 0.880 (best yet), E6 tongue
0.720/0.292/0.046; E1 conserves at 1.7e−15 with conversions exchanging
energy with ψ exactly.

### 8.4 Scorecard, final

| prediction | round 1 (transport field) | round 2 (amplitude field) |
|---|---|---|
| P1 fringes at computed loci | refuted — additivity | **confirmed** — V_norm 0.316, loci parameter-free |
| P2 single-slit envelope | confirmed | confirmed (a = 0.99, b = 0.98) |
| P3/P3b which-path kills fringes | moot | **confirmed** — one strong meter: V → 0.025 |
| P4 partial coupling → partial V | moot | confirmed in direction (0.316/0.312/0.282/0.025) |
| P5 conservation through detection | confirmed | confirmed (≤1.1e−15) |

The double slit first falsified the single-phase field sector, then
confirmed the two-plane amplitude repair the falsification demanded. The
fabric now interferes, records, and pays for what it knows — and the books
balance to the last digit.

---

## 9. The reality ladder, climbed (tiers 1–3; see ROADMAP.md §5 for full tables)

Same-day continuation under the standing criterion: *correspond to
reality; name every deviation.* Kernel additions: single-quantum claim
sampling (`n_quanta`), a second chirality component with a **lossless
unitary which-path tagger** at slit A, and analyzer-basis screen ledgers
with a mid-flight basis switch. Conservation ≤1.7e−15 throughout.

* **Tier 1 (Tonomura):** 500 quanta, one indivisible claim each (R3
  atomicity); dots at scattered (y,z); the click histogram rebuilds the
  fringes (V_discrete = 0.412 vs wave 0.466, shot-noise limited).
* **Tier 2 (duality + eraser):** weak tags sit on the
  Englert–Greenberger–Yasin *equality* (V/V₀ = 0.927 measured vs
  √(1−D²) = 0.929); strong tags go mixed and fall below the bound —
  the pure/mixed marker structure of the real relation. Eraser at ±45°:
  fringes +0.403 in one ledger, anti-fringes −0.292 in the other, flat
  total, exact ledger swap under basis reflection.
* **Tier 3 (delayed choice):** basis chosen at t = 40, after slit
  transit — erasure appears in exactly the post-choice proportion; the
  total distribution is bit-identical (max bin delta 0 — no-signaling
  exact at the ensemble level, by construction: the analyzer never feeds
  back).
* **Forks resolved:** ħ is *not* emergent (19× per-cycle action spread
  among identical locks — the kernel's transfers are rate-proportional;
  the route is CONSTRUCTION R1's deposit quantum δ, unimplemented).
  No-signaling holds exactly.

Remaining rungs: tier 4 (Hong–Ou–Mandel via pairwise joint-amplitude
registries) and tier 5 (fidelity: V→1, dispersion decision, per-quantum
collapse feedback). The double slit now matches reality at every level a
single-quantum wave theory can reach — one click at a time, erasable
markings, choice-time invariance — with its two honest gaps (ħ, HOM)
named and routed.
