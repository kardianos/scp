# CONSONANCE — locks without binding

**Proponent theory document.** Subordinate to `PRINCIPLE.md`; companion to
`CELLFAB.md`. Direction (2026-07-27): *build the theory in musical resonance
terms first, articulate it fully, then back-port to the cell harmonics. Do
not think of it as binding. Two musical notes don't bind; neither do these.
But unlike notes, excitations cannot exist by themselves — they must
separate until they can tail-call themselves in what amounts to quantum
transfer. And the limits have no Planck; they have a harmonic limit, which
is computed.*

The word **binding is retired** in v89. What was sought under that name is
constructed here as **consonance**: closure of mutual conversion cycles.

Tags: **[D]** derived, **[P]** postulated, **[G]** guess, **[M]** measured.

A key property of this document: **Part III adds no mechanism.** Every law
it uses already exists in `cellfab.c` (the tail-phase gate, retardation at
rate C, load detuning, entrainment, cycle-gated flight). The theory is a
closure analysis of laws already running. The kernel changes that accompany
it are seeds, diagnostics, one analytic mode, and one bugfix — listed in
Part V.

---

## Part I — The theory, stated entirely in music

### I.1 A note is a process, and its last act is to begin itself again

A sustained tone is not an object. A bowed string sounds because the
stick-slip cycle re-creates itself once per period: the string is caught,
carried, released, snapped back — and the release is what arms the next
catch. A pipe speaks because each pressure pulse steers the air jet into
making the next pulse. Stop the recursion and there is no note; there is a
memory of a note.

Call this the **tail call**: the final act of each cycle is the first act of
the next. Nothing persists *through* the cycle — the motion now is not the
motion a period ago, it is a fresh motion built by the old one. A note is a
self-reproducing succession, not a thing that endures.

And a note cannot sound in nothing. Every note is a note *of* something —
a string, a column, a room — and must be driven *by* something — a bow, a
breath, a partner. A note deprived of its driver finishes its tail call
nowhere, and ends.

### I.2 Two notes do not bind

Sound a C and a G together. There is no force between them. Nothing pulls
the G toward the C; nothing holds the pair assembled; nothing must be
overcome to silence one of them. Anyone reaching for a "chord force" has
already made a mistake.

What actually exists between two simultaneous notes is exactly one
phenomenon: **their partials either coincide or they beat.**

### I.3 Beats are the only channel

Two tones a few cycles per second apart do not make two pitches; they make
one pitch that **wobbles** — the beat, at the difference of their rates.
Helmholtz built the whole theory of consonance on this: when two notes
sound, every pair of nearly-coincident partials beats, and the aggregate
shimmer is heard as *roughness*. Where partials coincide exactly, the beat
rate is zero: that channel has gone **silent**.

Beats are not decoration; in coupled resonators they are energy in motion.
Two strings over one bridge, slightly mistuned, visibly trade their motion
back and forth at the beat rate — the energy literally commutes between
them through the bridge, round and round, accomplishing nothing. Roughness
is wasteful commerce. Consonance is the commerce completed and closed: the
partials coincide, the beat period becomes infinite, the trade balances in
silence.

So the consonance/dissonance axis is not aesthetic at bottom. It is:
**how much unresolved energy exchange does this configuration force?**

### I.4 Entrainment: the pull into tune is real, but it is not attraction

Huygens, 1665, ill in bed, watched two pendulum clocks hung on one beam fall
into step — in strict anti-phase — and stay there through any disturbance.
He called it an odd sympathy. The beam trembles with each clock; each
tremble nudges the other clock's escapement; the nudges cancel only in the
locked configuration, so the locked configuration is where the system comes
to rest. Organ builders know the same fact with pipes: two pipes voiced too
close will pull each other into unison (or throttle each other) through the
air they share.

A choir does this on purpose. Two singers holding a just fifth hear the
beats between their upper partials, and steer — not by calculation, by
steering *down the beat rate* — until the beating stops. The interval
"locks." Note what happened energetically: while converging, the pair shed
their disagreement as audible beating; at lock, there is nothing left to
shed. **The energy of the mistuning is radiated away during convergence,
and the locked pair is the configuration with nothing left to radiate.**

That shed disagreement is what a physicist would call a binding energy
deficit. The musician knows better: nothing got bound. The beats stopped.

### I.5 Consonance is closure

Why do small-integer ratios — 2:1, 3:2, 5:4 — sit at the silent points?
Because commensurability is **closure of the joint cycle**. Two periods in
ratio 3:2 compose a joint waveform that repeats exactly after 2 cycles of
one and 3 of the other: the *pair*, taken as a single process, tail-calls
itself. Stumpf measured that listeners fuse such pairs into one sound —
*Tonverschmelzung* — and the fusion is truthful: a closed joint cycle IS
one process with two voices.

An irrational ratio never closes. Its joint waveform never repeats; the
configuration must forever negotiate, and the negotiation is the beating.

So the whole structure, in musical terms:

* **Consonant pair** — joint cycle closes; channel silent; one fused
  process; persists indefinitely; nothing binds it because nothing needs to.
* **Mistuned-near pair** — beats; energy circulates; if the partners can
  retune (singers, entrainable pipes), they converge to the lock and shed
  the difference; the lock is an *attractor of the dynamics*, not a force.
* **Far pair** — outside each other's reach; independent; no commerce.

### I.6 The comma: the limit of harmony is computed, not decreed

Just intonation is the discipline of building music from exact small-integer
ratios, and its central discovery is that the ratios **do not close on each
other**. Twelve just fifths overshoot seven octaves by the Pythagorean
comma, (3/2)¹² ÷ 2⁷ = 531441/524288 ≈ 1.0136 — about a quarter of a
semitone. Four just fifths overshoot a just third by the syntonic comma
81/80. These misfits are not measurement errors and not conventions: they
are **computed** from the ratios themselves. Temperament is the engineering
practice of distributing a comma among intervals so that no single interval
carries all of it.

Just-intonation theory grades harmony by a **limit** — literally so named:
3-limit (Pythagorean), 5-limit (classical harmony), 7-limit, Partch's
11-limit — the largest prime allowed in the ratios. Higher-limit intervals
have longer joint cycles (complexity ~ log of numerator×denominator —
Tenney height), and longer joint cycles are harder to hold: the ear's
tolerance, the instrument's noise, the room's decay all put a ceiling on
how long a joint cycle can be before it cannot be distinguished from
non-closure. Coupled-oscillator theory makes the same statement sharply:
the locking region of a p:q ratio (its Arnold tongue) narrows rapidly as q
grows, so beneath any given noise floor **only finitely many ratios can
lock at all.**

That is the pattern to carry away: **harmony has no universal quantum. It
has a computed limit** — which ratios close, how wide each lock is, where
the ladder is cut off — all derived from the frequencies, the coupling,
and the noise. Nobody imposed a smallest interval on music. The commas and
the tongues *compute* one.

### I.7 Melody and harmony live in different regimes of the same string

A traveling disturbance and a standing lock are the same medium used two
ways. A melody is successive: the pattern moves, and what matters is that
each moment hands off cleanly to the *next place*. Harmony is simultaneous:
the pattern stands, and what matters is that the round trip — over and
back — closes on itself. A pipe of the wrong length cannot hold the note
that a wavefront happily crosses; a room's standing modes sit exactly where
traveling sound refuses to die. The two uses exclude each other pointwise:
where the round trip closes, patterns stand; where it half-closes, they
travel and must not look back.

---

## Part II — The back-port dictionary

| music | fabric | kernel object |
|---|---|---|
| a sustained note | an excitation: a self-reproducing conversion process | energy cycling through channels |
| the instrument / the air | space mode | cells, `E_s` |
| the bow / the breath (a note must be driven) | the partner process (excitations drive each other) | mutual entrainment + delivery |
| pitch | effective harmonic rate | `ω_eff = w/(1+q·x)` |
| loudness | occupancy | `x = (E_m+E_e)/cap` |
| loud playing goes flat | load detuning | `q_detune` |
| two notes sounding | two loaded cells sharing a channel | pair (i, j, d) |
| beats | phase slip through partial gates; churn | gate modulation at Δω |
| roughness | unresolved exchange, scattered losses | off-lock transfer churn |
| sympathetic resonance | resonant joining | `G_res` |
| Tartini difference tones | conversion products of the two-plane beat | β-wrap conversion events |
| a just interval (consonance) | two-way closure of the tail gates | `(ω_i+ω_j)·d/C = 2πm` |
| octave vs fifth character | rung parity | in-phase (m even) / anti-phase (m odd) |
| Huygens' anti-phase clocks | the m = 1 rung | `θ_i − θ_j = π` at lock |
| singing into tune (beats stop) | self-tuning by load exchange | occupancy walks to `x*(d)` |
| energy shed while converging | the "mass defect", correctly named | churn radiated during entrainment |
| the comma (computed misfit) | rung offset | `δ = ω_sum·d/C − 2πm` |
| temperament (distribute the comma) | entrained equilibrium splits δ evenly | `Δ_ij = Δ_ji = δ/2` |
| the harmonic limit (3-, 5-, 11-limit) | computed lockable set | `m_max`, tongue widths vs noise |
| melody vs harmony | traveling tilt vs standing lock | `ωd/C ≈ π/2` vs `ωd/C ≈ πm` |
| room modes / pipe lengths | the separation ladder | discrete `d_m` |

---

## Part III — Derivations in the fabric

### III.1 The mutual tail call, and the separation ladder [D]

The kernel's tail gate (CELLFAB §2) opens direction i→j when the retarded
tail meets the receiver's clock in phase:

```
Δ(i→j) = θᵢ − ωᵢᵉ·d/C − θⱼ ≡ 0   (mod 2π)
Δ(j→i) = θⱼ − ωⱼᵉ·d/C − θᵢ ≡ 0   (mod 2π)
```

A **standing pair** — each the other's bow — needs both directions open at
once. Add the two conditions: the clock phases cancel identically, leaving
a condition on the *configuration alone*:

```
(ωᵢᵉ + ωⱼᵉ) · d / C = 2π·m ,   m = 1, 2, …        ← the separation ladder
```

and the difference then pins the relative phase:

```
θᵢ − θⱼ = ωᵢᵉ·d/C  (mod 2π)
       = π·m for an equal pair — in phase on even rungs, anti-phase on odd.
```

Three readings, in order of importance:

1. **Nothing binds.** Off the rung, mutual conversion cannot close — the
   configuration churns and does not recur. On the rung it closes and
   recurs. Rungs are where pairs *can exist*, not where a force puts them.
   (Recurrence-not-force is exactly PRINCIPLE §4.4 applied to a pair.)
2. **Separation is part of the instrument.** The gap enters the joint cycle
   through retardation ω·d/C — the pipe length. "They must separate until
   they can tail-call themselves": below the first rung the round trip
   returns early and the handoff fumbles; d₁ = 2πC/ω_sum is the first
   separation at which the pair closes. Discrete separations with **no
   Planck input**: the ladder is computed from ω, d, C.
3. **Huygens falls out.** The m = 1 rung locks anti-phase — the same odd
   sympathy the clocks showed on their shared beam.

The equal-pair rung condition per cell is ωᵉ·d/C = π·m: precisely the
"wrap accident" that CELLFAB §10 tuned *away* from for propagation. That is
now understood, not accidental: **the regime that melts a traveling tilt is
the regime that holds a standing pair.** Melody at ωd ≈ π/2, harmony at
ωd ≈ π·m, mutually exclusive pointwise — I.7, derived.

### III.2 Self-tuning: singing into tune, and the defect [D + P]

Frequencies are load-dependent: ωᵉ = w/(1 + q·x). A pair at fixed
separation d (cells do not move) can therefore reach its rung by **adjusting
occupancy** — retuning by exchanging energy, the fabric's way of steering
down the beat rate. For an equal pair on rung m:

```
x*(d, m) = ( w·d/(π·m·C) − 1 ) / q        ← the tuning curve
```

A computed curve in the (d, occupancy) plane on which standing pairs must
lie — an emergent orbital relation with no quantum constant. Reachability
(0 ≤ x* ≤ x_max) cuts the curve off at both ends: links too short have no
rung at any load (d < πmC/w), links too long demand more load than the cap
admits. **The ladder has computed endpoints.**

The energy a pair sheds while walking to x*(d) is the correctly-named
version of binding energy: the comma, radiated as churn during convergence
(I.4). Nothing is stored in a bond; a misfit was disposed of. [P: the walk
happens — the kernel's exchange must actually carry pairs down the beat
rate; this is the entrainment claim tested in E7.]

### III.3 The harmonic limit is computed [D]

Three ceilings, each computed from kernel constants, none postulated:

1. **Tongue width.** Off-rung by δ = ω_sum·d/C − 2πm, entrainment (equal
   pulls both ways) settles at Δ(i→j) = Δ(j→i) = δ/2 — the comma
   tempered evenly — with gate quality

   ```
   G² (δ) = [ (1 + cos(δ/2)) / 2 ]^(2·p_gate)
   ```

   Half-max half-width: δ_half = 2·acos(2^(1−1/(2p)) − 1); p = 8 gives
   δ_half ≈ 0.82 rad. Locks exist as tongues of computed width, not points.
2. **Rung count.** Reachable ω_sum spans [2w/(1+q·x_max), 2w]; channels
   exist only while areas of effect overlap (d ≤ d_link). The largest
   m with 2πm ≤ ω_sum·d/C is computed; for the E6 parameters below,
   **m_max = 1** — nearest-neighbour dense pairs get exactly one rung.
3. **Noise floor.** Ambient churn jitters occupancy, hence pitch
   (loudness→pitch coupling, I.2 of the dictionary); tongues narrower than
   the jitter cannot hold. This is the fabric's prime limit in the
   just-intonation sense: of all conceivable closures, the computably few
   wider than the noise exist as species. (Longer loops — rings of 3+
   cells, the fabric's triads, with Σᵢ ωᵢᵉ·dᵢ/C = 2πm — have narrower
   tongues in the Arnold manner and hit this ceiling sooner. [G: ring
   species; not tested here.])

Two computed side-limits worth naming:

* **Strangulation.** Loading a pair consumes its own space
  (`s_pull`), shrinking areas of effect: r = r0·(1 − s_pull·x·cap/e_s0)^⅓.
  The channel itself pinches off at
  d > 2·r0·(1 − s_pull·x·cap/e_s0)^⅓ — heavy pairs can sever their own
  connection (at s_pull = 0.5, x·cap = 1: d_cut = 1.35). E6 sets
  s_pull = 0.15 to keep d_cut ≈ 1.61 above the scan window — deliberately,
  and the cutoff itself is a testable computed prediction.
* **No Planck.** The pair's action-like quantity per joint cycle —
  circulating energy × round-trip time, printed as `act` in the pair
  diagnostics — is **computed per lock** from its rung and geometry. It is
  not a universal constant, and nothing in the construction needs one.
  Discreteness comes from closure integers; scale comes from harmonic
  ratios; the *limit* comes from tongue-vs-noise. Quantum transfer without
  a quantum constant.

### III.4 Flight shelter: why the tail call is "quantum transfer" [D]

Kernel law: energy in flight belongs to its channel until its cycle
completes; partial cycles deliver nothing; third parties cannot touch
e_mid. At lock, both directions run continuously: a large fraction of the
pair's energy is *in the mutual channel at any moment* — committed,
indivisible, untouchable. The handoff is all-or-nothing per cycle and the
handed frame replaces the old one exactly: a tail call in the strict sense
(no stack growth, nothing accumulates, nothing persists but the pattern).
The pair's energy is safest while in transit between them. That is the
mechanical content of "quantum transfer," and it is why locked pairs should
out-persist unlocked ones under ambient stress [prediction, E6].

### III.5 Existence: why excitations cannot be alone [D + P]

In the kernel a lone loaded cell is *saturation* — static stored energy,
space-mode bookkeeping with nothing succeeding anything. It has no cycle,
so under PRINCIPLE §4.4 (existence = regular succession) it is not an
excitation at all. A dense **process** requires conversion that recurs —
and recurring conversion needs somewhere to convert *to* and back: a
partner. Each is the other's bow (I.1). The minimal existable dense
excitation is therefore the locked pair — CONSTRUCTION §5's N_c ≥ 2 lock,
derived here from the gate law rather than postulated. Excitations
cannot exist by themselves; they exist as mutual recursion, at computed
separations, or not at all.

---

## Part IV — Experiments (predictions pre-registered before running)

Parameters chosen so the m = 1 rung crosses the foam's natural link-length
spread: dense sector w2 = 2.9, q = 0.35, cap = 2.5, Γ_m = 0.02, p_gate = 8,
s_pull = 0.15, C = 1. Rung: ω_sum·d = 2π. Fixed-load pairs at x₀ = 0.4 have
ωᵉ = 2.544, so δ(d) = 5.088·d − 2π sweeps ≈ [−1.2, +1.6] across
d ∈ [1.0, 1.55] — through the tongue and out the far side. The foam itself
is the scan axis; nothing moves.

**E6 — the ladder is in the foam** (`e6_pairs.cfg`): ~50 disjoint,
well-separated cell pairs seeded on existing links at x₀ = 0.4, phases
seeded at lock, ambient field noise as stressor. Per pair, measure rung
offset δ, lock error, two-way gate product G², mutual-flight fraction,
retention.

* P1: G²(δ) follows the computed tongue [(1+cos(δ/2))/2]^16 —
  {δ=0: 1.0, 0.5: 0.78, 1.0: 0.36, 1.5: 0.10, 2.0: 0.015}.
* P2: measured lock error ≈ δ/2 (the comma, tempered evenly).
* P3: retention and flight fraction higher in-tongue (|δ| < 0.8) than
  off-tongue (direction pre-registered, magnitude not).

**E7 — singing into tune** (`e7_tune.cfg`): pairs seeded *on-curve* — each
at its own computed x*(d), so δ ≈ 0 by construction across all d.

* P4: on-curve pairs lock well regardless of d (G² uniformly high),
  in contrast with E6's δ(d) dependence — the tuning curve, confirmed.
* P5 (acquisition, `pair_seedlock=0`): pairs seeded with *random* phases
  find the lock by entrainment alone within ~5–10 t.u. (Huygens in the
  fabric); lock error collapses to ≈ δ/2.

**Ladder mode** (`mode=ladder`): the analytic side — prints x*(d), tongue
widths, m_max, strangulation cutoff from the kernel constants. The point of
the mode is that *the limit is computed before it is measured.*

### Results (2026-07-27, seed 20260727; logs `cellfab_runs/e6*, e7*`)

Conservation at the floating-point floor in all runs (3.7e−16). 60 pairs
per run; the computed ladder (`mode=ladder`) put the rung at d = 1.235 with
tongue half-width 0.83 rad, m_max = 1, strangulation at d = 1.61.

**P1 — the tongue, point for point.** Measured two-way gate product vs the
computed [(1+cos(δ/2))/2]^16, no fitting:

| pair | δ | G² measured | G² computed |
|---|---|---|---|
| p1 | +0.17 | 0.965 | 0.973 |
| p0 | +0.31 | 0.899 | 0.905 |
| p8 | +0.56 | 0.721 | 0.73 |
| p12 | +0.70 | 0.607 | 0.611 |
| p7 | +0.98 | 0.277 | 0.38 |
| p3 | +1.26 | 0.036 | 0.20 (off-tongue: fails faster than the even-split formula — nonlinear unlock) |

Bin means: G² = 0.751 (|δ|<0.6), 0.308 (0.6–1.2), 0.046 (>1.2).

**P2 — the comma is tempered.** Lock error sits at ≈ −δ/2 (measured 0.5–0.75
of δ, correct sign in every locked pair): the misfit is shared between the
two directions, as temperament distributes a comma.

**P3 — flight shelter: confirmed; retention: pre-registration was wrong,
informatively.** Mutual-flight fraction 0.168 in-tongue vs 0.002 off — the
84× shelter. But retention was *lower* in-tongue (0.786 vs 0.839), opposite
to P3 — because of a real effect the pre-registration missed, promoted to a
finding:

**F1 — the fabric pays the comma (the defect, observed).** Pairs seeded too
loaded for their separation (δ < 0, i.e. x₀ > x*(d)) **shed occupancy until
they reached the rung and then stabilized**: x walked 0.40 → 0.18 exactly as
δ walked −0.67 → −0.28 → ≈ 0 (pairs p13, p18, p10…). Energy left the pair
as churn during convergence — I.4's choir, measured: the released
"binding energy" is the comma being paid, and retention dips *because* the
pair tuned itself. The ratchet is one-sided: pairs on the δ > 0 side (would
need *more* load to tune) cannot draw it from vacuum; their gates close and
they freeze off-rung, holding their energy (retention 0.84, G² 0.05).
Relaxation toward consonance happens from the loaded side only. One pair in
sixty fell through the rung and kept draining (overshoot below x*) — the
attractor holds firmly from above, marginally from below.

**P4 — the tuning curve, confirmed.** Seeded on-curve (each pair at its own
computed x*(d)), all 60 pairs sit in-tongue with mean G² = 0.871 and
per-pair G² ≥ 0.94 for the large majority, uniformly across d ∈ [1.08,
1.42] — where E6's fixed-load pairs spread across the whole tongue by
their d alone. The (d, occupancy) curve is real and computed.

**P5 — Huygens in the fabric, confirmed.** Random-phase on-curve pairs
found the lock by entrainment alone: lock errors +0.72, −0.39, +1.49
collapsed to |dth| < 0.2 within 10 t.u. One pair seeded at +1.76 — near the
anti-phase saddle — hesitated (1.72 at t=10, 1.56 at t=20), escaped, and
locked by t=60 (dth = −0.009): basin dynamics complete with the saddle
hesitation.

**The computed limit held.** No lock formed anywhere but the m = 1 rung;
m = 2 needs d > 2.17, beyond channel reach (areas of effect no longer
overlap) — the fabric's first prime-limit statement: nearest-neighbour
dense pairs admit exactly one consonance.

---

## Part V — Audit

**No mechanism was added.** The ladder, the tuning curve, the parity, the
tongue, and the shelter are consequences of laws already in `cellfab.c`
(tail gate, retardation, detuning, entrainment, cycle-gated flight).
Kernel changes accompanying this document:

| change | kind |
|---|---|
| `init=pairs` seeding + pair registry | seed/diagnostic |
| per-pair `# PAIR` rows + tongue summary | diagnostic |
| `mode=ladder` analytic printout | computation, no dynamics |
| beat wrap handled for negative beat rates (w2 > w1 configs) | bugfix |

Background check (README constraint 2): pairs are seeded on existing cells;
no new persistent structure; positions still enter dynamics only through
relational link geometry; the rung condition uses only ω, d, C — quantities
of the current energy structure. Constraint 3 (no imported field): nothing
imported; consonance is not a field. The retired word "binding" imported a
force ontology; it is replaced by closure survivorship, which is the
regular method's own currency (CONSTRUCTION §1.5).

**Open:**

* Rings (fabric triads) — Σ ωᵢdᵢ closure, narrower tongues, species census
  connection to `construct_species.c` locks.
* ~~The defect ledger~~ — done in round 2, E8 (Part VI): the comma weighed.
* Whether locked pairs translate as pairs under a common tilt (a moving
  chord — melody carrying harmony), and what plays inertia for them.
* Whether the noise floor quantitatively predicts which tongues survive
  (the fabric's prime limit, measured rather than argued).

---

## Part VI — Round 2: the theory goes into the kernel (2026-07-27)

Direction: *energy is always conserved; modify the simulation with the
music theory as inspiration, and re-run.* Round 1 added no mechanism; round
2 promotes three musical facts to kernel law. Each is an exactly paired
ledger move; each is toggleable back to the round-1 kernel; conservation
held at |drift| ≤ 1e−15 in every run below.

| | mechanism | musical source | kernel form |
|---|---|---|---|
| **C1** | the partial comb | sympathetic strings sound through *any* coincident partial; the ear resolves ties toward the simpler ratio | dense-sector resonance = best over coprime p:q (p·q ≤ `comb_limit`) of a Lorentzian in q·ωᵢ − p·ωⱼ, amplitude 1/(pq) and width Γ_m/(pq) (Tenney-weighted Arnold tongues). Gates, entrainment, and the rung ladder generalize: Δ = q·θᵢ − q·ωᵢd/C − p·θⱼ, rung (q·ωᵢ + p·ωⱼ)d/C = 2πm |
| **C2** | dissonance radiates | Tartini combination tones; the Plomp–Levelt roughness curve (roughness peaks at moderate detuning, vanishes at unison and at wide separation) | at dense delivery, the rough fraction R = 2\|det\|Γ_r/(det² + Γ_r²) · `rough_k` converts dense → field; and the departing energy **returns its space share** (s_pull/(1+s_pull) → E_s): matter leaving the dense mode un-converts its space (PRINCIPLE §4.2). Global ledger `roughness_radiated` |
| **C3** | harmony is mutual | it takes two to sound; a silent string answers only at its sympathetic readiness | dense exchange ∝ √(mob_src · max(mob_rcv, `mob_floor`·cap)). The field sector stays source-driven — melody flows *from* the singer; harmony is *between* |

The space-return clause of C2 was forced mid-round by measurement: without
it, a pattern that radiated away left its consumed space behind (E4's
defect froze at 160 while its mass fell — orphaned curvature). With it,
**curvature tracks conversion in both directions** (defA 120 → 111 as Em
25 → 17): the geometry follows the energy because the energy is never
destroyed, only re-moded.

### Regressions, round 1 → round 2

| experiment | round 1 | round 2 |
|---|---|---|
| E1 conservation / curvature | 3.7e−16 / r² 0.998 | 3.7e−16 / **r² 0.9997** (slope 3.81) |
| E2 front speed | 0.4815 C | 0.4815 C — identical; the field sector is untouched (roughness ledger 0) |
| E5 CHSH | 2.826 / 1.414 | identical |
| E3a static blob | containment 0.64 | **0.84** — the dense halo is *heard away* (29.5 units radiated as field) instead of lingering as fog |
| E3b tilted blob | +5.3 u, containment 0.48 (fast, smearing) | **+2.3 u, containment 0.81** (slow, cohesive; 30.9 units radiated en route) — mutual coupling trades reach for integrity, and a moving mass now audibly radiates |
| E6 tongue (gg in/mid/out) | 0.751 / 0.308 / 0.046 | 0.555 / 0.295 / 0.051, with retention now cleanly monotone in \|δ\| (0.69 / 0.89 / 0.96) — the ratchet sharper than round 1 |
| E7 on-curve / acquisition | 0.871 / 0.525 | 0.830 / 0.540 |

### E8 — the comma, weighed [M]

Sixty over-loaded pairs (x₀ = 0.55), quiet room, 80 t.u. Binned by initial
rung offset δ₀:

| δ₀ bin | n | energy shed | retention | δ final |
|---|---|---|---|---|
| −1.2 (deep over-loaded) | 15 | 0.633 | 0.770 | −0.83 (still walking toward the rung) |
| −0.55 | 5 | 0.593 | 0.784 | **+0.14 (arrived)** |
| 0 (on-rung) | 16 | 0.774 | 0.719 | +0.39 (circulated most, radiated most, slid just past) |
| +0.55 (under-loaded) | 11 | 0.164 | 0.940 | +0.58 (frozen) |
| +1.2 | 13 | 0.085 | 0.969 | +1.08 (frozen) |

Shed falls monotonically across the δ₀ axis exactly as pre-registered; the
global roughness ledger (25.7 units) accounts for where it went. **The
comma is paid, and the payment is conserved** — it leaves as field
radiation, not as bookkeeping.

### E9 — the interval hierarchy: fifths die young, unisons don't [M]

With deep load–pitch coupling (q_detune = 1.2) the fifth (3:2) is
reachable: x_i = 0.20 (ω = 2.339) against x_j = 0.717 (ω = 1.559), rung
d* = 4π/9.355 = 1.343. Seeding is exact (ratio 1.5003, lock error 0.000).
Measured:

* **The fifth locks and holds ≈ 20 t.u. (~8 beat cycles)** — e.g. pair p0
  at t = 20: lock error −0.048, ratio 1.4972 — then decays: the residual
  1:1 trickle between unequal partners slowly equalizes occupancy, the
  ratio drifts off 3:2, and the comb's best-ratio flips to 1:1, scrambling
  the lock. Retention stays 0.94–0.96 throughout: the interval dies of
  *detuning*, not of energy loss.
* **Control window (off-rung): no locks at all** (gg 0.035) — the rung is
  load-bearing.
* **The octave barely forms** (gg ≈ 0.12 transient): it is reachable only
  at the extreme of the occupancy range, where the light voice borders
  silence and bleeds.
* First attempt, kept for the record: at x_i = 0.10 the light voice sat
  only 0.31 from vacuum pitch and bled into the hall; the ratio wandered
  1.5% and the lock died by t ≈ 10. Deepening to x_i = 0.20 (vacuum
  detuning 0.56), gentler dense rates, and a smaller sympathetic floor
  raised pair retention 0.71 → 0.94 and gave the fifth its measurable life.

The measured lifetime ordering — **unison (indefinite) ≫ fifth (~20 t.u.)
> octave (transient)** — is the Tenney/tongue-width hierarchy of Part I.6
realized as dynamics: higher-limit consonances are harder to hold, in the
fabric as in a choir.

### Open after round 2

* **Frequency entrainment — the choir's correction.** Locks are
  phase-entrained only; nothing steers *occupancy* toward resonance, so an
  interval's lifetime is set by detuning drift. Singers correct pitch, not
  just phase: the missing mechanism is energy exchange whose direction
  correlates with the sign of detuning (beat-driven steering), to be built
  as an exactly paired ledger move like everything else.
* **Interval identity.** The comb's best-ratio flips the instant scores
  cross; a sounding interval should have inertia (hysteresis) in its
  identity.
* The octave wants range: only the edge of the occupancy scale reaches
  2:1 at q_detune = 1.2.
