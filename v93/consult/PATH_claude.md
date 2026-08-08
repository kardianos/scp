# v93 research path — claude consultation (post faces A/B/C)

Consultant: Claude (Fable 5). Date: 2026-08-08. Brief: `PATH_REQUEST.md`.
Grounding: the brief §1–4, plus the v93 record (`RESUME.md`, `RING_DNLS.md`
outcomes as summarized in the handoff). Everything below is labeled
**[MEASURED]** (established in the record), **[INFERENCE]** (follows from
measured facts under stated assumptions), or **[SPECULATIVE]**.

---

## 0. The honest starting position

What Face C actually delivered [MEASURED]:

- Spontaneous condensation of a spread packet into long-lived dense hoards
  under the unitary channel + load-detuning (w2e = w2/(1+q_detune·x)), with
  NO lock, gate, or door. Stable to t=300. The additive law cannot do this.
- A long-tailed hierarchy (189 hoards, 0.5→7.2), with radiance TRUNCATING
  the spectrum at the fixed point Em*≈1.55 and deep self-traps resisting.
- Deep corner (qd=3.6, amp=2): Em_max≈7 (3×cap), frozen envelope.

What it did NOT deliver [MEASURED]:

- Tightness: PR≈250 cells — a condensed packet, not a soliton.
- Phase: it binds ENERGY, not PHASE (retention probes dead at every corner).
- Mobility: the corner that binds does not translate (kx arms stall).
- Endurance: t=300 is short against the programme's own standards (v91
  chords held t=5000; the QUENCH cloud half-life ~20k t.u.).

And the standing structural facts that constrain the whole path:

- Transport / binding / phase-retention are **pairwise antagonists** in the
  current channel [MEASURED — the RING_DNLS synthesis].
- The medium kills phase structure by contact dephasing in ~10 t.u.
  regardless of law, door handling, or hop schedule [MEASURED — faces A/B].
- Stable matter's space cycle closes at its own door: no space-carried far
  field from stable objects; only dying matter disturbs the medium
  [MEASURED — v91 ASTRO pad-11/g4 closure].
- There is currently **no charge** in the free-cell substrate. Pre-v89 U(1)
  charge lived on a fixed lattice and is banned as a template. Any
  charge-analog must be found, not imported.

The last two points are load-bearing for goals 3–5 and shape the whole
long-term arc (§3).

---

## 1. Short-term: the particulate reckoning (8 experiments + 1 debt)

Ordered. Each entry: apparatus / bar / what the outcome decides. All are
config-only (existing knobs: amp_nat, amp_logate, hop_order, q_detune via
the law table, seed_mw, k_rad, noise_amp, p1_meter) — **no kernel changes
needed for any of P1–P8**. Analysis instruments exist (analyze_melt.py hoard
census, analyze_winding.py / R2d, spec.awk from v91 ASTRO).

### P1 — Endurance: is t=300 a plateau or a slope?
**Apparatus:** rerun the Face C deep corner and one mid corner (qd=0.6–1.2)
to T=3000, then the survivor to T=10000. Snapshot cadence enough for a hoard
census time series (count, Em_max, mass spectrum, per-hoard centroid).
**Bar:** hoard count and spectrum statistically stationary over the last
half of the run; Em_max envelope drift < a few % per 1000 t.u.; individual
hoard identity traceable (centroid continuity) for ≥1000 t.u.
**Decides:** transient condensate vs meta-stable particulate. This is the
cheapest experiment with the highest information; everything else is
downstream. If hoards coarsen (merge into ever-fewer, ever-bigger — the
classical condensate fate) or evaporate, Face C is a *relaxation* phenomenon
and the "particulate" framing dies here honestly.

### P2 — Isolation: does ONE hoard exist on its own?
**Apparatus:** seed a single Gaussian Em packet at a measured hoard's mass
and width (take a census hoard from P1 as the template — self-seeding, no
imported profile), in vacuum (bath=0, k_rad=0, wf_on=0), deep-corner law.
Sweep seeded mass across the measured hierarchy (0.5 → 7).
**Bar:** stationary envelope (rms flat, PR flat) for ≥500 t.u. at some mass
range; a threshold mass below which it disperses.
**Decides:** whether the hoard is a genuine self-trapped solution of the
DNLS-like law (a soliton branch with an existence window — the analog of the
old Q-ball ω-window, but emergent) or a collective object that needs the
hierarchy/bath around it. A measured **existence window in mass** is the
first quantization-adjacent fact of the lane. [INFERENCE from DNLS
phenomenology: discrete self-trapping usually has a threshold norm; whether
this law's version does is exactly what this measures.]

### P3 — Monodispersity: does radiance make a mass UNIT?
**Apparatus:** the Face C sweep extended: (q_detune × amp × k_rad × p_rad)
grid, 3 seeds per corner. Meter: the hoard mass spectrum (histogram, peak
width, seed-robustness of the peak location).
**Bar:** a corner where the surviving spectrum collapses to a narrow,
seed-robust peak (dispersion, say, <10% of the peak mass) — or a clean
discrete ladder of peaks.
**Decides:** whether the fabric SELECTS a preferred mass. This is the single
most valuable measurement on the board: a reproducible preferred mass is (a)
the strongest particulate claim, (b) the in-fabric mass unit that the scale
correspondence (§2) needs, and (c) the first place a *discrete spectrum of
stable configs* (goal 5's shape) could show up. Radiance-as-mass-filter is
already half-measured [MEASURED: truncation at Em*≈1.55]; this asks whether
the filter has a passband, not just a cutoff.

### P4 — The honest medium: does the hoard survive the adopted law?
**Apparatus:** embed a P2-validated hoard (or run Face C to condensation,
then switch on) the full V3a surface: bath=1, wf_on=1, k_rad=0.05, door
live. Measure t_half of hoard mass and count vs the vacuum control.
**Bar:** t_half ≥ 10× the radiance-taxed hoard baseline of the additive era
(v91 RADIANCE measured t_half 80–140 for v90 hoards). Stretch bar: alive at
t=5000 (the v91 chord standard).
**Decides:** whether Face C matter is matter *under the law we adopted*, or
a vacuum-only curiosity. Note the asymmetry with faces A/B [INFERENCE]: the
medium kills **phase** in ~10 t.u., but the hoard binds **energy** — contact
dephasing should be much less lethal to it. If that inference holds, it is
itself a discovery: the first object whose stability does not route through
phase coherence, i.e. genuinely *incoherent* matter. If the medium kills the
hoard too, the antagonism table gains a row and the lane narrows sharply.

### P5 — Internal structure audit: what IS a hoard, inside?
**Apparatus:** per-hoard, over time: (a) th2 phase map + R2d about the hoard
centroid + m-spectrum (does it hold ANY internal winding?); (b) the w2e /
clock spectrum via spec.awk (does it have a LINE?); (c) radial Em profile
(core+skirt? flat-top at cap? self-similar across the hierarchy?); (d) its
Es footprint (does it dig the graded depression that v89 g1 measured for
additive-era mass?).
**Bar:** descriptive, but with two sharp sub-bars: linewidth vs the
spectrally-dark bath background (the v91 D-line was 30–60× over dark — is
the hoard bright or dark?), and Es-depression depth vs a control void.
**Decides:** three downstream forks at once. A clean internal line gives an
ω to hang ε=A₀ω/2π on → connects hoard mass to pitch (the E=Qω echo) and
feeds the scale gauge (§2). Internal m-mode capacity (even m=0-only) is the
fact that gates the substructure stage (§3, Stage IV). The Es footprint
gates the shell stage (§3, Stage III): shells need a potential well, and the
hoard's own space depression is the only lawful candidate.

### P6 — Two hoards: is there ANY interaction?
**Apparatus:** two P2 hoards at separations d = contact, 2, 4, 8 cells;
relative clock phase Δth2 ∈ {0, π/2, π} (the seeds control th2); also
relative pitch (slightly different masses → different loaded pitches).
Measure: separation vs t (approach/recede/static), mass exchange, merge
outcomes.
**Bar:** any reproducible, seed-robust separation drift distinguishable from
the ~1e-3 COE floor; any phase-dependence of the contact outcome.
**Decides:** whether nuclei are even conceivable in this lane. The programme
has NEVER measured a mutual force between stable objects [MEASURED: v91
COMBINE tri2 — 2× luminosity, NO mutual force; ASTRO — footprints
contact-local]. Expectation [INFERENCE]: contact-range interaction only
(merge or exchange), nothing at range — because nothing carries the books at
range (the g4 closure). A null at range here is not a failure of Face C; it
is the measurement that forces Stage II of the long arc (§3) and justifies
the eventual kernel discussion. A phase-DEPENDENT contact outcome would be
gold: the first place hoard phase does work.

### P7 — Mobility frontier: can a bound thing move at all?
**Apparatus:** the (q_detune, amp, kx, hop_order) surface, focused where
Face C stalled. Strang (hop_order=1) is new since the DNLS route-B stall
arms ran [MEASURED: face A symmetrization is a real dynamical change —
vacuum ring retention 5×]. Meter: envelope rms (binding) AND fd dense
current + centroid track (translation) simultaneously.
**Bar:** any corner with rms growth <10% over the run AND centroid
translation ≥5 cells (route B managed ~1.3 cells before stalling).
**Decides:** whether the transport/binding antagonism is fundamental to the
channel or partly a sweep-schedule (Trotter pinning) artifact. In DNLS
physics, self-trapped modes pin on the discreteness barrier; face A already
showed the sweep order is a symmetry-breaker the schedule can repair.
[INFERENCE] Strang may shift the pinning threshold; if it does, "moving
particulate" opens with zero kernel changes. If no corner exists, record it
as a wall: particulates are pinned in the current channel, and motion must
come from elsewhere (the space sector, or door-mediated hop-along).
**Do not** chase exotic schedules beyond the three that exist; one clean
null on the existing surface is worth more than a schedule zoo.

### P8 — Creation robustness: is condensation generic or corner-tuned?
**Apparatus:** vary the initial condition class (spread packet → beam →
quench-style transient → noise_amp-lit bath) at fixed deep-corner law; 3
seeds each. Census as in P3.
**Bar:** condensation (yes/no + spectrum shape) across ≥3 qualitatively
different initial conditions.
**Decides:** whether Face C is an attractor of the dynamics (matter forms
whenever energy is dense enough — creation as a law) or a prepared-state
phenomenon. This is the difference between "the fabric makes matter" and
"we made a matter-shaped state." It also connects to the standing QUENCH
tension (#22): a second, independent creation mechanism changes how much
weight the B7-forward mechanism has to carry.

### P0 (debt, schedule alongside) — the armed L1-C + tolerance abx
The ρ_coh≈0.77 armed coherent bath anti-ignition bar is UNRUN, and no C≡Go
A/B exists at amp_nat>0 [MEASURED status]. Neither blocks P1–P8, but BOTH
block any *adoption* of amp_nat>0 into a law table. Run L1-C before the P4
results are used to argue for adoption; build the relative-tolerance abx
before any Face C number becomes a ratchet bar. This is ratchet hygiene, not
physics.

**Sequencing:** P1 → (P2, P5) → (P3, P4) → (P6, P7) → P8, with P0 threaded
in. P1 gates everything: if endurance fails, only P8 (creation-as-relaxation
is still interesting) and the writeup survive.

---

## 2. Scale correspondence (goal 2) — the only lawful gauge

This can be stated precisely, and most of it is a *negative* statement that
protects the programme from a category error.

**The no-background constraint forbids an absolute scale.** [INFERENCE from
§1, and standard dimensional analysis.] The simulation is a dimensionless
dynamical system; SI units are a background par excellence — a permanent
external ruler. There is no in-fabric observable that "sets the scale" in
SI, and looking for one is the wrong question. What the fabric CAN do is
define its **own** unit system from measured invariants, and then make
**dimensionless predictions** — pure numbers — which are the only things
that can ever be compared to reality without importing constants.

**The internal unit system (all three candidates already exist):**
- **Action:** the atom ε = A₀ω/2π [MEASURED, ħ-linearity ~1e-8 in the v89
  battery]. This is the fabric's ħ̂.
- **Speed:** c_eff, the field-sector wave speed [MEASURED instrument-grade
  since v89]. This is the fabric's ĉ, and sets the length↔time exchange.
- **Mass:** the Face C preferred hoard mass M* — IF P3 finds a seed-robust
  peak. This is why P3 is the scale experiment, not just a particulate
  experiment. (Fallback: the radiance fixed-point mass Em*≈1.55 is already a
  measured distinguished mass, but it is a *filter parameter* of the law,
  not a *selected state* — weaker.)

From (ε, c_eff, M*) every other quantity is a pure number: the hoard's
"Compton" ratio λ̂ = ε/(M*·ĉ) vs its measured envelope width; binding-energy
per mass of a P6 bound pair; ratios of masses in the hierarchy ladder;
line-frequency ratios in P5 spectra. **The correspondence programme is:
accumulate dimensionless numbers, and only at the very end ask whether any
subset matches reality's dimensionless numbers** (mass ratios like m_p/m_e ≈
1836, coupling ratios like α ≈ 1/137, spectral ratios like Rydberg series).
A match of even one nontrivial pure number would fix the SI gauge for free
(three unit choices then follow trivially); until then, publish everything
in fabric units and refuse the question.

**Concrete short-term deliverable (fold into P3/P5, no extra runs):** a
one-page table of the fabric's measured pure numbers — λ̂/width ratio,
spectrum peak ratios, Em*/cap, x̂*, the P2 existence-window edges in units
of M*. Cheap, honest, and it converts "what scale is this?" from a
philosophical embarrassment into a growing list of falsifiable numbers.

**[SPECULATIVE] flag:** any hope that these numbers land near real-particle
ratios is pure speculation at this point. The deliverable is the *method*
(the gauge protocol), which is sound regardless of whether the numbers ever
match.

---

## 3. Long-term arc (goals 3–5), with the load-bearing assumptions marked

The arc is NOT the pre-v89 stage table revived. It reorders around what v89+
actually measured: no far field from stable matter, no charge sector, cap
refusal as the exclusion seed, and the antagonism triangle. Five stages,
each with its failure point named.

### Stage I — Confirmed particulate (P1–P8 above)
Exit criterion: a seed-robust, isolable, medium-surviving hoard with a
measured mass window (ideally a preferred mass). Everything below assumes
this; if Stage I fails, the lane's residue is "creation-adjacent relaxation
physics" (still a result — record and pivot back to the retention wall).

### Stage II — Interaction at range: the space-sector completion
**The wall [MEASURED]:** nothing reaches out. Stable matter has no far
field; forces do not act at distance because no books are carried at
distance. P6 will almost certainly confirm contact-only interaction.
**The lawful route [INFERENCE]:** the missing far field was already
diagnosed in v89 (g4): it "awaits a stable particle's internal space
cycle." A confirmed Stage-I particulate is exactly the object that could
run a standing internal Es cycle (intake=outtake through its own door — the
flux-machine picture the user set in v89). The experiment: measure whether a
P4-surviving hoard in the honest medium runs a steady Es throughput, and
whether TWO of them bias each other's throughput at range (the v91 ASYM
"two eaters grow a channel between them" acceptance surface — registered,
never met). First attraction = two hoards' space books coupling.
**Where a kernel change becomes unavoidable [SPECULATIVE but likely]:** if
the space sector's transport (s_k/s_disp diffusion) provably cannot carry a
standing 1/r-like graded structure, the space sector may need its own
coherent channel — a "pass S" cousin of passes F and U. That is a major
face, to be designed like v93 was (consultation first). It must preserve:
default-0 byte-inertness against the 87-bar surface, conservation as a
theorem of the update, no per-cell immortal identity, and the measured g1
footprint at 0. **Do not open this before P6 + the throughput measurement
force it.**

### Stage III — Shells and exclusion (goal 3), WITHOUT importing charge
**Drop the v75 three-fabric plan.** It is a pre-v89 construct: imported
opposite-charge species on fixed fibers — banned twice over (§1 no imported
species; v93 §IV.2 no fixed fiber). The lawful ingredients on the board:
- **The well:** the hoard's own Es depression [g1, MEASURED for
  additive-era mass; P5(d) measures it for hoards]. Field waves refract in
  graded space [MEASURED: dark-void lensing, seed-robust]. So: **shells as
  standing field modes trapped in the particulate's space footprint.**
  Discreteness comes from the mode boundary condition — no charge needed
  for the FIRST shell result.
- **The occupancy rule:** cap saturation refusal is pitch-blind and exact
  [MEASURED: PAULI-0, v90]. A trapped field mode condensing at the hoard
  through the quantized door (ε atoms) inherits a discrete occupancy; the
  refusal at cap is the exclusion-adjacent mechanism. [SPECULATIVE in this
  combination — but every ingredient is separately measured, which is rare
  for a stage-3 plan in this programme.]
**First experiment when Stage I closes:** park a weak coherent field wave
on a P4 hoard; look for discrete trapped-mode resonances (transmission dips
/ standing-wave nodes) in its footprint. A single measured discrete bound
field mode = "shell 0."
**Failure point:** the hoard's Es depression may be too shallow/contact-
local to bind any mode (ASTRO measured footprints <1 cell for additive
matter). Then shells wait on Stage II's space completion — the well and the
far field are the same missing physics, which is economical: one wall, not
two.

### Stage IV — Substructure / flavor (goal 4)
The constraint is explicit [v93 §IV.2/IV.7]: substructure must be **mode
structure that dies with the cell** — internal harmonics, not fixed fibers.
The candidate: internal m-modes of ψ_m within one hoard (the m-spectrum
instrument already exists). The measured obstacle: hoards bind energy, NOT
phase — today they cannot hold internal circulation at all. So Stage IV has
a sharp, honest gate: **the phase-retention wall must fall first.** The two
lawful cracks in that wall found so far: topological support (closed cycles
hold winding in vacuum — face A) and schedule symmetrization. A hoard whose
envelope contains a closed high-|ψ| annulus (a self-trapped RING rather
than a ball — does the DNLS law have ring solutions? P2's seeding machinery
can ask cheaply) would combine both cracks. [SPECULATIVE.] If phase
retention inside dense matter never works, "flavor" in this substrate is
dead, and goals 4–5 need a different carrier — record the wall and consult
again rather than force it.

### Stage V — Atoms and elements (goal 5)
Only definable after III+IV. The honest formulation available NOW: an
"element" is a *discrete, reproducible spectrum of stable
hoard+trapped-mode configurations*, identified by dimensionless signatures
(mass ratios, line ratios — the §2 table, and the v91 spectroscopy
instruments spec.awk/prof.awk are ready for exactly this). Species-by-light
was already demonstrated once [MEASURED: v91 ASTRO — UUD vs D by line
ratios alone]. Stage V is therefore mostly *measurement infrastructure that
already exists*, waiting on physics from III/IV.

**Where no-background bites hardest across the arc:** (1) the scale gauge —
answered only by dimensionless ratios (§2); (2) anything at range — no
potential on a stage, so every force must be carried by the medium's books,
and the books provably stay home for stable matter (Stage II is the arc's
critical path); (3) identity — "the same hoard" over time is a centroid-
continuity convention, not an ontological fact (the v91 IDENTITY lane's
lesson: identity preserves, does not create — do not let particulate claims
quietly assume gid-style identity).

---

## 4. Recommendations

**Highest-value next moves (in order):**
1. **P1+P2+P4 as one bundle** (endurance, isolation, honest-medium
   survival). Cheapest possible reckoning of whether Face C is matter. Every
   week spent elsewhere before this is at risk of building on a transient.
2. **P3, the radiance mass-selection sweep.** A seed-robust preferred mass
   is simultaneously the strongest particulate claim, the fabric's mass
   unit, and the entry to the only lawful scale correspondence (§2). It
   converts the metaphysical goal-2 question into a table of measurable
   pure numbers.
3. **P6+P7 (interaction + mobility).** These two nulls-or-not determine the
   entire long arc's shape: contact-only + pinned → the space-sector
   completion (Stage II) becomes the programme's critical path and the next
   consultation topic; any interaction at range or any moving bound corner
   → the current channel carries further than expected, no kernel change
   needed yet.

**Do NOT yet:**
- **No kernel changes** — every short-term experiment is config-only. The
  pass-S/space-channel discussion is premature until P6 + the Es-throughput
  measurement force it, and it should get the full v93 treatment
  (design doc + external review) when it comes.
- **No v75 three-fabric revival, no imported charge.** Banned-era shape.
  The charge-analog must be found (winding sign, pitch sign, or something
  the fabric offers), and Stage III's first shell result does not need it.
- **No SI mapping, no matching to real particle masses** until the §2 pure-
  number table has entries. Premature correspondence claims are the fastest
  way to make the programme unfalsifiable.
- **No adoption of amp_nat>0 into a law table** before P0 (armed L1-C +
  tolerance abx) is green. Ratchet discipline is the programme's spine.

**The one-sentence path:** confirm the hoard is matter (P1/P2/P4), let
radiance mint the mass unit (P3), measure whether hoards touch and move
(P6/P7) — and let those answers, not ambition, decide whether the next face
is the space sector's coherent channel, which is where the road to shells,
forces, and elements most likely runs.
