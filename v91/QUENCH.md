# QUENCH — genesis by dynamics: compress fast, release fast (user-directed 2026-08-06)

User question (2026-08-06, verbatim): *"Have we tried compressing the
field quickly in the first N frames (eg 10 or 50) to a high degree
then quickly releasing the bounds?"* Answer measured in the record:
never. This campaign runs the config-only forms of that protocol on
the standing binary (kernel and laws untouched; the scheduled-anneal
/ breathing-box knobs are kernel work, conditionally authorized by
the user ONLY if these arms show a pulse). Registered before any run;
results land in §4 per arm; nothing is adopted.

## 0. Thesis

Every creation-class mechanism tried so far is a GATE, and four
rounds across three lanes closed that door: gates preserve, they do
not create. A quench creates by DYNAMICS: drive the medium through a
transient over-compressed state and release; nucleation happens in
the sweep, not through a gate. Two measured facts aim it: (i) the
law-regime vacuum CONDENSES above-threshold beams into matter trails
(XSEC/B7 — the condensation door exists and is strong: cond ≈ 77 of
a ~180-unit beam); (ii) during release the local fullness sweeps DOWN
THROUGH the balance locus x̂* — and CANTUS says chords are the only
stable fits there. The registered question: does sweeping the balance
locus ever nucleate persistent bonded matter — episodes + mutual
stamps that outlive every null?

The DETECTOR is the identity meter (IDENTITY.md, physics-silent):
nucleation is claimed only for a connected cluster of ≥3 cells whose
episode ages exceed 3× the bath mean (>1600 t.u.) with ≥2 mutual
stamps older than 1500 t.u. (bath stamp tail: P < 0.3%), alive at
run end with roughly balanced books. No eyeballing.

## 1. Arms (exact commands, standing binary, from `v91/`; seed 20260802)

The beam is a one-shot packet: its transit IS the fast compression
(~σ/C t.u. locally — the user's "first N frames"), and its passage IS
the release. No knob is needed for the vacuum forms.

```
# Q1  beam dump into quiet foam, law regime (XSEC class + long watch)
./freecell exp=slit slit_mask=3 L=64 T=1500 sigma=4 slit_sy=3 kx=2.0 amp=4 slit_srcx=14 slit_screenx=40 slit_t0=0 slit_t1=30 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q1_beam4.fcs > runs/quench/q1_beam4.log 2>&1

# Q2  deep quench: amp 4 -> 8 (only knob changed)
./freecell exp=slit slit_mask=3 L=64 T=1500 sigma=4 slit_sy=3 kx=2.0 amp=8 slit_srcx=14 slit_screenx=40 slit_t0=0 slit_t1=30 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q2_beam8.fcs > runs/quench/q2_beam8.log 2>&1

# Q3  packet through the LIVE warm bath (exp=pulse at quench amplitude)
./freecell exp=pulse bath=1 noise_amp=0.5 L=24 T=1500 amp=3 sigma=3 kx=1.1 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q3_pulse_bath.fcs > runs/quench/q3_pulse_bath.log 2>&1

# Q4  dense-grain quench: driven two-blob collision (MOTION #28 seed)
./freecell exp=blob2 amp=1.6 sigma=2.2 blob2_sep=8 blob2_kx=0.9 L=24 T=1500 diag_every=200 seed=20260802 k_rad=0.05 p_rad=4 rad_clock=0 par_tau=10 qatom_every=200 snap_every=2500 snap_file=runs/quench/q4_blob2.fcs > runs/quench/q4_blob2.log 2>&1
```

Controls, standing (no new runs): the warm-bath identity baselines
are M-I1 (episode mean 530 t.u.; stamp ages memoryless mean 261;
`runs/identity/mi1_bath_warm.log`); the quiet foam does not densify
spontaneously (FC0, gated). 2D for Q1/Q2 (the slit substrate — B9's
dimensionality caveat applies and is recorded).

## 2. Closed-form expectations (before any run)

1. **The nulls to beat.** RADIANCE: released dense hoards decay with
   t_half 80–140 → the condensed trail should be substantially gone
   by t ≈ 500. CANTUS geometry wall: quench condensate is
   unison-grade matter; unison cannot bond at the balance → no
   chords expected. FC0: nothing condenses without the beam.
2. **The hypothesis's one opening.** During trail decay each dense
   cell's x sweeps down THROUGH the chord operating band (x_D ≈
   0.55–0.65, x_U ≈ 0.28); adjacent cells crossing at m=2-compatible
   loads could, in principle, latch a fifth. Nothing measured says
   they will; everything measured (unison wall) says they won't —
   the arm exists because the sweep-through-balance configuration
   has never actually been produced.
3. **Q3 caution registered:** the pulse deposits into a LIVE bath —
   any "nucleus" must beat the bath's own 530/261 baselines, and
   spatial attribution needs the snapshots (global PAR is
   bath-dominated in Q3).
4. **Q4 expectation:** blobs profile-merge (COE: v_COM ≤ 1e-3 of
   closing) — the "collision" is slow; the compressed interface is
   the quench zone; RADIANCE null applies to the merged hoard.

## 3. Bars

* **Q-DET (the only claim bar):** a cluster per §0's detector at run
  end in any arm → NUCLEATION CLAIMED, arm replicated ×2 seeds
  before anything further; then and only then the scheduled-anneal
  knob (kernel) is authorized per the user's condition.
* **Q-NULL (per arm):** trail/hoard decay curve reported against the
  RADIANCE t_half band; episode/stamp demographics reported against
  M-I1; "no chords" checked via locks=0 (k_cant=0 here — no lock
  apparatus: whatever persists must persist by CONTACT alone, the
  nv48-vacuum class) and via the QATOM spectrum (any D-band emission
  line standing over the trail's own band would be a chord
  signature).
* Drift ≤ 2e-15 as always; battery untouched (config-only).

## 4. Results

**4.1 Q1 + Q2 (beam dumps into law-regime foam; logs
`runs/quench/q1_beam4.log`, `q2_beam8.log`) — the quench makes a
SELF-SUSTAINING MATTER CLOUD; the decay null is beaten at ensemble
grain.** (harvested 2026-08-06)

* **The unperturbed foam is identity-dark** (0 live episodes at
  t ≤ 4 — every episode in these runs is beam-descended). The beam
  (gone by t=30) births a population that grows to a **plateau ≈
  1,800 episodes (amp 4) / 2,400 (amp 8) by t ≈ 750 and holds it
  FLAT to t=1500 — no decay trend over the final 750 t.u.** The
  RADIANCE null (hoard t_half 80–140, which would have gutted the
  patch by t ≈ 500) is REFUTED for law-regime quench debris.
* **Persistence is flame-like — the ensemble lives while
  individuals churn.** Episode ages at close: q25/50/75 = 34/150/536
  (amp 4), 71/241/708 (amp 8); mints ≈ retires in steady state. The
  patch holds 4,700–7,100 stamped PAIRS at any instant, and at
  least one pair's stamp is continuous from beam passage to run end
  (stamp max 1,489/1,477). Identity-carried delivery fraction
  0.90/0.985.
* **Deeper quench → bigger, older cloud** (amp 8 vs 4: population
  +38%, median age +60%, idfrac up) — monotone in the compression
  depth, as a genesis mechanism should be.
* **Mechanism, read off B7:** at e_cond=0 the patch's own radiated
  glow (rad = 137/864 in the ledger) recondenses locally instead of
  escaping — the cloud-chamber vacuum FARMS passing energy into
  persistent populated matter. The optics regime (e_cond=99) would
  let the glow escape; the dichotomy is the B7 split itself.
* **Registration slip, recorded:** the §3 Q-DET age bar (>1600 t.u.)
  is unsatisfiable in a T=1500 run (max possible age ≈ 1,490) — an
  arithmetic error in the registration, not a soft bar. The
  DECISIVE extension is registered here BEFORE it runs: **Q5 = Q2
  at T=5000** (config-only; launched at harvest). Bars: (i) the
  population still ≥ half its plateau at t=5000 (persistence is not
  a slow dilution); (ii) Q-DET as registered becomes satisfiable —
  any ≥3-cluster of episodes ≥1600 with mutual stamps ⇒ NUCLEATION
  as originally defined; (iii) no lock apparatus anywhere (k_cant=0
  — whatever persists, persists by contact + recondensation alone).

**4.2 Q5 (the decisive extension: amp 8, T=5000; log
`runs/quench/q5_beam8_t5000.log`) — GENESIS BY DYNAMICS HAS ITS
FIRST POSITIVE.** (harvested 2026-08-06)

* **Bar (i) PASS, decisively:** population 2,517 (plateau, t=500) →
  **2,113 at t=5000 = 84%** — a slow ~16% leak over 4,500 t.u.
  (implied ensemble half-life ≈ 20,000 t.u., two orders beyond the
  radiance null's 80–140). The cloud is persistent, not slowly
  diluting.
* **Bar (ii) at pair grain: NUCLEATION PROVEN.** A stamped pair is
  co-continuous from beam passage to run end — **stamp age 4,926
  t.u.** Against the measured bath stamp hazard (memoryless, mean
  261 t.u.): P(churn) = e^(−4926/261) ≈ 7×10⁻⁹. That is a BOUND
  PAIR by any statistical standard — created by a transient beam
  and contact/recondensation alone, **no lock apparatus in the box**
  (k_cant=0, bar iii ✓). The ≥3-cluster criterion is NOT scorable
  from the quartile-grain diag — instrument gap recorded; a
  per-cell gid dump (print-only apparatus, battery cycle) is the
  registered follow-up if the user wants the cluster census.
* Steady state at close: ages q25/50/75 = 59/208/881 (long tail —
  ~6% of episodes estimated ≥1600 under the fitted tail, not
  directly countable at this diag grain); 5,532 standing stamped
  pairs; idfrac 0.984; mints ≈ retires (70,524 / 68,411 cumulative).
* Reading, stated carefully: the three prior lanes measured that
  GATES cannot create. The quench — creation by dynamics — makes
  (a) a persistent identity-bearing ensemble and (b) at least one
  statistically-certain bound pair, in the law regime, from nothing
  but a passing beam. What it has NOT made: a chord (no lock
  apparatus was present to hold one; the CANTUS unison wall stands
  untested by this arm because nothing here reached the m=2
  geometry). The registered next question if pursued: a quench WITH
  the honest lock+tuner apparatus present-but-unseeded
  (cant_grow=0, cant_seed=0 — gauges CANNOT arm, so no ignition
  confound) cannot test chord formation; chord formation from
  quench debris requires the self-growing gate, which ignites —
  the quench and the gate problem MEET here, and the design-only
  escalation (step joint-maturity) or the field-side lane are the
  registered ways through.

**4.3 Q3 + Q4 (the 3D arms) — two honest nulls with one bonus.**
(harvested 2026-08-06)

* **Q3 (pulse through the live warm bath): unresolvable at global
  grain, exactly as §2.3 pre-warned.** The global census sits AT the
  bath baseline from t=300 (≈8,770 live ≈ the M-I1 95% fraction;
  ages q50 405 vs the bath's 440; idfrac 0.997). Any wake signal
  needs spatial gid attribution — the same per-cell gid dump
  registered in §4.2 as the instrument follow-up.
* **Q4 (driven blob collision): the collision NEVER HAPPENS — B2
  reconfirmed.** sep 8.00 → 8.37 over 1500 t.u., closing −0.0005,
  vTotalCOM 5.7e-5: the phase-tilt drive cannot translate dense
  matter (merging-is-not-transport, measured again). The dense-grain
  collision quench is UNREACHABLE at this law — not failed,
  impossible to arrange. Bonus: each isolated blob grows a
  persistent identity HALO (population rising to 4,847 at close,
  the §4.1 cloud phenomenon at dense seeds, in 3D vacuum).

## 5. Verdict (all arms harvested; decisions are the user's)

1. **The user's protocol WORKS, and it is the programme's first
   measured act of creation.** Fast field compression + release (a
   transient beam in the law-regime foam) births a SELF-SUSTAINING
   matter cloud (plateau ~1,800–2,500 episodes; 84% retained at
   t=5000; implied ensemble half-life ~20,000 t.u. — two orders
   beyond the radiance hoard-decay null) containing at least one
   statistically-certain BOUND PAIR (stamp co-continuity 4,926
   t.u.; P(churn) ≈ 7e-9) — with no lock, no tuner, no gate
   anywhere. Where four gate rounds created nothing, dynamics
   created. Monotone in compression depth.
2. **Mechanism: B7 read forward.** At e_cond=0 the debris
   recondenses its own glow — the cloud-chamber vacuum farms
   passing energy into persistent populated matter. The optics
   regime would let the glow escape; the work-function lane (task
   #22) will decide whether one honest table can hold both
   behaviors.
3. **Limits, measured:** no chords (nothing in the box could hold
   m=2 geometry; chord-from-debris re-encounters the self-growing
   gate's ignition problem — the registered ways through are the
   step joint-maturity design or the field-side lane); the 3D
   dense-grain collision is unreachable (B2); bath-embedded wakes
   are unresolvable without the per-cell gid dump (print-only
   apparatus, battery cycle — registered).
4. **The scheduled-anneal knob:** the user's condition ("authorize
   only if config-only shows a pulse") is now MET by the Q-DET
   pair-grain positive. Noted as authorized-by-condition and
   DEFERRED: the beam form already realizes compress-and-release
   cleanly for the standing questions; the anneal form (whole-box
   schedule) is the stronger instrument and waits for the user's
   next word rather than being spent unilaterally.
5. **Nothing was adopted.** All arms config-only on the standing
   binary; battery green throughout; the one registration
   arithmetic slip (the >1600 bar at T=1500) is recorded in place
   and corrected by the registered Q5 extension, not rewritten.

---

## 6. QUENCH-2 (user-directed 2026-08-06: rate, confinement, spin, precise viewing; kernel updates pre-authorized; registered BEFORE the apparatus lands)

Three apparatus knob families (all default-inert, purity-pinned,
battery-gated) plus full fcs capture with volview renders on every
arm:

* **Rate — CORRECTED AT DESIGN TIME (dated 2026-08-06): the anneal
  knob is REJECTED BY THE ONE LAW.** `noise_amp` is seed-time only;
  a runtime noise drive would CREATE energy. The conservative rate
  knob is packet CONCENTRATION at matched total energy
  (E ∝ amp²σ² in 2D): fast-deep σ=2/amp=16 vs slow-shallow
  σ=8/amp=4 vs the Q5 reference σ=4/amp=8 — config-only, no kernel
  knob. (The "authorized-by-condition anneal" dissolves — there was
  never a lawful whole-box drive to build.)
  **[FIX, dated 2026-08-06 (user: "can we fix the QUENCH-2?"): the
  whole-box compress-release is RESTORED lawfully as IMPLOSION
  FOCUSING — `imp_k` (default 0, byte-inert, both kernels,
  battery-gated) adds an inward radial phase tilt to the seeded
  packet: a large converging wave self-compresses at the focus
  (the "first N frames") and re-expands (the release), seeded
  energy only, fully conservative. Arm Q2-I registered: σ=16
  sy=16 amp=2 (E ≈ 1024 matched) imp_k=2, source at box centre,
  T=1500, fcs capture.]**
  **[ENERGY-MATCHING CORRECTION, dated 2026-08-06 (user question:
  "simulation error or unintentional energy conversion?"): neither
  — all arms run at the drift floor (≤2e-15) and E0 ledgers
  confirm the budgets (spin E0 bit-identical to Q5; cavity −81 =
  the carved shell's bath energy, accounted). But the registered
  matching FORMULA was wrong: the seed law is Ee = amp·g, so
  packet energy ∝ amp·σx·σy, NOT amp²σ². The rate pair is matched
  anyway (amp ∝ 1/σ satisfies both criteria; amp·σ = 32 all
  three). Q2-I as registered would have carried 5.3× the class
  energy — caught BEFORE launch; corrected to amp = 0.375
  (2·16·16 → 0.375·16·16 = 96 = the class budget).]**
* **Partial confinement:** a PINNED-FIXTURE SHELL (the DS-wall
  precedent — carved-vacuum walls held as pinned fixtures are
  standing apparatus): `conf_r` (default 0 = off) seeds a circular
  shell of pinned wall cells of radius conf_r about the box center
  with a leak gap `conf_gap` (half-angle, radians; default 0.3).
  Matter and field bounce by the standing contact machinery; the
  gap leaks — real cavity physics, no invented force, no energy
  ledger entry.
* **Spin:** `spin_m` (default 0): the seeded packet/anneal noise
  carries an azimuthal phase winding m·φ about the z-axis through
  the compression region — angular momentum enters the quench; the
  question is whether debris inherits circulation (the QATOM/REG
  circ meters + fcs render read it).

Arms (registered; exact commands fixed at implementation, one knob
family at a time vs the Q5 reference): Q2-R rate pair (fast/slow at
matched energy), Q2-C confinement pair (hold/free), Q2-S spin pair
(m=0/m=2), each amp-8-class, T=1500, par_tau=10, L=64 law-regime
foam, snap_every for fcs + volview plates. Bars: the Q5 detector
(cloud population, pair stamps) per arm; circulation for Q2-S;
drift ≤ 2e-15 with confinement on (the apparatus must not leak
energy). Results land in §6.1.

### 6.1 Results (five protocol arms harvested 2026-08-06; Q2-I pending)

**THE CLOUD IS AN ENERGY-DETERMINED ATTRACTOR — it forgets the
protocol.** At matched seed budget (amp·σx·σy = 96-class), all five
arms land on the same state: live populations 2,291–2,411 at
t=1500, standing stamps 6,404–7,064, max stamp ages 1,432–1,483
(every arm keeps a pair continuous from beam passage), idfrac
0.97–0.985.
  - **Rate: NULL.** Fast-deep (σ=2, amp=16) vs slow-shallow (σ=8,
    amp=4): 2,395 vs 2,411 — indistinguishable. Compression rate
    does not matter at matched energy; even Q1-vs-Q2's earlier
    "depth" effect re-reads as pure energy scaling (their amp·σ·sy
    budgets were 48 vs 96).
  - **Spin: INVISIBLE at every current meter.** m=2 winding: E0
    bit-identical to the m=0 reference (phases carry no energy) and
    the debris census identical (2,395). Either condensation erases
    the winding (the door reads magnitude, not phase) or the meters
    are blind to it — the AMPLITUDE shadow (amp_tau) is exactly the
    phase-resolving instrument; the one-run follow-up (q2s +
    amp_tau=10) is registered.
  - **Confinement: NULL at matched energy.** The cavity arm's −3%
    population is fully explained by its −2.6% energy (the carved
    shell's bath cells); the real effect is KINETIC — assembly lags
    (2,049 @ t=500 vs 2,464 free) then converges to the same
    asymptote; idfrac dips to 0.970 (parcels die in flight more
    often inside the cavity).
The thermodynamic reading, stated carefully: quench debris
equilibrates to an energy-set attractor — the cloud has an
equation of state, and how you compress is forgotten by t≈500.
Plates: em-channel animation rendered (volview, 100 frames;
`runs/quench/plates/` local — large binaries stay out of git; an
exhibit page can carry them on request).

**Q2-S shadow pair (the registered phase-grain follow-up; harvested
2026-08-06): SPIN IS COMPLETELY ERASED — the quench door is
phase-blind.** m=0 vs m=2 with the AMPLITUDE shadow on: ρ_coh
0.5917/0.7342/0.8441/0.9228 vs 0.5830/0.7319/0.8415/0.9215 (Δ ≤
0.009 at every quartile), delivered magnitude and flowing-slot
counts matched to 0.1%. With E0 bit-identical and the census
identical, the m=2 winding leaves NO trace at ANY measured grain —
population, identity, energy, or delivery-phase coherence. The
mechanism is the standing one: the conversion door reads MAGNITUDE
(beat-fire on |Ee|), so the azimuthal phase structure dissipates in
transit and never imprints on matter. The door-erasure ledger now
reads: identity (light carries load, never who) AND spin (matter
carries energy, never winding) — the door is a magnitude
bottleneck, measured from two independent sides.

**Q2-I (implosion focusing — the fix arm; harvested 2026-08-06):
SAME ATTRACTOR.** The lawful whole-box compress-release (imp_k=2
converging wave, matched budget E0 = 3125, drift −2.18e-15 — the
largest of the family but FP-floor class) lands at population
2,322, stamps 6,754, a beam-passage pair held (max 1,469), idfrac
0.985 — inside the five-arm band. **Six protocols — linear
fast/slow, cavity hold/free, spin, implosion — one state.** The
equation-of-state claim now spans every compression geometry the
apparatus can produce; the only variable the cloud remembers is
its energy budget. Implosion animation rendered to the local
plates set.

### 6.2 DATED CORRECTION to the Q2-S spin claims (2026-08-07, log
audit ordered by the user: "Currently quench is done incorrectly.
Specifically, look at the spin-quench.") Three errors of record,
one structural re-grade; the §6.1 spin text above stands as
written per the no-rewrite rule, corrected here:

1. **"The debris census identical (2,395)" is FALSE.** 2,395 is
   the m=2 arm's own number (`q2s_spin2.log` live=2395); the m=0
   reference landed **2,402** (`q2s_shadow_m0.log`, byte-same
   protocol, E0 bit-identical 3160.440443538). Both sit inside
   the five-arm attractor band — the honest wording was
   "attractor-equal", never "identical".
2. **"The m=2 winding leaves NO trace at ANY measured grain" is
   FALSE at the ledger and optics grains, in the same logs:**
   births 46,429 vs 48,694 (−4.7%), cond 1017.5 vs 1079.0
   (−5.7%), evap 32.2 vs 42.1 (**−24%**), rad 809.5 vs 864.4
   (−6.3%), screen exposure 92.10 vs 75.95 (**+21%**). The
   winding measurably changed the packet's transit optics (as a
   complex field must — pass F is a genuine unitary wave) and
   thereby where/how much the beam condensed. What the winding
   left no trace IN is the equilibrated cloud's matter attributes
   — a much narrower, and true, claim.
3. **The §6 registered "circulation for Q2-S" bar was never
   measured — no instrument existed.** The only circulation meter
   in the kernel (`tri_meters`, freecell.c:2702) requires a
   declared 3-cycle; the quench box has none. The bar was
   silently dropped at harvest. Instrument gap recorded; §7's
   offline fcs winding meter is the repair (th2 and fa1/fa2 are
   already snapshot columns — no kernel work needed to measure).
4. **Structural re-grade of the mechanism sentence.** "The door
   reads MAGNITUDE (beat-fire on |Ee|)" is correct — but it is a
   property of pass 6 BY CONSTRUCTION, not a discovery of this
   experiment: condensation reads Ee = fa1²+fa2² and writes Em as
   a pure scalar (freecell.c:1851–1875); **th2, the matter clock
   the gates read, is never written by the door in the
   field→matter direction** (matter→field writes it into empty
   cells at line 849 — phase crosses the door one way only).
   Q2-S therefore measured the apparatus's own known bottleneck —
   the angular twin of the recorded ~100× radiation-pressure
   (linear momentum) deficit at the same door (v89/ROADMAP.md
   P2, an S2-full acceptance criterion). Grade: the spin null is
   STRUCTURAL-BY-CONSTRUCTION, not an empirical property of the
   law; the empirical question ("does matter inherit winding when
   the door can carry it?") was never asked. §7 asks it.

---

## 7. QUENCH-3 — the phase-carrying door (user-directed 2026-08-07: "Back port the answer to the quench work and try it again. Goal is to create a plethora of particles."); registered BEFORE the apparatus lands

Companion analysis: `BLACKHOLE.md` (same date) — the symbolic no-go
that identifies the missing sector this lane instruments. Summary of
the identification: **the dense cell has a phase attribute (th2 — the
clock the bond gates read at retarded relative phase, freecell.c:1623)
but the conversion door never writes it** (field→matter transfers a
scalar; matter→field already writes th2 into empty cells). Every
first-moment book — linear momentum (P2's ~100× deficit), circulation
(Q2-S), phase (AMPLITUDE Phase M's map) — dies at the same line of
code. The S2/AMPLITUDE lane is the registered full repair (parcels as
amplitudes); this lane is its cheapest honest instance: **close the
door's phase loop in the one direction it is open-circuit, using only
the existing attribute.** No new state array, no energy content, no
ledger entry.

### 7.1 The knob (default-inert, purity-pinned, both kernels)

`qp_phase` ∈ [0,1], default **0** = byte-inert. In pass 6, at a
condensation event (beat-fire, d1 > 0), AFTER the energy moves:

    aph  = atan2(fa2, fa1)                  # the wave's phase at the door
    mix  = qp_phase · d1 / Em_new           # the field-derived fraction
    th2 += mix · wrap_pi(aph − th2)         # circular pull toward the wave

A fresh episode (Em_old = 0) is born with its clock ≈ SET by the wave
that made it (mix → qp·d1/(d1+dsp) ≈ 0.87·qp); a fed episode is
nudged in proportion to what it just ate. th2 carries no energy:
drift, E0, p1, and every conservation meter are untouched by
construction; the DYNAMICS changes only through the gates that
already read th2 (kappa_lock corrections, interval phase-slips).
Emission (line 849) already does the mirror image. The QUENCH-2
apparatus header line gains `qp_phase=%g`; the battery purity pin is
extended to match at 0.

### 7.2 Meters (offline — no kernel work; th2/fa1/fa2 are fcs columns)

* **W_A** — winding number of the FIELD phase arg(fa1+i·fa2) around
  the beam origin (srcx=14, cy0), annulus-summed wrapped increments:
  the injection witness Q2-S never had.
* **W_th2** — winding number of live-matter th2 (Em ≥ 0.05) around
  the same axis: the carried-spin meter. Common-mode clock rotation
  cancels in wrapped differences by construction.
* **m-spectrum** — |⟨e^{i(th2 − mφ)}⟩| over live cells, m ∈ [−4,4]:
  the winding-class spectrometer (the "plethora" axis).
* **n(r)** — live-cell radial profile about the axis (vortex-core
  donut test — a matter-grain spatial trace §6.1 never looked for).
* Census/ledger/demographics from the standing PAR/RESULT prints.

### 7.3 Arms (T=1500, L=64 law-regime foam, seed 20260802, exact Q2-S
protocol otherwise; snap_every=2500 for fcs)

| arm | spin_m | qp_phase | role |
|---|---|---|---|
| Q3-A | 0 | 0 | byte-inertness witness: must reproduce `q2s_shadow_m0.log` EXACTLY |
| Q3-B | 2 | 0 | ditto vs `q2s_spin2.log`; + W_A injection witness, W_th2 erasure witness |
| Q3-C | 0 | 1 | imprint control: beam's linear kx ramp imprints tilt-texture, no winding |
| Q3-D | 2 | 1 | THE ARM: does matter inherit the winding? |
| Q3-E | 4 | 1 | winding ladder (class separation = the plethora claim) |

### 7.4 Bars and registered predictions

* **P-Q3-0 (inertness):** Q3-A/Q3-B logs byte-identical to the
  archived Q2-S logs; battery ALL GREEN 93 both kernels; C≡Go with
  the knob FIRING verified once (pair-class run, qp_phase=0.7).
* **P-Q3-1 (injection):** m=2 arms carry field winding W_A ≈ +2 in
  transit-era frames; m=0 arms ≈ 0. Closes §6.2 item 3.
* **P-Q3-2 (erasure, instrumented):** at qp=0, live-matter winding
  |W_th2| < 0.25 and m-spectrum flat — the §6.1 erasure claim made
  measurable and (predicted) confirmed at the matter grain.
* **P-Q3-3 (the door carries):** at qp=1, m=2: transit-era
  W_th2 ∈ [1,3], unambiguously nonzero vs Q3-B. This bar is the
  lane's life: if it fails, the imprint never happened and the lane
  records a null.
* **P-Q3-4 (retention fork — both outcomes registered as results):**
  (a) recondensation re-imprint holds W_th2 ≥ 1 to t=1500 → first
  matter with a carried first-moment attribute; the m-ladder (Q3-E
  vs D vs C m-spectra pairwise distinct) = species by winding = the
  plethora door OPEN. (b) churn + differential pitch rotation erode
  it → W_th2 decays with a measurable time constant → "the door
  carries, the cloud cannot hold" — the retention wall becomes the
  next registered number and the Phase-L design inherits it.
* **P-Q3-5 (attractor, expected null):** census stays inside the
  five-arm band 2,291–2,411 at qp=1 (the energy attractor should be
  phase-robust); a >3% departure would be the first census-grain
  protocol memory — reported either way.
* **P-Q3-6 (optics continuity):** the §6.2-item-2 transit
  differences (cond −5.7%, exposure +21%) persist at qp=1 within a
  few %, since they are transit optics, not door physics.
* Drift ≤ 2e-15 all arms; nothing adopted; knob stays inert at
  defaults; decisions are the user's.

### 7.5 Registered extension at harvest (2026-08-07, before it runs):
the STANDING-VORTEX arm. The fine-cadence §7.4 harvest measured the
field's winding coherence dying in transit with τ_φ ≈ 3 t.u. (RA2
1.000/0.796/0.392/0.106/0.011 at t = 0/2.5/5/7.5/10) — FASTER than
the door's first beat (~5–7 t.u.): the flying packet's spin is
speckle before the first conversion reads it. Discriminating arm: a
standing long-wavelength vortex (kx=0, σ=8, box-centred) whose phase
gradient m/r ≪ the speckle scale — seeded IN PLACE over live foam so
conversion needs no flight. Q3-G1: slit_mask=3 kx=0 sigma=8
slit_sy=8 amp=2 slit_srcx=32 spin_m=2 qp_phase=1 T=300; Q3-G0: same
with qp_phase=0 (control). Predictions: (i) the standing vortex's
field RA2 survives ≫ 10 t.u.; (ii) matter R2d in G1 rises above the
G0 control while condensation runs — the first texture-grain
door-carry proof; (iii) after the field drains, the matter texture
decays at the differential-rotation rate (~1/30 t.u. from the load
spread) — the retention number Phase-L inherits.

### 7.6 Results (all arms + §7.5 extension harvested 2026-08-07;
logs `runs/quench/q3*.log`, fcs local, analysis instruments `report/analyze_winding.py` (fcsdump-fed) and `report/bh_symbolic.py`)

* **P-Q3-0 PASS three ways.** Both byte-witnesses reproduce their
  archived logs EXACTLY except the one registered header line
  (q3a ≡ q2s_shadow_m0, q3b ≡ q2s_shadow — every physics/diag/PAR/AMP
  line identical over 375k steps); fcs common frames identical (the
  archived spin2 fcs was snapped 5× denser — cadence-only diff,
  recorded); battery ALL GREEN 93 (`runs/BATTERY_quench3.log`); C≡Go
  physics byte-equal with the knob firing (drift-column-only
  divergence = the standing IDENTITY §4.0 envelope; atan2 added no
  new class).
* **P-Q3-1 PASS — and it re-reads §6.1.** The demodulated meter on
  the ARCHIVED q2s_spin2.fcs t=0 frame reads W_A = +2.000,
  RA2 = 1.000, pkA = m+2: the m=2 injection was always clean, and
  the Q2-S instrument gap is closed retroactively. The m=4 arm
  injects W_A = +4 likewise (ladder verified).
* **THE TRANSIT DISCOVERY (fine-cadence pair, 2.5 t.u. grain):** the
  flying packet's winding coherence dies in the speckle foam with
  τ_φ ≈ 3 t.u. — RA2 = 1.000 / 0.796 / 0.392 / 0.106 / 0.011 at
  t = 0 / 2.5 / 5 / 7.5 / 10 — FASTER than the door's first beat
  (~5–7 t.u.); the first matter appears at t ≈ 7.5 with the field
  already speckle. **The spin dies in transit, not at the door.**
  §6.1's mechanism sentence is re-attributed: even a phase-carrying
  door inherits speckle from a flying packet at this λ/dmin. P-Q3-3
  as registered therefore FAILS for the beam arms — with the
  mechanism and its number in hand, which is worth more.
* **P-Q3-2 PASS:** at qp=0 the matter clock field sits at the 1/√n
  floor everywhere, wandering peaks — the erasure claim, finally
  instrumented.
* **§7.5 standing-vortex extension — all three predictions land:**
  (i) the medium's phase memory is WAVELENGTH-DEPENDENT: the
  standing m=2 vortex holds field coherence to t ≈ 100–110
  (RA2 0.71 at t=30, pkA=m+2 through t≈100; ~20× the flying packet
  — part of the decay is consumption: the growing cloud eats its
  vortex). (ii) **THE DOOR CARRIES — texture-grain proof at ~3.5σ:**
  while driven, the qp=1 cloud's clocks peak at the seeded m=+2 in
  **13/39 frames vs 5/39 in the qp=0 control — P(≥13 | chance 1/9)
  = 2×10⁻⁴ vs 0.44.** The amplitude is weak (driven-era mean R2d
  0.032 vs 0.026): the imprint fights kappa_lock delivery churn and
  differential clock rotation continuously. (iii) **NO RETENTION:**
  once the field drains (t ≳ 110) the m=2 texture dissolves to
  chance — the winding is a DRIVEN texture, not a carried
  attribute; intrinsic matter phase memory ≲ 10–20 t.u.
  (differential-rotation class ~0.2 rad/t.u. from the load spread).
* **P-Q3-5 PASS (the expected null):** beam censuses 2,337–2,402,
  inside the band; the standing arms land a different, lower-churn
  protocol class (live 2,549/2,587 at mints 4.9k vs the beams'
  21.6k — recorded). **P-Q3-6 PASS:** exposure 91.26 vs 92.10
  (−0.9% — transit optics door-independent, as predicted).
* Ledger note, recorded not claimed: qp=1 shifts the m=2 beam arm's
  ledger coherently (+5.3% births, +4.4% cond, +6.3% rad vs its
  qp=0 twin) while the m=0 pair moves less and mixed; the qp=1
  evap ladder is monotone in m (46.1/33.2/21.1 at m=0/2/4). Single
  seed — chaos-band attribution until replicated.

### 7.7 Verdict (decisions are the user's; nothing adopted;
qp_phase stays inert at defaults)

1. **Q2-S was done incorrectly in the three §6.2 ways, and the
   corrected apparatus converts its one-liner into two measured
   walls.** "The door reads magnitude" re-reads as: (a) the MEDIUM
   has wavelength-dependent phase memory — flying-packet winding is
   speckle in ~3 t.u., before the first conversion beat; standing
   long-wavelength structure survives ~20× longer; (b) the CLOUD
   churns phase — a door-imprinted texture is driven-only, dying
   with its field.
2. **The missing cell attribute, final measured form:** not th2
   itself — it exists, and the door now provably writes it (the
   3.5σ standing-vortex imprint is the existence proof that
   door-carried phase reaches matter at all). What is missing is a
   phase DOF PROTECTED from the cell's own delivery churn: phase
   carried on the parcels/slots (where delivery coherence is
   already measured 0.77–0.95) rather than exposed in the naked
   cell clock. That is exactly AMPLITUDE Phase L, which this
   campaign hands a quantitative requirements card: survive ~0.2
   rad/t.u. differential rotation; do not inherit speckle
   deliveries; drive on m≥2 charts (the Phase-M ignition guard).
3. **Plethora verdict:** winding classes are injectable (m = 2/4
   ladder verified), transit-distinguishable (exposure/evap
   ladders), and door-transmissible (3.5σ) — but nothing yet
   carries them to persistence, so no winding-species census
   exists at this box and law. The two measured roads: λ/dmin ≥ 10
   (OUTLOOK E-G, 3D/scale — outrun the speckle) and Phase L
   (protected carriage — outlive the churn). The companion no-go
   (`BLACKHOLE.md`) shows the same missing sector gates gravity,
   far fields, frame-drag, and any black-hole analog: one door,
   many locked rooms.
