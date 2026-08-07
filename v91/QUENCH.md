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
