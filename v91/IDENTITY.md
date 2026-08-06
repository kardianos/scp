# IDENTITY — the parcel-carried ontological-identity lane (user-opened 2026-08-06)

User directive (2026-08-06, verbatim): *"first regrade B4, then open
parcel-carried ontological-identity lane."* The B4 re-grade is
executed (v90/REALITY.md §B4, commit before this file). This document
is the lane's pre-registration and record, written and committed
BEFORE any kernel edit and any run. Kernel work is explicitly
authorized by the directive; the ratchet governs: every knob defaults
inert (= the committed pre-lane kernel byte-exactly), the full battery
must read ALL GREEN 93 at defaults before any results commit, C and Go
kernels mirrored operation-for-operation, arms registered before they
run, bars sharpen and never soften. Nothing here adopts a law;
laws_V2g stays the table.

## 0. The thesis (inherited, now measured from three sides)

REGISTRY §5 closed every statistical alternative: identity inferred
from rates has no slot-grain gap to hide in, ignites the bath when
made a law (three measured ignitions), and is impotent as maintenance
(152 vs 140). ASTRO §6/COMBINE sharpened the demand: config-level
assertion (frozen gauges) keeps exactly one object class alive,
preserves skeletons but never contents (26/48 gauges rock-stable
through a starvation), and is null on the flagship. The registered
conclusion, now to be built and measured: **identity must be
ONTOLOGICAL — born with the matter, carried by its exchanges, dying
with it — never inferred from statistics.** Concretely (REGISTRY
§5.5, verbatim design): a gid stamped at the door when a voice's atom
fires, carried through flight, registered at arrival; born at
nucleation, dying with the cell (no permanent index, no background);
sgg growth gated on gid-continuity. Acceptance target, measured in
advance: **match or beat ret 0.5 @ t=5000 on the UUD chord with
self-carried identity, bath dark, battery green.**

Scope note from the kernel's actual transport structure (read before
design): the DENSE channel carries flight on slots (`slem[slot][dir]`
— discrete, directional, flushed-to-source on slot death) — exact
parcel carriage is representable there with no assignment rule. The
FIELD channel is a per-cell wave amplitude (`fa1/fa2`, `field_inject`)
— it has no discrete carrier; field-side identity is the S2
"amplitude mode label" formulation and stays DESIGN-ONLY this
campaign (same status REGISTRY §3.3 gave parcels). Consequence worth
registering: dense-flight identity runs with or without radiance —
identity carriage is radiance-independent even though the field
door's metabolism is not (COMBINE CO-AR).

## 1. The design

### 1.1 State (all of it; where it lives; no background)

| field | size | lives on | born | dies |
|---|---|---|---|---|
| `cgid_[i]` | NC | the cell | minted when x crosses par_hi upward (and at seed where x ≥ par_hi) | set 0 when x falls below par_lo |
| `gid_next` | 1 | the run | 1 at init | with the run (a label spring — nothing anywhere is indexed BY gid) |
| `cbirth_[i]` | NC | the cell | sim_t at mint | with the episode |
| `slgid_[2s+dir]` | 2·NLMAX | the slot | stamped at each dense deposit = depositor's cgid | with the slot |
| `sborn_[s]` | NLMAX | the slot | sim_t at `slot_new` | with the slot |
| `pstampa_[s], pstampb_[s]` | NLMAX | the slot | stamped at first mutual identity-carried delivery | cleared on endpoint-gid change/death; with the slot |
| `pstampt_[s]` | NLMAX | the slot | sim_t at stamping | with the stamp |
| `pdp_[s], pdm_[s]` | NLMAX | the slot | 0 at `slot_new` | with the slot |
| `pdel_[2s+dir]` | 2·NLMAX | the slot (per-step scratch) | 0 | consumed each step |

`cgid_` = the matter-EPISODE identity: x = (Em+flload)/cap crossing
par_hi upward mints a fresh gid; falling below par_lo retires it. The
identity attaches to the episode of being matter, not to the array
address — a re-filled cell is NEW matter with a new gid. (Measured
basis for the thresholds: drained voices touch x = 0.000–0.001 in the
standing runs — a fully drained voice IS dead matter under ASTRO
6.5.6's matter/topology split; lean-but-alive voices run x ≈ 0.01–0.15.)
`slgid_` = the parcel label: the depositor's episode-gid, stamped at
departure, read at arrival (design decision I-D1 below). `pdp_/pdm_` =
low-passed delivered rate per direction counting ONLY identity-carried
energy (parcel label == the source endpoint's CURRENT gid — stale
flight from a dead episode does not count). `pstamp*` = the bond's
identity pair, recorded once and never inferred. No energy is stored,
moved, or created by any of this; rule-α flushes are returns, not
deliveries, and are not stamped (same convention as the registry).

### 1.2 Keys (printed on their own header line, pinned by battery purity at inert defaults)

* `par_tau` (default **0** = lane OFF: no gid state, no stamps, no
  prints; the step is byte-identical to the pre-lane kernel). When
  > 0: the identity-ledger memory time-constant.
* `par_gate` (default **0**). When 1: the cantus gauge growth target
  is gated by the identity-continuity r_id (form at I-G1).
  `par_gate=1` without `par_tau>0` is a config error (r_id ≡ 0,
  gauge dark everywhere; harmless, printed).
* `par_form` (default **0**): 0 = I-A binary continuity; 1 = I-B
  flow-shaped (menu below).
* `par_lo` (default **0.002**), `par_hi` (default **0.02**): episode
  thresholds (hysteresis pair; inert when par_tau=0).
* `par_mature` (default **400**): the stamp maturity time — r_id
  scales by min(1, stamp_age/par_mature). The anti-ignition clock
  (closed-form case in §2).

### 1.3 Kernel sites (C first, Go mirrored operation-for-operation)

1. **Episode update** (once per cell per step, after loads are
   current): mint/retire `cgid_` by the hysteresis pair. Mint =
   `cgid_[i] = gid_next++; cbirth_[i] = sim_t`.
2. **Deposit stamp** (the dense departure site where slot in-flight
   gains energy): `slgid_[2s+dir] = cgid_[src]` (last-depositor
   label, I-D1).
3. **Arrival registration** (the delivery site, `take > 0` — the same
   site the registry stamps): if `slgid_[2s+dir] != 0` and equals the
   CURRENT `cgid_` of the source endpoint, `pdel_[2s+dir] += take`.
   Deliveries only; flushes untouched; atoms/credit machinery
   untouched.
4. **Pass H parcel block** (before the cantus block, all non-FREE
   slots): `kpar = min(1, dt/par_tau)`; low-pass `pdp_/pdm_` from
   `pdel_/dt`; zero the scratch. Stamp logic: an UNSTAMPED slot whose
   two endpoints hold live gids and whose pdp_ and pdm_ are both > 0
   stamps `(pstampa,pstampb) = (cgid_[sli], cgid_[slj])`,
   `pstampt = sim_t`. A STAMPED slot clears its stamp (and nothing
   else) the step either endpoint's cgid no longer equals its stamp.
5. **The continuity variable** (computed where used):
   `cont = stamped && cgid_[sli]==pstampa && cgid_[slj]==pstampb && pdp_>0 && pdm_>0`;
   `age = sim_t − pstampt`; `mat = min(1, age/par_mature)`;
   **I-A:** `r_id = cont ? mat : 0` ·
   **I-B:** `r_id = cont ? mat · 2·min(pdp_,pdm_)/(pdp_+pdm_) : 0`.
6. **The gate** (`par_gate=1`, inside the existing cantus growth,
   same action point as reg_gate — I-G3 carries R-G3's registered
   choice verbatim): the gate multiplies the GROWTH TARGET
   (`sgg += ktau·(gg·r_id − sgg)`), never the decay, never the
   lock/tune forces. `reg_gate` and `par_gate` both set is a config
   error (printed; par wins nothing — the run refuses law arms with
   two gates).
7. **Diagnostics** (only when `par_tau > 0`): `# PAR t= ...` at diag
   cadence — live-gid count, mint/retire rates, episode-age
   quartiles, slot-age quartiles by class (tag-pair vs bath),
   stamped-slot count and stamp-age quartiles, identity-carried
   delivery fraction (Σpdel_matched/Σdelivered), and for tri arms the
   triad's three gid ages; a final `# RESULT par ...`. QATOM lines
   gain a `gid=` field only when `par_tau > 0` (byte-purity at
   defaults preserved).

### 1.4 Decision points

* **I-D1 (parcel label mixing)** — in-flight can accumulate deposits
  across steps before arrival; the label is LAST-DEPOSITOR. The meter
  measures the approximation honestly (M-I4's identity-carried
  fraction; deliveries from a slot whose label went stale count as
  foreign). If the stale fraction exceeds 10% on bond-class slots,
  design round before any law arm.
* **I-G1 (form)** — I-A if it separates (bath r_id q90 < 0.2 with
  bond median > 0.6 in the meter arms), else I-B. Simplest wins; no
  new constant enters unless forced.
* **I-G2 (par_tau)** — menu {10, 30, 100}; carry REGISTRY R-G2's
  measured prior (τ=10 selected there); select by the widest
  bond-vs-bath separation of the pdp_/pdm_ ledgers.
* **I-G2b (par_mature)** — menu {200, 400, 800}, selected AGAINST
  M-I2's measured bath slot-age distribution: the smallest value
  ≥ the bath q90 slot age. If bath slot ages reach into the
  thousands, the maturity clock cannot separate and the lane reports
  the failure honestly (no threshold fishing).
* **I-G4 (episode thresholds)** — if M-I1 shows episode flutter
  (mint+retire on the same cell > 1/100 t.u. sustained), widen the
  hysteresis ×2 and re-run the meter once; both runs reported.

## 2. Closed-form expectations (written before any run)

1. **The naive gate is the FOURTH IGNITION — derived, not feared.**
   From the standing numbers: ~450k slot birth/death events over
   5000 t.u. across ~26k live slots → mean bath slot life ≈ 300 t.u.
   Bath episodes (cell x rarely crosses 0.002 downward in a warm
   bath) are long. So WITHOUT the maturity clock, a bath slot holds
   cont=1 for ~300 t.u. ≫ the cant_tau=50 arming time — gid-continuity
   alone would ignite the bath exactly as flow-reciprocity did. The
   maturity clock is the structural difference: **before arming the
   gauge exerts no force, so nothing the gauge field does can extend
   a bath slot's life to par_mature** — the lock cannot manufacture
   stamp AGE. Prediction: L-I2 bath stays dark at par_mature ≥ bath
   q90 slot age; if it ignites anyway, the thesis takes real damage
   (report as such).
2. **Maintenance can mature; self-assembly deadlocks.** Seeded-at-birth
   gauges (cant_seed=1, an instrument, not an assertion of identity —
   the identity is carried by the matter) hold the bond alive while
   the stamp matures honestly → identity-gated maintenance should
   hold the chord where flow-gated maintenance failed (152 vs 140).
   Self-assembly (cant_grow=1, no seed) cannot arm: an unlocked chord
   edge dies at t≈36–48 (COMBINE), far below par_mature — the bond
   cannot survive to maturity without the lock it is trying to grow.
   Registered as the lane's measured design tension; if L-I3
   self-arms anyway, the deadlock prediction is falsified — either
   way a result.
3. **The triad's gids should be the programme's first "same matter"
   statement**: three gids minted at seed, alive at t=5000, carrying
   the whole metabolism through their door events. Cool/lean arms
   sit near the retire threshold — the U voices of the cool chord
   (x ≈ 0.008–0.010) are ABOVE par_lo=0.002 but a drained voice
   (x → 0.000–0.001) correctly loses its identity: matter that
   empties is dead; what refills is new.
4. **Identity-carried delivery fraction ≈ 1 on bonds.** Bond flight
   times are ~d/C ≈ 1.5 t.u.; episodes live thousands — stale labels
   should be rare on bonds (M-I4 validates I-D1). Bath somewhat
   lower (churn).

## 3. Pre-registered arms and bars

### 3.1 Verification bars (before any physics arm)

* **P-ID0 (byte purity):** at `par_tau=0` the new binary's full diag
  stream is byte-identical to the committed pre-lane kernel on the
  reference configs (chord recipe T=200; bath L=16 T=100), −O0 both.
  Battery ALL GREEN 93 at defaults. Header purity line carries
  `par_tau=0 par_gate=0 par_form=0 par_lo=0.002 par_hi=0.02 par_mature=400`.
* **P-ID1 (C≡Go):** with the meter firing (`par_tau=10`, chord recipe
  T=200), C and Go diag rows byte-equal (same discipline as the
  registry apparatus).
* **P-ID2 (meter physics-silence):** `par_tau=10` vs `par_tau=0` on
  the chord recipe T=1500: physics rows byte-identical once `# PAR`
  and `gid=` fields are stripped. The meter must not touch physics.

### 3.2 Meter arms (par_gate=0; the standing recipes, one knob added)

* **M-I1 (bath demographics):** warm bath L=16 T=1500 par_tau=10.
  Report: live-gid fraction, mint/retire per t.u., episode-age
  quartiles. Bar: none (measurement); I-G4 trigger as registered.
* **M-I2 (slot ages):** same run — slot-age quartiles at death,
  bath-class vs the chord arm's bond-class. Feeds I-G2b.
* **M-I3 (the triad's identity):** the frozen chord reference recipe
  + par_tau=10, T=5000. Bars: the three triad gids minted at seed
  are ALIVE at t=5000 (first "same matter" statement); bond stamps
  mature (age@5000 ≥ 4000); identity-carried delivery fraction on
  the three bonds ≥ 0.9 (I-D1 validation, else design round).
* **M-I4 (purity under lean):** the cool chord recipe + par_tau=10,
  T=5000. Bars: U-voice gids survive the starved-lean mode (x ≈
  0.008–0.010 > par_lo); delivery fraction reported. If lean voices
  flutter their episodes, I-G4 fires.

### 3.3 Law arms (registered now; RUN ONLY after a design round
publishes I-G1/I-G2/I-G2b selections against the meter data)

* **L-I1 (identity-gated maintenance — the acceptance arm):** UUD
  chord, cant_seed=1, cant_tau=50 (honest decay, NOT frozen),
  k_cant=1 k_tune=0.2, par_tau=I-G2, par_gate=1, par_form=I-G1,
  T=5000, seed 20260802 (+314159 replicate). **Acceptance bars
  (REGISTRY §5.5 verbatim): ret@1500 ∈ [0.36,0.66] and ret@5000 ≥
  0.35 with NO frozen assertion — matching the frozen bound by
  SELF-carried identity; bath dark (a_max off-object < 0.5, nlock =
  the object's own slots only); battery green.** References: frozen
  0.5101/0.4975; flow-maintained FAILED at 152-vs-140 grade; ctl
  dies.
* **L-I2 (flammability null):** bath-only (no object), cant_grow=1,
  par_gate=1, T=1500, warm. Bars: nlock=0 (or ≤ 3 transient), a_max
  < 0.5, conversion economy within 3% of the no-cant control — the
  bath CANNOT manufacture stamp age (closed-form §2.1). This is the
  first self-growing form in four rounds predicted NOT to ignite.
* **L-I3 (self-assembly probe — deadlock on record):** tri seed,
  cant_grow=1, cant_seed=0, par_gate=1, T=1500. Registered
  prediction (§2.2): the bond dies before its stamp matures — NO
  self-arming. If tri slots arm and the chord survives ≥ 1500, the
  deadlock is falsified and the lane has its miracle; report either
  way.

### 3.4 Out of scope this campaign

Field-channel (mode-label) identity — design-only, S2's lane. Any
law-table adoption. Any change to atoms/door/space/contact machinery
beyond the registered sites. Winding-sector bookkeeping. Escalation
(per REGISTRY: only if the meter arms strand the thesis, and then
design-only first).

## 4. Results

**4.0 Verification (2026-08-06; apparatus landed both kernels).**

* **P-ID0 PASS.** At defaults the new binary's full output is
  byte-identical to the committed pre-lane kernel on both reference
  configs (chord recipe T=200; bath L=16 T=100), −O0 both, drift
  column included — only the new header line differs. Battery **ALL
  GREEN (93 bars)** with the identity purity line pinned
  (`runs/BATTERY_identity.log`).
* **P-ID2 PASS, both kernels.** C at T=1500 and Go at T=200: with
  `par_tau=10` the physics stream is byte-identical to `par_tau=0`
  once the meter's own lines (`# PAR`, `# RESULT par`, the QATOM
  `gid=` field) are stripped. Trajectory-grain proof besides the
  byte proof: per kernel, event and channel-birth counts are
  bit-identical with the meter on vs off (C: 4426 QATOMs / 24646
  births both ways; Go: 4398 / 24557 both ways). The meter moves
  nothing.
* **P-ID1 — amended by measurement, cause demonstrated.** The
  registered "C and Go diag rows byte-equal" is NOT achievable on
  the warm-chord config at T=200 for ANY apparatus: the unmodified
  HEAD pair already diverges there (C 4426 vs Go 4398 QATOM events;
  births 24646 vs 24557; drift column differs from t=4 — a 1-ulp
  C/Go arithmetic difference that the chaotic medium amplifies into
  trajectory divergence within ~60–200 t.u.). Demonstrated directly
  on the pre-lane binaries; the lane's C↔Go deltas are IDENTICAL to
  the pre-lane deltas. What holds and is verified: the apparatus
  fields (gid stream, PAR lines) are byte-equal across kernels on
  the agreeing window, and each kernel is byte-silent within itself
  (P-ID2). Recorded as the pair's standing FP-divergence envelope —
  the bar's original formulation was impossible, not failed.
* Meter first-light facts (shakedown T=200, warm chord recipe): the
  triad mints gids **1, 2, 3** at birth and holds them; bath live
  fraction 93%, mint/retire ≈ 16/3.5 per t.u. at warm-up; slot
  deaths run at median 2–3 t.u. with q90 ≈ 150 (the §2.1 mean ≈ 300
  was tail-dominated — churn is fiercer than the mean); ~20k slots
  stamp easily (as expected); identity-carried delivery fraction
  0.961 at T=200 rising to **0.993** by T=1500 — the I-D1
  last-depositor label is near-exact at steady state. Note for the
  design round: bath STAMP survival is endpoint-episode-churn
  limited, so I-G2b must select par_mature against the measured
  stamp-age distribution, not slot ages alone.

**4.1 M-I1 + M-I2 (warm bath demographics; stationary by t≈600;
log `runs/identity/mi1_bath_warm.log`).** Live episodes 2555–2585 of
2703 cells (95%); churn balanced at 4.85 mints ≈ 4.83 retires per
t.u. → per-cell episode hazard 0.00188/t.u., mean episode life ≈
530 t.u.; live-age census q25/50/75 = 170/440/941. No episode
flutter — I-G4 does not fire. Slot deaths: median ≈ 2 t.u., q90 ≈
100–140 (the churny fringe; the standing core is long-lived).
Identity-carried delivery fraction 0.9936 (I-D1 near-exact).
**The decisive measurement — bath STAMP ages are memoryless with
hazard = 2 × episode hazard:** predicted median 261·ln 2 = 184 vs
measured 181 (and ~21,600 standing stamps). Consequence, computed
BEFORE any law arm: P(stamp age > par_mature) = exp(−M/261) → at
the registered M=400, **~22% of bath stamps (~4,700 slots) are
mature at any instant**; at the menu's top M=800, ~4.6% (~1,000).
The §2.1 closed-form claim ("the maturity clock closes ignition
structurally") is REFUTED by measurement: a warm bath whose
episodes live 530 t.u. manufactures stamp age simply by persisting.
The registered honest-failure clause (I-G2b) is armed; the design
round (§4.3) must reason from the measured hazard, not the slot-age
rule.

**4.2 M-I3 + M-I4 (the triad at identity grain; logs
`runs/identity/mi3_chord_warm.log`, `mi4_chord_cool.log`).**

* **The first complete "same matter" statement — and it localizes
  in the D voice.** Warm: the D engine holds **gid 3 from seed to
  t=5000** (age 5000; ret 0.4975 byte-matches the ASTRO reference —
  the meter touched nothing). Cool: same — gid 3, age 5000, in the
  ring-killing medium. The registered M-I3 bar (triad gids alive at
  5000) is MET by the D voice and FAILED by the U voices, and the
  failure is the discovery:
* **U matter is replaceable fuel.** Warm U voices end as gids
  21472/21759 (ages 780/722) — each U episode died and was reborn
  ~6 times over the run as its load dipped below par_lo during deep
  lean excursions (the t2far/t2near x→0.000–0.001 dips, now read at
  identity grain). Cool: one U identity-dead at close (x = 0.004 —
  inside the hysteresis dead band, exactly as designed), one a
  35-t.u. rebirth. **The flux machine is ONE persistent identity
  plus consumable matter at the intake positions.** The M-I4 bar
  ("U gids survive lean") FAILS as registered — honestly: matter
  that empties IS dead under the lane's own ontology; par_lo stays
  0.002 (no threshold fishing to make empty matter carry identity).
* Bond-grain: tagged-pair slot deaths over 5000 = 4 (edge blinks,
  ages 0.1–11.5); 2 of 3 bond stamps current at close (U rebirths
  clear stamps; re-stamp follows re-mint within ~par_tau); stamp
  max-age 4101 (a bond stamp carried >4000 t.u.). idfrac
  0.9959 warm / 0.9364 cool — I-D1 bar (≥0.9) PASS everywhere.

**4.3 DESIGN ROUND (registered 2026-08-06, after the meter arms,
BEFORE any law arm — selections argued from measured numbers).**

* **I-G1 = I-A** (binary continuity). Bond pdp/pdm are robust and
  the flow-shaped form adds a constant with no measured need.
* **I-G2 = par_tau 10.** Re-stamp latency after a U rebirth ≈
  par_tau; smallest menu value = fastest honest recovery, and the
  REGISTRY prior (τ=10 selected) carries.
* **I-G2b — the maturity clock SPLITS by arm class (the measured
  stamp hazard forces it):**
  - **L-I1 (maintenance; cant_grow=0, gauges seeded, nothing else
    ever grows): par_mature=0 — no ramp.** The clock's
    anti-ignition job does not exist in this arm (bath gauges are
    never seeded and cannot grow), and any ramp would convert every
    U-rebirth stamp reset into a kill window. The gate is bare
    identity continuity: the thing REGISTRY measured the lock
    cannot manufacture.
  - **L-I2 / L-I3 (self-growing): par_mature=800 (menu top), with
    the REVISED prediction on record:** the measured bath stamp
    hazard (memoryless, mean 261 t.u.) leaves ~4.6% of ~21.6k bath
    stamps mature at any instant (~1,000 slots) — **ignition is
    now EXPECTED**, not excluded; §2.1's structural claim is
    already refuted by the meter, and the arm tests the revised
    arithmetic. The design-only escalation if L-I2 ignites as
    revised-predicted: joint maturity on min(stamp age, both
    endpoint episode ages) at ~2400 t.u. (bath joint fraction
    ~1e-4) — outside the registered menu, therefore NOT run this
    campaign; recorded for the user.
* **L-I1 expectation, quantified from the meter:** during a U dead
  window the bond's gate target is 0 and the seeded gauge decays as
  e^(−Δ/50); the U dead windows are the lean-floor excursions (rare
  — ~6 per U per 5000 t.u.; duration set by refill from the D's
  shine). If dead windows run ≳50 t.u. the gauge dips below hold
  and the bond faces the unlocked death clock (~36–48 t.u.) — the
  arm measures exactly whether identity-gated maintenance rides
  through its own honest gaps. Acceptance unchanged (§3.3):
  ret@1500 ∈ [0.36,0.66], ret@5000 ≥ 0.35, bath dark, battery
  green.

**4.4 L-I2 + L-I3 (self-growing arms; logs
`runs/identity/li2_flam.log`, `li3_selfasm.log`).**

* **L-I2 — THE FOURTH IGNITION, measured as revised-predicted.**
  nlock = **2595** at t=1500 (bar ≤ 3), a_max 0.996, tune_total
  32,086 — the warm bath crystallizes under the identity gate at
  par_mature=800 exactly as §4.1's stamp-hazard arithmetic said it
  would (~1,000 mature stamps to seed the cascade). Two causes, now
  separated: (i) mature bath stamps are ABUNDANT (the 261-t.u.
  memoryless hazard), and (ii) **the maturity RAMP leaks
  pre-maturity force** (r_id = age/M gives partial growth targets
  everywhere from t=0). The identity-continuity VARIABLE remains
  unmanufacturable — what ignited is the ramp's leak plus the
  bath's own longevity. Registered bar FAILS as revised-predicted;
  §2.1's original structural claim stands refuted twice (meter,
  then law arm).
* **L-I3 — the deadlock prediction is FALSIFIED, and the
  falsification is CONFOUNDED.** The unseeded chord self-armed
  (a_tag 0.97, two edges alive) and lived to 1500 (ret 0.3740,
  conn 1.0) — but inside a box that had itself ignited (nlock
  2598). The mechanism is the same ramp leak: pre-maturity partial
  force bridged the t≈36–48 unlocked death window, for the chord
  AND for the bath. A self-armed chord in a crystallized medium is
  not self-assembly by the programme's honest-medium standard
  (CANTUS/REGISTRY). Recorded as: ramped maturity = a self-assembly
  bridge and an ignition leak — the SAME mechanism, inseparable at
  this design point.
* **Design-only escalation, recorded for the user (not run —
  outside the registered menu):** a STEP maturity (no force below
  M) on the JOINT clock min(stamp age, both endpoint episode ages)
  at M ≈ 2400 would close the leak (bath joint fraction ~1e-4) —
  and with the leak closed, self-assembly re-deadlocks (no bridge).
  The coherent reading of all four self-growing rounds across three
  lanes: **identity preserves; it does not create.** Genesis of
  coherence needs a different door (S2 §2 / the field-side mode
  label), not a cleverer gate.

**4.5 L-I1 (identity-gated maintenance — the acceptance arms; logs
`runs/identity/li1_maint.log`, `li1_maint_s314159.log`).**

**Acceptance FAILS as registered, both seeds, both windows:**
ret@1500 = 0.2364 / 0.2656 (bar [0.36,0.66]); ret@5000 = 0.2739 /
0.1831 (bar ≥ 0.35). Self-carried identity at honest decay
(cant_tau=50, no frozen assertion) does NOT match the frozen bound
(0.4975 / 0.4114).

What it measurably DOES — a new state class, reported precisely:
  - **Alive at t=5000, without assertion, in both seeds**: bonds
    locked at close (nlock=2, a_tag 0.65/0.73), conn 0.667/1.0,
    matter oscillating 0.09–0.47 through deep dips WITH RECOVERIES
    (0.099@1000 → 0.397@2000; 0.093@4500 → 0.183@5000). The
    references: flow-gated maintenance faded to death-grade (152 vs
    140); no lock dies by t=560. **Identity-gating is the first
    maintenance that keeps the chord alive to 5000 with honest
    decay.** The bath stayed dark by construction (nlock = the
    object's own slots only).
  - **The shortfall's mechanism, at identity grain:** every fuel
    (U) episode death opens a stamp gap; each gap is a gauge-decay
    window (target 0 for the dead span + re-stamp latency ≈
    par_tau); the chord survives its gaps but bleeds retention
    through them — maintenance inherits the fuel's mortality. Seed
    20260802's engine held **gid 3 for the full 5000** (same
    matter, maintained, unfrozen — the statement the lane was built
    to make); seed 314159's rough state killed even the engine's
    episode (D gid dead at close, age 0) — under honest
    maintenance, nothing is guaranteed the engine's persistence
    either.

## 5. Verdict (all arms complete; drift ≤ 2e-15; DECISIONS ARE THE USER'S)

1. **The lane exists, exact, inert, verified.** Episode gids +
   parcel labels + identity stamps in both kernels behind
   par_tau/par_gate/par_form/par_lo/par_hi/par_mature; defaults
   byte-identical to the pre-lane kernel (−O0, drift included);
   meter physics-silent at byte grain in both kernels; battery ALL
   GREEN 93 with the purity line pinned. The C↔Go pair's
   pre-existing FP divergence envelope is now measured and recorded
   (P-ID1 amended, cause demonstrated on HEAD binaries).
2. **Identity is cheap and near-exact to carry.** The last-depositor
   parcel label is 99.4–99.6% exact at steady state; the whole
   apparatus is bookkeeping riding on deliveries.
3. **The first "same matter" statements — and matter has a class
   structure.** The chord's D engine keeps one gid for its whole
   life (warm, cool, and under honest maintenance in one of two
   seeds); its U intake positions are REPLACEABLE FUEL (episodes
   die and refill ~6×/5000 t.u.). The warm bath is 95% live
   episodes with mean life 530 t.u. — a medium of short-lived
   identities.
4. **Identity-continuity is unmanufacturable — and still not a
   self-assembly gate.** The binary continuity variable cannot be
   forged by the lock (gids are birth events). But the maturity
   RAMP leaks pre-maturity force: fourth ignition (nlock 2595, as
   revised-predicted from the measured 261-t.u. stamp hazard), and
   the same leak bridges L-I3's deadlock (confounded
   self-assembly in a crystallized box). Step-maturity on the joint
   episode+stamp clock (~2400) would close the leak and re-deadlock
   self-assembly — design-only, recorded. **Across four
   self-growing rounds in three lanes: identity preserves; it does
   not create.**
5. **The acceptance target is NOT met — and the miss is
   informative.** Identity-gated maintenance holds LIFE to 5000
   without assertion (first ever) at roughly half the frozen
   retention, failing the registered bars because maintenance
   inherits fuel mortality. The frozen bound remains
   assertion-dependent. What the whole arc now says: config
   assertion holds shape; ontological identity preserves engines
   and enables assertion-free survival; NEITHER reaches the frozen
   bound honestly, and nothing creates. The remaining registered
   door is field-side identity — the S2 amplitude mode label —
   where the parcel is the wave and the identity gap cannot open.
6. **Nothing was adopted.** All knobs default inert; laws_V2g
   verbatim; the battery gates; every arm registered before it ran;
   the two §2 closed-form errors (ignition arithmetic; the deadlock)
   are corrected by measurement in place, not rewritten.
