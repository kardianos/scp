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

## 5. Verdict

(after all arms; decisions are the user's)
