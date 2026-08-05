# REGISTRY — the exchange-registry lane: identity-carrying transfers

Opened 2026-08-05, **user-directed**: "Use the CANTUS experiment as a
mechanism of discovery. Ensure the discoveries from that are
documented and other concepts clearly adjusted inline. **Do open
exchange registry lane and do go for a UUD chord.**" This authorizes
the registry kernel candidate (inert by default, ratchet governing)
and names the goal object. Predecessors: `CANTUS.md` (the measured
need), `carried/S2_CHANNEL.md` §2 item 3 + addendum (the design slot),
`THEORY.md` §2.2 (the standing constraints).

## 0. The thesis

CANTUS measured that coherence cannot self-assemble on closure
statistics in this substrate: the churn bath is a persistent-closure
medium (66% of cells hold time-averaged two-sided gate > 0.5 at zero
force), so every closure-gated growth rule ignites a bath-wide
Kuramoto transition. The discriminator between a BOND and CHURN must
be **identity**. The registry thesis, one sentence:

> **A bond is a remembered reciprocal exchange between two continuous
> identities — and coherence may only grow there.**

The minimal realization needs no new tags at all, because the
substrate already carries the identity structure: a **slot** lives
exactly as long as both endpoints exist and stay in contact
(`slot_new` resets everything at rebirth — the identity relationship
is BORN and DIES with the slot; no background). Every pass-5 delivery
across slot s is therefore a parcel from a known continuous identity
to its known continuous partner. "Identity-carrying transfers" at
this grain = **a per-slot ledger of reciprocal deliveries**: the
registry remembers, per link, whether energy has actually been
exchanged BOTH WAYS, recently, between the same two voices.

Why this can discriminate where closure could not (the prediction
this campaign exists to measure, not assume): bath closure is phase
alignment and is free (66%); but bath TRANSPORT is conductive — energy
transits a churn slot predominantly one way (down the local gradient)
or in isolated nucleation/collapse bursts, while a bond at a rung is
a STANDING two-way exchange (the flux-machine interior: CQ7/CQ8
forced persistent internal circulation; the nv=48 object runs books
0.87–0.96 at 100–200× the sealed-door metabolism). Phase is cheap;
reciprocal standing exchange with one continuous partner is exactly
what "bond" has meant since v89.

**Escalation branch (registered now):** if the meter finds the bath
is ALSO reciprocity-rich at slot grain (no usable gap), then
slot-grain identity is insufficient and the design escalates to
parcel-carried identity (atoms stamped at fire with a birth-counter
gid, registered at the receiving door) — a substantially bigger law
change that will NOT be implemented in this campaign without that
measured need.

## 1. The design

### 1.1 State (all of it; where it lives; no background)

| field | size | lives on | born | dies |
|---|---|---|---|---|
| `rfp_[s]` | NLMAX | the slot | 0 at `slot_new` | with the slot |
| `rfm_[s]` | NLMAX | the slot | 0 at `slot_new` | with the slot |
| `rdel_[2s+dir]` | 2·NLMAX | the slot (per-step scratch) | 0 at `slot_new` | consumed each step |

`rfp_` = low-passed delivered rate in direction i→j (slot dir 0);
`rfm_` = j→i (dir 1). Units: energy per t.u. No energy is stored,
moved, or created by the registry — it is bookkeeping riding on
deliveries, like the credit register.

### 1.2 Keys (new law-candidate keys, printed on their own header line, pinned by battery purity at inert defaults)

* `reg_tau` (default **0** = registry OFF: no stamps, no prints, the
  step is byte-identical to the pre-registry kernel). When > 0: the
  memory time-constant of the delivery ledger.
* `reg_gate` (default **0**). When 1: the cantus gauge growth target
  is gated by the registry match r_s (form selected at R-G1 below).
  `reg_gate=1` without `reg_tau>0` is a config error (r_s ≡ 0 —
  the gauge field stays dark everywhere; harmless, printed).

### 1.3 Kernel sites (C first, Go mirrored operation-for-operation)

1. **Pass 5 stamp** (the delivery site, `take > 0`): when the
   registry is on, `rdel_[2s+dir] += take`. Deliveries ONLY — the
   rule-α flush is a RETURN to source, not an exchange, and is not
   stamped. β-credit/atoms machinery untouched.
2. **Pass H registry block** (runs BEFORE the cantus block, s
   ascending 0..NSLOT, all non-FREE slots — including pinched
   `sA=0` slots, so the memory decays through lens blinks instead of
   freezing):
   `kreg = min(1, dt/reg_tau); rfp_ += kreg·(rdel_[2s]/dt − rfp_);
   rfm_ += kreg·(rdel_[2s+1]/dt − rfm_); rdel_[2s]=rdel_[2s+1]=0.`
3. **The match** (computed where used):
   `gross = rfp_+rfm_; rho = gross > 0 ? 2·min(rfp_,rfm_)/gross : 0`
   ∈ [0,1]. 1 = balanced standing exchange, 0 = one-way or silent.
4. **The gate** (`reg_gate=1`, inside the existing cantus growth):
   `sgg += ktau·(gg·r_s − sgg)` where r_s is the R-G1-selected form.
   With `cant_grow=1` this is the identity-gated SELF-ARMING law:
   bonds ignite their own gauges through their own reciprocal
   exchange; the bath stays dark **by construction** (r ≈ 0), not by
   threshold. `cant_seed` remains available for instrument arms only.
5. **Diagnostics**: `# REG t= ...` at diag cadence and a final
   `# RESULT reg ...` — per class (tag-pair slots = both endpoints
   tagged; bath slots = neither), n, ρ quartiles, gross median,
   fraction-with-flow. Prints only when `reg_tau > 0`.

### 1.4 Decision points

* **R-G1 (match form)** — selected by the meter arms, menu registered
  now: **F-A** r = ρ; **F-B** r = ρ·gross/(gross+f0) with f0 = the
  geometric midpoint of the measured bond/bath gross medians (guards
  against two-microscopic-parcel ρ=1 artifacts); **F-C** r = ρ².
  Selection rule: the simplest form whose bath tail stays below 0.2
  while the bond median stays above 0.6 in the meter arms. If F-A
  satisfies it, F-A wins (no new constant enters the law).
* **R-G2 (reg_tau value)** — meter arms run at reg_tau ∈ {10, 30,
  100}; select the value with the widest (bond median − bath q90) ρ
  separation. Physical prior: a bond parcel period is d/C ≈ 1.2–1.9
  t.u. (C=1), so τ=30 spans ~15–25 delivery periods; τ must also sit
  well below object lifetimes (control UUD t_half 140).
* **R-G3 (what the gate multiplies)** — registered choice: the gate
  multiplies the GROWTH TARGET (gg·r), not the decay, and not the
  lock/tune forces directly. A bond that stops exchanging loses its
  gauge at 1/cant_tau — coherence is maintained by exchange, which is
  the thesis. (Gating the forces directly would let armed gauges act
  while silent — rejected at design time.)

## 2. Closed-form expectations (written before any run)

* At a fed bond (pair at rest on its rung): deposits flow both
  directions every step; deliveries land every ~d/C ≈ 1.2–1.9 t.u.
  per direction ⇒ rfp ≈ rfm ⇒ ρ → near 1 on a τ ≫ 2 t.u. ledger.
* At a conductive churn slot in a cascade: flow direction follows the
  local gradient; if direction persistence exceeds τ, ρ ≈ 0; if the
  foam flips direction faster than τ, ρ rises — **this is the honest
  risk of the design**, and reg_tau selection (R-G2) exists exactly
  to find the window between bond-breathing (fast) and gradient-flip
  (slow). If no window exists at any tested τ, P-REG1 fails and the
  escalation branch triggers.
* The registry meter at reg_gate=0 moves NO physics — every non-`#`
  output line byte-identical to the same run without reg_tau.
* At defaults (reg_tau=0) the kernel is byte-identical to the
  committed pre-registry kernel, binary-diff level.

## 3. Pre-registered predictions, arms, and bars

Ambient/harness carried from CANTUS: k005 radiance stack
(`k_rad=0.05 p_rad=4 rad_clock=0`), bath=1 noise_amp=0.5 (0.15 for
ring6@0.47), seed=20260802, controls verifiably run.

### 3.1 Meter arms (reg_gate=0; no cantus force; k_cant=k_tune=0)

| arm | config core | measures |
|---|---|---|
| M1 bath null | bare bath, T=600 | bath slot ρ/gross distributions at τ ∈ {10,30,100} |
| M2 UUD | `exp=tri tri_xU=0.28 bath=1` T=600 | tag-pair bond ρ/gross vs its own bath |
| M3 pair | `exp=pair pair_x0/1=0.47 doff −0.15` T=600 | the two-voice ledger |
| M4 i5 ring | `exp=ring ring_n=6 ring_x=0.47 noise 0.15` T=600 | multi-edge ledger incl. blink behaviour |

* **P-REG1 (the gap — THE campaign hypothesis):** tagged-bond median
  ρ exceeds the bath q90 ρ at some tested τ, with bond median ≥ 0.6
  and bath median ≤ 0.3. PASS ⇒ R-G1/R-G2 select the form and τ;
  FAIL ⇒ escalation branch (parcel gids), no law arm this campaign.
* **P-REG2 (meter honesty):** M1 with reg_tau=30 vs reg_tau absent —
  all non-`#` lines byte-identical; battery at defaults ALL GREEN 93
  after the meter lands; C≡Go byte-equal with the meter ON.

### 3.2 Law arms (reg_gate=1, cant_grow=1, NO cant_seed; k_cant=1 k_tune=0.2 cant_tau=50 — the CANTUS-selected coupled point)

Run order is the registered order; a failed guard STOPS the lane
(record, no adoption path, no goal-arm spin).

* **P-REG3 (the vacuum guard — the bar that killed v1 and v1.1):**
  bare bath, FULL stack live, T=480: cond and rad within 10% of the
  k005 control (2203.97 / 2016.22); tune_total ≤ 1000 (v1.1 bath:
  13k–175k); nlock ≤ 2× the null-meter level. **The medium must stay
  a medium.**
* **P-REG4 (self-arming where it should):** in M2-class UUD bath at
  the same point, the three internal bonds arm themselves (sgg ≥ 0.5
  within ~5τ) through their own exchange — no seeding.
* **P-REG5 (THE GOAL — the UUD chord):** `exp=tri tri_xU=0.28
  bath=1` T=5000, full stack: **t_half ≥ 1400 (10× the k005 control
  140)**. References: the frozen-gauge family bound measured 1160
  (8.3×) WITHOUT re-arming; the self-arming law can re-arm reborn
  slots, so exceeding the bound is possible and would be the first
  measured benefit of identity over frozen memory. Report alongside:
  retention curve shape (the bound showed ret RECOVERING to 0.587),
  books (in/out), internal gauge levels, d vs the chord rung.
* **P-REG6 (secondaries):** pair arm ≥ 50% retention at T=1200 (the
  P-Hi4 bar it failed at 25%); i5 ring recorded (winding is a known
  separate wall — no bar, measurement only).
* **P-REG7 (ratchet):** battery ALL GREEN 93 at defaults after every
  kernel change; C≡Go byte-equal with the FULL stack firing on a
  bath arm; drift ≤ 1e-12 on every arm.

### 3.3 What is explicitly out of scope this campaign

Parcel-carried gids (unless P-REG1 fails — and then design-only),
winding-sector bookkeeping (registered CANTUS finding, own lane),
any table adoption (laws_V3 remains the user's decision), any change
to atoms/door/space/contact machinery (S2 §2 item 5 invariants).

### 3.4 Design round 2 (2026-08-05 — registered AFTER the meter arms M1–M4, BEFORE any law arm)

**P-REG1 verdict: FAIL as registered — and the failure is a
measurement, not a mishap.** At the warmed ledgers (stationary by
t≈300): bath ρ median 0.49 (τ=10) and 0.68 (τ=30), both ≫ the
registered 0.3 ceiling; bond medians clear 0.6 only marginally and
not uniformly (§4.1 table). Two discoveries force a redesign of the
match form, not of the thesis:

1. **The churn bath is a persistent-EXCHANGE medium.** 92–93% of
   live bath slots carry delivered flow within a τ=10–30 memory;
   slot turnover ~100+ t.u.; and the ρ distribution RISES with the
   ledger window (median 0.49 → 0.68 from τ=10 → 30): **per-line
   balance, absent instantaneously (the v90 result), EMERGES
   statistically at long windows.** The bath is a population of
   long-lived, low-current proto-bonds — closure (CANTUS) and now
   balance are both free there. Identity cannot be read from rate
   SYMMETRY.
2. **Balance is a unison instrument; the chord's signature is
   FLOW.** The resting unison pair reads ρ = 0.94–0.97 — but the
   flavored (chord) bonds read ρ 0.5–0.9 because they carry a
   standing directed circulation — the flux-machine interior current
   (CQ7/CQ8's forced conclusion), here OBSERVED PER LINK for the
   first time. Meanwhile the ALIVE bond's gross delivered flow is
   15–300× the bath median (pair 0.046/t.u. at ρ 0.94; i5 edges
   0.008–0.013 at ρ 0.90–0.97; UUD ~0.002 early-life; bath
   stationary 0.7–1.4e-4). The weak "2–5×" endpoint numbers first
   harvested were dying-object artifacts — the RESULT-line summary
   is survivorship-biased for mortal bodies; mid-life series are the
   honest bond statistics.

**Amended match menu** (implemented behind `reg_gate`, both kernels,
A/B byte-verified; `reg_f0` is a new printed key, purity-pinned at
0): **F-B** (`reg_gate=1`): r = ρ·gross/(gross+f0) — balance ×
flow-saturation; f0=0 reduces F-B to the original F-A exactly.
**F-D** (`reg_gate=2`): r = s/(s+f0) with s = 2·min(rfp,rfm) =
ρ·gross — the RECIPROCAL FLOW rate; a one-way conductor scores 0 at
any magnitude; a balanced trickle scores ~0; a sustained two-way
current saturates to 1. F-D is the thesis stated in flow units.

**Amended selector (replaces the statistical gap as the gate
decision):** the registered purpose of P-REG1 was ignition
protection. Ignition is a PHYSICAL event and the vacuum guard
P-REG3 measures it directly — a warm but sparse, non-percolating
tail of gate-passing bath slots is harmless, and only the guard can
tell. Therefore: **the P-REG3 guard becomes the form/f0 selector.**
Menu: {F-B, F-D} × f0 ∈ {2e-4, 5e-4} (the geometric midpoint of the
measured alive-bond vs bath reciprocal-flow scales, and its
conservative double), all at τ=10 (the least balance-contaminated
ledger; τ=100 registered expectation: bath ρ median ≥ τ=30's — if
it comes in LOWER the emergence story is wrong and R-G2 reopens).
Guard bars unchanged from P-REG3 + one addition: bath `a_max`
(# CANT) < 0.5 throughout. Among guard-passing combos, a short UUD
arming probe (T=300, self-arming, no seed) selects the combo with
the highest internal a_tag at t=200 ("arm the object, not the
world"); tie → F-D (fewer factors).

**Instrument lessons (this round):** (i) a missing RESULT line means
UNFINISHED, not dead — verify a run's exit by PID before requeueing;
I misread six in-flight arms as killed and double-launched second
writers onto their log paths (two processes truncating/interleaving
one file = both datasets lost; the six arms were stopped and rerun
cleanly, single-writer). Corollary verified in passing: `make` over
the running binary is BENIGN on this platform (ld replaces via
unlink+rename; the six original arms sailed through a rebuild
untouched). (ii) RESULT-line summaries are survivorship-biased for
mortal objects — read mid-life series.

## 4. Results

(§4 is written only from executed runs; nothing above this line may
be edited after the first run — corrections land as dated addenda.)

### 3.5 Design round 3 (2026-08-05 — registered AFTER the §4.2 guard verdict, BEFORE the goal arms)

**The guard rejected the entire §3.4 menu** (all four combos ignite
the bath — §4.2). The mechanism, read off the logs: **the lock
manufactures its own gate evidence.** Locking a neighborhood
equalizes loads and pitches → opens gates → dense wants flow BOTH
ways between equalized neighbors → reciprocal flow r rises → gauges
grow → the lock spreads. Closure (CANTUS v1.1) and reciprocal flow
(this round) are both DYNAMICAL quantities the lock improves; any
slot-grain statistical gate built from such a quantity is flammable.
**The self-assembly door therefore needs a gate variable the lock
cannot manufacture** — identity carried from BIRTH (parcel gids
stamped at fire, registered at the receiving door): the §0
escalation branch is now the measured design requirement, recorded
design-only per §3.3.

**The goal arms restructure to the honest-instrument frame** (the
CANTUS v1.1i shape — sanctioned in §3.2: seeding is instrument-only),
which is immune by construction: `cant_grow=0` means NOTHING arms in
the bath, ever; seeded gauges on the object's own bonds are then
MAINTAINED by the registry (growth target gg·r on armed slots —
coherence sustained by identity-continuous exchange, fading when the
exchange stops). Form F-D, f0 = 2e-4, τ = 10 (the §3.4 tie rule and
the measured bond scale: alive UUD bonds s ≈ 5e-4–2e-3 ⇒ r ≈
0.7–0.9; a dying bond's s falls through f0 and its gauge honestly
fades).

Arms (all T=5000 except pair; k005 ambient; controls same T, same
binary): **uud_ctl** (no cantus); **uud_seed_reg** (cant_seed=1
cant_grow=0, k_cant=1 k_tune=0.2 cant_tau=50, reg_gate=2) — THE
measurement; **uud_seed_plain** (same, reg_gate=0 — the v1.1i
reference: does registry maintenance beat plain gate-tracking?);
**uud_frozen** (same, cant_tau=1e18 — the family bound at this T);
**pair_seed_reg** vs **pair_ctl** (T=1200; the P-REG6 50% bar);
**i5_seed_reg** vs **i5_ctl** (measurement only — winding wall).

Bars: **P-REG5′** uud_seed_reg t_half ≥ 1400 (10× control 140);
report against seed_plain (maintenance value), frozen (bound), and
the CANTUS T=2000 bound 1160. **P-REG6′** pair_seed_reg retention
≥ 50% at T=1200. Medium honesty needs no bar — it is structural
(cant_grow=0; the §4.1 P-Hi1 precedent: all-zero gauges are a
physical no-op, bath byte-identical).

### 4.1 Meter arms M1–M4 (reg_gate=0; logs `runs/registry/m*_t*.log`)

Apparatus verification first: defaults byte-identical to the
committed pre-registry kernel (strong form: -O0 builds of both
sources byte-equal INCLUDING every drift value; at -O3 the drift
diagnostic column alone reassociates at ~1e-16 — compiler artifact,
recorded); meter physics-silent (P-REG2 PASS: reg_tau=30 vs off, all
non-registry lines byte-equal); C≡Go byte-equal with the meter and
with BOTH gate forms firing (class-D drift column and the
pre-existing truss ±0 print exempt as standing); battery ALL GREEN
93 at defaults (`runs/BATTERY_reg0.log`).

**The bath (stationary by t≈300, n≈25,700 live slots):**

| τ | ρ q25/q50/q75 | gross med | frac. with flow |
|---|---|---|---|
| 10 | 0.095 / **0.49** / 0.84 | 6.6e-5 | 0.925 |
| 30 | 0.287 / **0.68** / 0.906 (q90 0.977) | 1.4e-4 | 0.925 |
| 100 | 0.473 / **0.79** / 0.947 (q90 0.987) | 2.7e-4 | 0.925 |

(τ=100 came in ABOVE τ=30 as registered in §3.4 — the balance-
emergence story stands; τ=10 remains the R-G2 selection. The τ=100
ledger also reads the long-lived bodies cleanly: pair ρ 0.952 at
t=600, i5 edges ρ med 0.90 — long windows suit long-lived balanced
bonds and the bath alike, which is exactly why they discriminate
nothing.)

**The objects (alive-phase ledgers, τ=10 unless noted):**

| body | bond gross (E/t.u.) | bond ρ | reading |
|---|---|---|---|
| unison pair 0.47 | 0.014–0.062 sustained | 0.24–0.94 (τ10, chunky); 0.94 (τ30) | standing exchange 200–900× bath median; balance visible on the τ30 ledger |
| UUD triad | 1–4e-3 early, decaying | 0.32–0.49 both τ | **standing DIRECTED current — net ≈ 0.75·gross**: the flux-machine circulation observed per link; first UD bond goes SILENT at t≈40–50 (the D-drain), triad ledger dark by t≈100–140 (= control t_half) |
| i5 ring (unison) | 0.005–0.036 | med 0.6–0.95, q25 dips ~0 | quasi-balanced with one-way episodes |

### 4.2 The vacuum guard as selector (§3.4 menu) — ALL FOUR COMBOS REJECTED

Bare bath T=480, full stack live (k005 + k_cant=1 k_tune=0.2
cant_tau=50 cant_grow=1, reg_tau=10), logs `runs/registry/g_*.log`;
control h4_ctl: cond 2203.97, rad 2016.22:

| combo | cond | rad | nlock | a_max | tune_total |
|---|---|---|---|---|---|
| F-B f0=2e-4 | 862 (−61%) | 351 (−83%) | 2637 | 0.992 | 9615 |
| F-B f0=5e-4 | 874 (−60%) | 366 (−82%) | 2618 | 0.988 | 9096 |
| F-D f0=2e-4 | 831 (−62%) | 313 (−84%) | 2636 | 0.993 | 10553 |
| F-D f0=5e-4 | 862 (−61%) | 347 (−83%) | 2635 | 0.991 | 9910 |

Every bar failed (economy −60/−83% vs the ±10% band; tune 9–10.5k vs
≤1000; nlock at the v1-crystallization level ~2640; a_max 0.99 vs
< 0.5). The f0 doubling moved nothing — the warm tail is not a
static population the threshold can cut; the lock RAISES bath flows
to bond scale as it spreads (the §3.5 flammability mechanism). This
is the third measured ignition of a self-growing coherence rule
(v1 closure/cell-borne, v1.1 closure/link-borne, registry
flow/link-borne) — the family verdict is in §3.5 and the lane's §5.

**Verdicts:** P-REG1 **FAIL as registered** (bath ρ median 0.49–0.68
≫ 0.3 at every warmed τ; the two discoveries and the amended
protocol are §3.4). P-REG2 **PASS exactly**. The reciprocal-flow
scale s = ρ·gross separates alive bonds from the bath MEDIAN by
15–500× (pair/i5) and 15–60× (UUD), but the bath's warm TAIL
(ρ q75 0.84 at gross tail) overlaps the weakest UUD bonds — the
statistical gap is real at the median and leaky at the tail, which
is exactly the situation the amended guard-as-selector (§3.4) was
registered for. **Sharpened arming-race warning for P-REG5:** the
UUD's first bond dies at t≈40–50, before any cant_tau=50 gauge can
arm (~3τ = 150) — the goal arm's registered secondaries (fast-tau,
and the seeded-instrument ceiling reference) carry the burden if
the self-arming primary loses the race.
