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

## 4. Results

(§4 is written only from executed runs; nothing above this line may
be edited after the first run — corrections land as dated addenda.)
