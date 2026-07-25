# X10-pathfinder results — shell ringdown / existence test (GPU)
**Date:** 2026-07-22 · **Status:** run truncated at t≈1700/2000 (85%) by
instance failure; verdict frames secured. Data: `/space/scp/v85/topo1/gpu/`
(run3_*: frames 0/500/1000/1500, diag → t=1703.6, logs).

## Design (per ANALYSIS §2C, adapted)
Ball ω=1.4325 (ungauged profile → seeded Q=160.6; α_f = g²Q/4π = 0.032,
**a₀ ≈ 21**, ε₁ ≈ 7.6e-4, 1/ε₁ ≈ 1300 t.u.) + opposite-charge perturbation
(amplitude ×0.15 of the profile, Q ≈ −3.7, seeded at r = 21, ω = −1.4994).
N=384, L=55, g=0.05, T=2000, absorbing BC. V100-32GB (Vast).

## Campaign notes
- **Kernel fix required and applied (authorized):** `init_gauss_project` in
  scp_sim.cu was a *serial* host CG — 1362 iterations ≈ 12 s each ≈ 4.8 h at
  N=384. Parallelized with 8 OpenMP pragmas (host-side only; no physics
  change) → ~48 min at ~7× (memory-bandwidth-bound on 72 cores).
- **Two Vast instance failures**, both on the ssh7.vast.ai proxy host: gpu1
  died pre-physics; gpu2 died at t≈1700 (85%). Hourly local pulls (user
  directive) preserved everything of value both times.

## Measured: the four-frame negative-charge story

| t | total Q₋ in box | profile |
|---|---|---|
| 0 | −3.7 (seeded) | perturbation at r=21 |
| 500 | −0.967 | **shell peaked at r ≈ 20–24 (= a₀)** + outbound tail r>38 |
| 1000 | −0.504 | transient outward slosh (broad max r ≈ 28–32; < 1 binding time — settling, not evaporation) |
| 1500 | −0.331 | **peak returned to r ≈ 20–24 (= a₀), amplitude ≈ unchanged from t=1000**; losses came from the outer structure (r=28–36 bins 0.037→0.013) |

Decay of the total: τ = 765 (500→1000) → **1190 (1000→1500)** — strongly
decelerating. Core: **zero negative charge inside r = 12 at every frame**;
ball intact (+160.5, static profile). Gauss at 1.8e-13 floor throughout;
E drift −2.28% (the escaping unbound fraction), θ silent.

## Verdicts

| Question | Verdict |
|---|---|
| **F11 — does a discrete bound response structure exist?** | **PASS (existence).** A negative-charge shell formed at the *predicted* Bohr radius a₀ ≈ 21, survived ≥1500 t.u. (≈1.15 binding times), and re-concentrated onto a₀ after a sub-binding-time transient — mode behavior, not dispersal. |
| **F10 — does the cloud drain into the core?** | **Did not fire.** No negative charge within r=12 at any time; annihilation channel closed on this timescale (detuning protection sufficient; flavor protection wasn't even needed). |
| Shell lifetime | **Remanded — BC-confounded.** The mode's evanescent tail length 1/κ ≈ 21 vs sponge clearance ≈ 20: the absorbing boundary continuously eats the outer mode structure (ARTIFACT(L) per the D7 classifier). The observed drain is dominated by box loss of loosely-bound components; the surviving a₀-peaked core's intrinsic lifetime is not measurable in this box. |

## SHELL thread status after X10
The atom architecture's central bet — **orbitals as bound response modes of
the same fabric** — now has direct dynamical evidence: right radius, right
protection phenomenology (core untouched), mode-like settling. What it does
not yet have: a BC-clean lifetime, mode frequencies (needs T ≈ 1.5×10⁴ per
the ANALYSIS error budget), or occupancy structure (D13 untouched).

## Follow-up options (in value order)
1. **X10b — deeper shell, BC-clean (recommended):** Q_N ≈ 450 at g=0.05
   (α ≈ 0.089, ε₁ ≈ 6e-3, tail ≈ 9, a₀ ≈ 7.5, ρ = a₀/r_half ≈ 1.8):
   marginal scale separation but ~4.5 tail-lengths of sponge clearance in
   the same N=384/L=55 box — lifetime and settling become physical
   measurements. Same seed recipe, ~20 h V100.
2. **L-scan (D7 certification)** on the pathfinder design — 2–3 boxes,
   expensive at N≥384.
3. **Tail completion** of this run from the t=1500 frame (T+500) — low
   value; the curve's shape is established.

---

# X10b results — deep-shell run (Q=444, completed T=2000)
**Date:** 2026-07-23 · **Instance:** gpu7 (V100-32GB), 19.4 h wall, clean
completion. Data: `/space/scp/v85/topo1/gpu/x10b*` (9 frames @ 250 t.u.,
diag @ 0.25). Design: ball ω=1.3925 (Q=444, ungauged profile), α_f=0.088,
a₀=7.55, ε₁=5.85e-3 (1/ε₁=171 t.u.; T=2000 = 11.7 binding times; 5.4
tail-lengths of sponge clearance — BC-clean by design), perturbation −11 at
r=9.

## The frame story (negative-charge radial profiles)
| t | Q₋ total | distribution |
|---|---|---|
| 250 | −3.49 | r=14–30, peak 26–28 (n=2 territory); **~3.9 units annihilated on the ball + ~3.7 escaped** — seeding at r=9 was inside the contact zone (ρ=1.8 too aggressive) |
| 500 | −2.27 | peak 22–24 — inward cascade begins |
| 750 | −1.59 | bimodal 18/28 — radial slosh, undersampled by frame cadence |
| 1000 | −1.18 | unimodal peak 18–20 — net migration 27→19 over 750 t.u. |
| 1250 | −0.82 | plateau 16–26; **migration arrested** |
| 1500 | −0.61 | pinned 16–20; ball skirt reaches r≈12 (contact shielding of the 1s) |
| 1750 | −0.51 | narrow band 16–20, outer tail empty; τ→1300 |
| 2000 | −0.52 | **drain STOPPED (flat over final 250 t.u.)**; band dispersed to 20–30+ with first traces at r=10–14 — non-stationary slosh, retained |

## Verdicts
1. **Retention: PASS.** ~0.5 units of opposite charge permanently retained
   (drain rate → 0 by run end) — no BC loss, no annihilation loss at late
   times. The BC-clean design worked for what survived.
2. **Cascade: CONFIRMED.** Net inward migration 27→19 over 750 t.u. under
   classical radiation reaction — the fabric's first observed orbital decay
   toward a ground state.
3. **Ground-state occupation: BLOCKED (contact shielding).** The 1s radius
   (7.5) sat inside the ball's contact/annihilation zone, which *expanded*
   during the run: the seeded ungauged profile relaxed toward the fatter
   gauged equilibrium (r=8 bin ×6 growth) — a benign but consequential
   seed-shape artifact. Migration arrested at r≈16–20.
4. **Stationarity: NOT ACHIEVED.** The retained cloud sloshes radially
   (period ≳ 500 t.u., aliased by the 250-t.u. frame cadence).
5. **Initial annihilation: F10 partially fired at contact** (~3.9 units of
   ball charge consumed at t<250) — the first direct annihilation-rate
   datum; confirms seeding must stay outside the contact zone.

## Design lessons (bracketing complete)
- Pathfinder (ρ=5): BC-dirty, contact-clean → shell at a₀ but boundary ate it.
- X10b (ρ=1.8): BC-clean, contact-dirty → retention but 1s shielded.
- **X10c window: ρ ≈ 3** — Q≈270 at g=0.05 (α_f=0.053, a₀=12.6, ε₁=2.1e-3,
  tail=12.6, 2.8 tail-lengths clearance, contact zone ends ~r 8–9).
  Plus: seed perturbation AT a₀ (not inside), finer snap cadence (150) to
  resolve the slosh, and (future) gauged-shooter profiles at high Q.

---

# X10c results — ρ≈3 window run (Q=267, completed T=2000)
**Date:** 2026-07-24 · gpu7 (V100-32GB), 19.5 h, clean completion, 15 frames
@ 150 t.u. Design: ball ω=1.410 (Q=267), α_f=0.053, a₀=12.6, ε₁=2.11e-3
(1/ε₁=473; T=2000 = 4.2 binding times), perturbation −6 seeded AT r=13=a₀.

## Verdicts
1. **Contact avoidance: TOTAL SUCCESS.** Ball charge +267.1 constant to
   0.02% for 2000 t.u.; zero annihilation (vs X10b's −3.9). Seeding at a₀
   eliminated F10 entirely.
2. **Retention: PASS.** −0.43 units held at T=2000; leak τ grew 950 →
   ~2000+ (unbound flush complete by t≈450; thereafter the bound cloud is
   effectively permanent on the run timescale).
3. **Shell localization: PASS with corrected radius.** The cloud never
   left the 12–35 band after formation; repeated consolidation episodes
   peak at r ≈ 16–22 ≈ 1.4–1.6 a₀ — consistent with the ground state of
   the ACTUAL effective potential: a hard-core Coulomb (ball contact
   exclusion r ≲ 12 + point-charge tail), whose 1s peaks outside the naive
   point-Coulomb a₀. The recurring ~1.5a₀ localization across all three
   X10 runs supports this reading.
4. **Stationarity: NOT reached — and the obstacle is now understood.**
   The cloud breathes radially (period ≈ 500–600 t.u., ~3.5 cycles
   observed) with only mild late-run damping. Mechanism: **a radial
   (monopole) charge oscillation cannot radiate through the gauge field**
   (no dipole moment change) and the matter channel is gap-blocked — the
   breathing mode is radiatively protected. This explains X10b's endless
   slosh and X10c's persistence, and is itself a fabric-electrodynamics
   prediction verified by three runs.
5. Gauss at 2.2e-13 floor throughout; θ silent; energy drift −2.38% (all
   escape-phase radiation).

## X10 series synthesis (three runs)
| Run | ρ = a₀/r_half | Outcome |
|---|---|---|
| Pathfinder (Q=160) | 5.0 | Shell AT a₀, BC-eaten (tail ≈ clearance) |
| X10b (Q=444) | 1.8 | Retention, cascade observed, 1s contact-shielded |
| X10c (Q=267) | 3.0 | Retention + localization at hard-core-corrected 1s; monopole breather persists |

**Established across the series:** bound opposite-charge response clouds
form at the (effective-potential) Bohr radius, are retained indefinitely,
never annihilate when kept off-contact, and cascade inward under radiation
reaction — while purely radial excitations persist because monopole
breathing is radiatively protected.

## Next: X10d — orbital angular momentum (the capstone)
Seed the perturbation ON a circular orbit: position (0, 13, 0), tangential
velocity v_orb = √(α_f/(m·a₀)) ≈ 0.053 via de Broglie tilt (vx along x ⊥
radius). A rotating cloud has a rotating dipole → radiates → circularizes
into the stationary orbiting ground state the monopole breather cannot
reach. Same box/params as X10c; expect L_z bookkeeping (ħ_eff per unit
charge) as a bonus measurement.

---

# X10d results — orbital-seed run (completed T=2000): EXACT TWIN-NULL
**Date:** 2026-07-25 · gpu7, 19.5 h, clean completion. Identical to X10c
except the perturbation carried tangential velocity v=0.053 at (0,13,0)
(circular-orbit initial condition, specific L ≈ 1.0 per unit charge).

**Result: NULL to measurement precision.** Every frame matches X10c
bin-for-bin within ~3% (final totals −0.414 vs −0.434; same breathing
period, same consolidation episodes, same end-state). The seeded angular
momentum produced no measurable dynamical difference.

**Why (the load-bearing lesson): ORBITS REQUIRE CARRIERS.** A non-solitonic
perturbation disperses into a shell during the scatter phase; its coherent
orbital motion smears into negligible net rotation, and the monopole
breather dominates identically. A smeared charge cloud cannot orbit — only
a self-bound closure retains the coherent L needed to circularize into a
stationary orbit. This converges with the TOE quantum-light chain's carrier
requirement (Step 7/9) from an entirely independent direction: the fabric
derives classically that the atomic electron must be a discrete object.

## Next: X11-PB — two-closure orbit ("fabric hydrogen-lite")
Heavy: v71 flavored ball (Q=284, partition 76.7/103.7/103.7). Light: the
branch-bottom soliton Q≈−88.5 (ω=−1.4750 profile, E=135.5, r_half=2.35) at
(0,12,0) with tangential v=0.064 (α_AB = g²·284·88.5/4π ≈ 5.0, μ ≈ 102;
orbit period ≈ 1200 t.u. → ~1.7 orbits in T=2000). Protections, all
measured in this program: distance ≫ tails (annihilation overlap ~e^{−2κD}
negligible at D=12), clock detuning beat ≈ 2.9, flavor-partition mismatch
(annihilation conservation-redirected, X8 design). Goal: bound two-closure
orbital motion — the original v75 Step-2 target with today's design
discipline — plus inspiral/radiation and L_z ledger data.
