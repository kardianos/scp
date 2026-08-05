# RADIANCE — the R campaign: an interior flux fixed point from steep sub-cap emission

Opened 2026-08-04 by the first v91 agent, per `README.md` NEXT STEPS.
Format precedent: `carried/FORGE_EVIDENCE.md` (pre-register, then run).
The law candidate under test is THEORY.md §2.1 (candidate A, implemented
inert in both kernels at k_rad=0). **§1–§3 of this document were written
and committed to disk BEFORE the first k_rad>0 experiment ran.** §4+
records what happened.

The deficiency being repaired (measured, FORGE E1): laws_V2g's
emission-vs-fullness curve is a cap-wall step — outtake EXACTLY 0 below
cap at every fullness, thin det-throttled trickle above — so in > out
everywhere sub-cap, the only fixed point is the wall, and nothing can
hover: no stable mass, no size selection (B1's cause as a diagram).
Candidate A adds one guarded term: on each beat, a graded sub-cap
emission demand k_rad·cap·x^p_rad·comp routed exactly as evaporation.

---

## 1. Apparatus (E1-f protocol, replicated exactly, one knob turned)

One probed cell in a live noisy 3D bath; partner at x=0 far-separated
(in-run control + clean attribution); fixed-cell compensation pins the
picker.

```
./freecell exp=pair bath=1 freeze_geo=1 convtag=1 noise_amp=0.5 \
    T=120 diag_every=200 seed=20260802 \
    pair_x0=<x> pair_x1=0 pair_doff=<7.0 - dstar(x)> \
    k_rad=<k> p_rad=<p> rad_clock=<c>
```

with dstar(x) = 2π·C/(w2/det(x) + w2), det(x) = 1 + q_detune·x
(= 2.16662·det/(1+det) at V2g constants). All other keys at defaults
(laws_V2g verbatim; L=16, dt=0.02, jam 800/0.1). The doff compensation
holds the picker target d0 = 7.0 at every x, so the SAME bath pair
(v90: i=475, j=1565) is probed at every x — and, because the radiance
keys do not enter bath construction or picking, the same cell across
**every arm of the grid**: foam-pocket variance is removed from the
entire campaign, not just the x sweep. Runs land in `runs/radiance/`,
named `r1_<arm>_x<xxxx>.log`.

**Meters.** Intake = ct_cond of `# RESULT convtag` (field captured at
the tagged pair). Outtake = ct_evap + ct_rough. Below cap, real
evaporation cannot fire (needs Em+Ee > cap), so **ct_evap sub-cap IS
the radiance channel, exactly**; above cap the two share the door and
the credit register (R-D4 — recorded, not claim-bearing here). Global
`RESULT conv rad=` counts all cells' radiance; rad_global − rad_tagged
= the vacuum glow. Time-resolved CONVTAG every 4 t.u.; QATOM sampler
prints every 200th fired atom; Em_tag diag column tracks the probed
pair's dense holdings.

**x grid (17 points):** 0, 0.10, 0.21, 0.28, 0.35, 0.4167 (fifth),
0.48, 0.52, 0.56, 0.60, 0.65, 0.70, 0.75, 0.8333 (octave), 0.92,
1.00, 1.05. The four added interior points (0.52–0.70) densify the
predicted crossing zone; the v90 grid is a subset (byte-check arm).

**Arms (coarse, T=120, seed 20260802):**

| arm | k_rad | p_rad | rad_clock | purpose |
|---|---|---|---|---|
| k000 | 0 | – | – | baseline; byte-check vs `v90/runs/forge/e1_f*.log` |
| k001 | 0.01 | 4 | 0 | near-wall marginal band |
| k002 | 0.02 | 4 | 0 | crossing band, low edge |
| k005 | 0.05 | 4 | 0 | predicted selection neighbourhood |
| k010 | 0.10 | 4 | 0 | crossing band |
| k030 | 0.30 | 4 | 0 | strong arm; door-ceiling probe |
| c1k005 | 0.05 | 4 | 1 | ablation: ride the slowing beat |
| p2k005 | 0.05 | 2 | 0 | steepness down (shallow) |
| p6k005 | 0.05 | 6 | 0 | steepness up |

153 runs. Then a refinement stage at the selected point: T=480,
x ∈ x*±0.10 in 0.02–0.03 steps, seeds {20260802, 111, 314159}
(measured x*, slope, seed-robustness) before R2.

**Instrument lessons carried (README rule 5):** picker pinned by doff
(above); sub-grain demands accrue credit silently — a zero can mean
"below one atom over this T" (threshold table in §2); ledger claims
only (no single-seed angular claims); sector meter not used in R1.

### 1b. Refinement stage R1b (pre-registered 2026-08-04 after the k-ladder, BEFORE R1b ran)

Arms k000/k005/k010; T=480 (4× grain resolution), diag_every=100
(CONVTAG window = 2 t.u.); x ∈ {0.52 … 0.80 step 0.04}; seeds
{20260802, 111, 314159}; 72 runs (`run_r1b.sh`). Deliverables: the
flux curve F(x) = Δ(ct_cond − ct_evap)/Δt pooled over windows binned
by instantaneous x; its zero x̂* and slope; Em(t) relaxation; seed
robustness. Pre-registered: x̂*(k005) ∈ [0.60, 0.72]; x̂*(k010) ∈
[0.50, 0.60]; dF/dx < 0 (restoring) at both; k000 has NO zero in the
band (F > 0 everywhere sub-cap: V2g fills monotonically).

*Amendment (dated 2026-08-04, before R1b ran):* trajectory analysis
of the coarse runs exposed a partner blur — the x1=0 partner
acquires ~0.3–0.4 Em from the churn, and Em_tag (a SUM over the
tagged pair) then overstates the probe's holdings (measured: the
first radiance atom at x0=0.8333 fires ~3× later than the demand
integral predicts unless the partner carries ~0.35). R1b therefore
uses the SYMMETRIC pair (pair_x1 = x0, the kernel default): per-cell
flux = ledger/2, per-cell Em = Em_tag/2 exactly, two independent
pocket samples per run. The picker target stays 7.0 (same pinned
pair). Per-cell predictions above are unchanged.

## 2. Closed-form expectations (derived from the implemented term, before any run)

Beat rate of a cell at fullness x: f_beat = (w2−w1)/(2π·det)
= 0.19894/det per t.u. Radiance demand per beat: k_rad·cap·x^p·comp.

* **rad_clock=0 (comp = det): the det cancels exactly.** Demand rate
  per unit time is the pure power law

      D(x) = k_rad · cap · x^p · (w2−w1)/2π = 0.49736 · k_rad · x^p

  — this is the design content of R-D3 (a fuller body radiates more
  per unit TIME, on a clock that does not slow).
* **rad_clock=1 (comp = 1): D(x)/det** — suppressed by 1+1.2x
  (×0.74 at x=0.35 … ×0.45 at x=1.0). The ablation does not abolish
  sub-cap emission; it re-introduces the slowing-clock throttle. The
  full V2g pathology (outtake ≡ 0 sub-cap) is the k000 row only.
* **Grain:** each fired atom is ε(x) = A0·w2e/2π = 0.53084/det
  (0.474 at x=0.10 … 0.235 at x=1.05). Credit carries between beats
  (two-atom register), so fired totals track demand within ±1 atom
  over any window. **One-atom threshold at T=120:** first fire needs
  0.49736·k·x^p·120 ≥ ε(x); for p=4 that is x ≥ 0.54 (k=0.01),
  0.46 (0.02), 0.37 (0.05), 0.31 (0.10), 0.24 (0.30). Zeros below
  threshold are silence, not absence of the law.
* **Door bandwidth ceiling:** the register caps at 2ε per beat ⇒
  outtake rate ≤ 2ε·f_beat = 0.21124/det² per t.u. (0.0436 at x=1).
  The ceiling binds when k_rad·x^p·det² > 0.42479 — for p=4 at x=1,
  k_rad > 0.088. The k030 arm near cap is predicted to BEND below the
  power law (door saturation); k005 and below never touch it sub-cap.
* **Intake reference (v90 e1_f, per 120 t.u., committed logs):**

  | x | 0 | 0.35 | 0.4167 | 0.65 | 0.8333 | 0.92 | 1.00 | 1.05 |
  |---|---|---|---|---|---|---|---|---|
  | in | 0.000 | 0.515 | 0.601 | 0.170 | 0.159 | 0.460 | 0.444 | 0.447 |

  Non-monotone: fifth peak, a 3.5× dip at 0.65–0.83, recovery near
  cap. (The FORGE prose quotes "0.44–0.60" and does not mention the
  dip its own logs contain; the dense grid re-measures it — the
  k000 arm is the honest in(x).) NOTE: at the fifth (0.4167) and
  octave (0.8333) the probed cell's pitch is comb-consonant with the
  bath and dense BOND exchange opens (Em_tag swings ±0.3 in the v90
  logs); in the gate-closed band x ∈ [0.48, 0.80] (ratio strictly
  between 3:2 and 2:1, p·q ≤ 6 admits nothing) the Em budget is
  ledger-only — balance claims are cleanest there.

## 3. Pre-registered predictions (BEFORE first k_rad>0 run)

* **P-R1 (baseline byte-check):** the k000 arm reproduces the eight
  v90 `e1_f*` logs byte-identically (same picked pair i=475/j=1565,
  same CONVTAG totals) at the shared x points.
* **P-R2 (sub-cap emission exists):** every k_rad>0 arm fires
  ct_evap > 0 for all grid x above its §2 one-atom threshold — the
  cap-wall step is abolished by the candidate term.
* **P-R3 (the power law, quantitative):** in the gate-closed band,
  measured outtake totals match 0.49736·k·x^p·120 within ±1 atom
  (rad_clock=0). Deviation beyond grain noise falsifies the det
  cancellation (R-D3) or the implementation understanding.
* **P-R4 (interior crossing exists — THE claim):** out(x) crosses
  in(x) at some x* < 1 for k_rad ∈ [0.02, 0.30], p=4. Two branches,
  by what in(x) turns out to be:
  (a) the v90 dip is real → the intake cliff after the fifth pins
  x* into [0.43, 0.65] across the whole k_rad decade (a comb-shaped
  intake quantizing the fixed point);
  (b) in(x) is smooth ≈ 0.42–0.60 → x* ≈ (0.00840/k_rad)^{1/4}:
  0.81 (k=0.02), 0.64 (0.05), 0.54 (0.10), 0.41 (0.30).
  Either way at least one arm yields an interior crossing. k001 is
  predicted marginal (x* ≈ 0.95, near-wall).
* **P-R5 (stable slope):** at every measured crossing,
  d(out−in)/dx > 0, dominated by out′ = p·out/x ≈ 4·in(x*)/x*.
* **P-R6 (ablation ratio):** c1k005 outtake = k005 outtake × 1/det(x)
  pointwise within grain noise — the beat-stretch measured directly;
  its crossing shifts wallward accordingly.
* **P-R7 (vacuum stays dark; steepness is the guardian):** the bath
  never crosses its own one-atom threshold at p=4: rad_global −
  rad_tagged = 0 exactly for k ≤ 0.05; ≤ a few atoms at k030.
  Ordering across steepness arms: glow(p2) ≥ glow(p4) ≥ glow(p6),
  with p6 exactly dark. (Radiance must not cook the vacuum — this,
  not analogy, is why emission must be STEEP.)
* **P-R8 (radiance is comb-blind):** out(x) is smooth through the
  fifth and octave (no resonance structure beyond grain steps) — it
  fires on the cell's own beat, gated by nothing. Contrast: in(x)
  keeps its fifth structure. (Consonance gates transport; radiance
  is thermal.)
* **P-R9 (conservation + purity):** |drift_rel| ≤ 1e-13 every run;
  no kernel or table change in R1, battery untouched and green.
* **P-R10 (the balance integrates):** in the gate-closed band, the
  probed cell's final Em − seeded Em is positive for x < x* and
  negative for x > x* (the static precursor of R2's perturbation
  return). Not claimed at the fifth/octave (bond exchange dominates
  there).

### 3.1 Addendum (2026-08-04, post-probe, pre-sweep — dated, nothing above rewritten)

The single k000 probe at x=0.65 (run before any k_rad>0 experiment, as
the timing/byte-check) exposed a v90 apparatus defect: the v90 f-arm's
doff at x=0.65 and x=0.8333 was mis-tabulated by +0.009 (target 7.009,
not 7.000), flipping the picker to a DIFFERENT pair (564/1650) at
exactly the two anomalous dip rows; the six rows that picked 475/1565
all read 0.44–0.60. The v90 "intake cliff" is therefore suspected
cross-pair contamination, and the FORGE prose range "0.44–0.60" was
unknowingly quoting the uncontaminated subset. Consequences, recorded
before the sweep runs:

* **P-R1 (amended):** byte-identity is claimed against the six
  pair-475/1565 logs (f000/f035/f042/f092/f100/f105) only; at 0.65
  and 0.8333 the v91 pinned-pair k000 runs REPLACE the baseline.
  (v91 uses doff = 7.0 − d*(x) exact at every x; at a constant
  target the picker cost is x-independent, so one pair everywhere.)
* **P-R4 weight shifts to branch (b)** (smooth in(x) ≈ 0.44–0.60;
  x* ≈ (0.00840/k_rad)^{1/4}). Branch (a) revives only if the
  pinned pair still shows the dip.
* **P-R11 (new):** pinned-pair in(0.65) and in(0.8333) land in
  [0.35, 0.70] per 120 t.u. — the dip was apparatus. If they land
  at 0.16–0.17 instead, the dip is physics and branch (a) stands.

**Pre-registered selection rule for (k_rad, p_rad):** the point with
(i) an interior crossing inside the gate-closed band [0.48, 0.80]
(radiance-owned, not comb-owned), (ii) ≥ 3 fired atoms on each side
of x* over the window (grain-resolved slope), (iii) the smallest
k_rad satisfying (i)+(ii) (minimal new-law footprint), (iv) vacuum
dark at that k (P-R7). Algebraic expectation: (k_rad, p_rad) ≈
(0.05, 4), x* ≈ 0.55–0.65. The selection is made by the measured
curves, not by this expectation.

---

## 3.2 R2 protocol pre-registration (written after R1 coarse, BEFORE any k>0 bath-object run)

Vacuum arms (already run as cavity probes): a k>0 object in true
vacuum is a CLOSED CAVITY — its radiated field can only recondense
on itself; predicted and observed (ring6: Em ret 0.82 plateau over
20000 t.u., rad 29.7 vs recond 25.0, all 6 links retained vs the
frozen k0 control's 6→4 link death; UUD: ret 0.91, D-bloat and dead
fifth UNCHANGED — radiance does not do the coherent channel's job).
The one vacuum leak is the space backsplash: bs/rad = s_pull/(1+s_pull)
= 13% of fired radiance → cavity lifetime ≈ Em/(0.13·rad_rate) ≈
2.5e4 t.u. at k=0.05.

Bath arms (run_r2.sh, T=1500–5000): the seeded voice loads (0.28,
0.325) sit BELOW x* ≈ 0.6; under radiance-balance voices should FILL
toward x*, which moves every bond rung (unison d* = πC/w grows ~29%
from x=0.28 → 0.60). Pre-registered outcomes, either informative:
(a) re-seat — the object swells onto the x* rungs and persists
(size selection acting on a real object); (b) rung mismatch breaks
the bonds faster than k0 dissolution (radiance–composite
incompatibility at this ambient). Because the churn decays, x*(t)
falls slowly; a balanced object TRACKS it. R2 bars therefore:
(i) topology retention (live links/gates at T, vs the k0 control);
(ii) Em tracking the moving fixed point — neither wall-pinning nor
drain-to-zero; (iii) windowed ledger balance (intake ≈ outtake at
the object, the flux-machine criterion); (iv) conservation ≤ 1e-13.

## 6. CAMPAIGN VERDICT (all runs complete; ~470 runs total, drift ≤ 1.35e-15 throughout)

**R2 final table (t_half to ret=0.5; terminal ret):**

| object (bath, noisy ambient) | k005 t_half | k0 t_half | k005 terminal | k0 terminal |
|---|---|---|---|---|
| ring6 @0.325 | 80 | 510 | **0.298 @5000** | 0.242 @5000 |
| composite nv=6 | 90 | ~490 | 0.153 @5000 | 0.076 @1500 |
| UUD | 140 | 260 | 0.170 @2000 | 0.192 @2000 |
| blob (loose) | >2000 | >2000 (grew 1.16 peak) | 0.607 @2000 | 0.966 @2000 |
| ring6 @0.62 (x*-built) | 140 | 410 | 0.130 @~5000 | 0.127 @~5000 |

The campaign, in five measured sentences:

1. **The door balance exists and is exactly the designed law.**
   Interior fixed point x̂* = 0.62 ± 0.02 (k=0.05, p=4), quarter-power
   in k AND p across the swept decade, beat-stretch ablation verified
   to 3%, V2g control at literal zero in 24 runs, conservation
   ≤ 1.35e-15 with the term firing globally, C≡Go byte-faithful
   under load.
2. **The substrate acquired a thermal medium.** The churn nucleates
   proto-matter that GLOWS (supply-throttled, ~constant in k;
   steepness p ≥ 4 keeps the low-x bath dark); a measured
   radiate→recondense economy runs everywhere the law regime runs.
3. **Radiance taxes hoards.** Everything built under V2g's sealed
   door — growing blobs, D-pile composites, parked rings at any x —
   is pumped down to the glow-economy floor 3–6× faster than the
   V2g leach; but the ruins then TRACK the dim ambient's fixed
   point and HOLD a floor the k0 objects sink through (ring k005
   0.298 > k0 0.242 at t=5000). The balance governs everything; it
   preserves no v90 structure.
4. **Structures cannot yet HOLD the balance — measured from both
   sides.** In bath, the Em-pile (Jensen on x⁴) + flight-loaded det
   + roughness amplify a bonded object's outtake ~2.3× past the
   single-cell law and pump it out; in vacuum AT the fixed point,
   the thermal atom-kick jitter (~10% of holdings per fire) breaks
   magnitude-gated bonds (all channels die while 91% of energy
   persists). **The flux-machine interior — the S2 coherent channel
   that can hold an internal current and flatten the pile — is
   therefore MEASURED as load-bearing for stable mass.** The two
   v91 modifications are a coupled pair; radiance alone builds
   balance points and cavities, not particles.
5. **Forging follows density under the new law and the forge
   glows** (byte-exact E3 replication at k0; −44/−67% void
   starvation, +12% glow-lift, region recirculation at 2×) — but
   loose forged matter is a through-flux zone, and the star-loop's
   finite stable product awaits the coupled construction.

**Recommendation (agent's, for the user to accept or reject):** do
NOT create laws_V3r now. Acceptance-surface items 2–5 (object
perturbation-return, 10× lifetimes, size selection, forging
saturation) are measured to be blocked on the coherent channel, and
the ratchet forbids adopting a table that cannot carry its own
acceptance surface. Keep k_rad wired inert (defaults 0 — battery
untouched and pure), keep this document as the standing record, and
open stage 2 (S2 coherent channel, `carried/S2_CHANNEL.md`) with
radiance as its test harness: S2's acceptance surface should now
INCLUDE "a coherent object sits at x̂* with intake = outtake and
survives ≥ 10× the v90 ceilings" — the two halves adopt TOGETHER as
laws_V3 when both pass. The §5 reckoning table stands ready if the
user instead directs adoption of radiance alone.

## 5. R5 — the ratchet reckoning (data complete; DECISIONS ARE THE USER'S)

Pure battery after the harness `-extra` flag: **ALL GREEN (93)** —
the flag off is byte-inert, certified. Full battery with
`-extra "k_rad=0.05 p_rad=4 rad_clock=0"` (`runs/BATTERY_r5x.log`,
logs `runs/r5x/`): **82/93 green with the candidate live** — every
drift/conservation bar passes; C≡Go byte-equality holds on every
masked arm with radiance firing globally. The 11 movers, classified:

| class | bars | what moved | proposed disposition (user decides) |
|---|---|---|---|
| A. flag-day markers | conserve c/go law-purity (2) | header now prints k_rad=0.05 | at adoption: purity line re-pins to the laws_V3r defaults (mechanical) |
| B. sealed-door physics, law working as intended | ring6 edge_dev 0.153 vs ≤0.05; ring6 min gg 0.000 vs ≥0.9 | the vacuum ring is now a WARM CAVITY (thermalized: bonds breathe, gates flicker, links + energy retained — §3.2 cavity result) | re-gauge to cavity observables (link retention, Em plateau) |
| | blob retention 0.453 vs ≥0.78 | hoard-tax: the bath blob pump-radiates (same physics as §4.7) | re-gauge or re-protocol at x* seed |
| | pauli0 at-cap admission 0.31/0.018 vs =0 (2) | the cap is SELF-RELIEVING under Stefan balance — full matter sheds then re-admits | re-protocol: instantaneous-on-beat refusal vs time-integrated; the identity-blindness control needs a live-door design |
| | xsec headroom net 6.70 vs [6.9,7.65]; pure-cond evap 0.38 vs =0 (2) | the absorber is WARM — re-emits ~5% of what it absorbs (Kirchhoff-adjacent; the reality-correct direction) | re-gauge band + replace "pure cond" with an emissivity bound |
| C. threshold graze | p1 v_COE 5.16e-5 vs ≤5e-5 | 3% over a bound pinned at the V2g floor; no-transport physics intact | re-pin to the k005-measured floor |
| D. port bookkeeping | abx blob2+p1 byte-equal (1) | drift-column-ONLY divergence (libm pow last-ulp via x^p_rad; every physics digit byte-equal — verified by diff) | mask drift on blob2 like the other arms, OR integer fast-path for integer p_rad in BOTH kernels (kernel change — requires explicit authorization), OR leave red as a standing marker |

**No bar moved for a conservation, determinism, optics (ds/ds1 all
green), FCQ, p2lc, or cross-kernel-physics reason.** The moves are
the new law being a law: doors that seal became doors that breathe.

## 4. Results

### 4.4 R1b refinement — the gated balance curve, and the two-timescale truth (72/72 runs)

**The gated result (early window t ∈ [0,120] = the churn plateau,
known-x0 attribution, 3 seeds pooled, per cell per t.u.):**

| x0 | k000 in/out | k005 in/out | k010 in/out |
|---|---|---|---|
| 0.52 | 0.0022 / **0.0000** | 0.0028 / 0.0023 (+) | 0.0035 / 0.0031 (+) |
| 0.56 | 0.0023 / **0.0000** | 0.0025 / 0.0021 (+) | 0.0027 / 0.0044 (−) |
| 0.60 | 0.0015 / **0.0000** | 0.0031 / 0.0029 (+) | 0.0023 / 0.0044 (−) |
| 0.64 | 0.0016 / **0.0000** | 0.0024 / 0.0032 (−) | 0.0028 / 0.0060 (−) |
| 0.68 | 0.0016 / **0.0000** | 0.0034 / 0.0041 (−) | 0.0036 / 0.0063 (−) |
| 0.72 | 0.0018 / **0.0000** | 0.0031 / 0.0045 (−) | 0.0022 / 0.0067 (−) |
| 0.76 | 0.0020 / **0.0000** | 0.0029 / 0.0048 (−) | 0.0039 / 0.0079 (−) |
| 0.80 | 0.0017 / **0.0000** | 0.0030 / 0.0059 (−) | 0.0041 / 0.0081 (−) |

* **k000: outtake exactly zero in all 24 runs** — the V2g cap-wall
  step at refinement strength (pre-registered no-zero PASS).
* **x̂*(k005) = 0.62 ± 0.02 (pooled)** — inside the pre-registered
  [0.60, 0.72]. **x̂*(k010) = 0.54 ± 0.02** — inside [0.50, 0.60].
  Ratio 1.148 vs quarter-power prediction (k010/k005)^{1/4} = 1.189:
  3.5%, within bin resolution. Restoring slope at both crossings
  (dF/dx ≈ −0.025 per unit x per t.u. at k005).
* Per-seed honesty: a single seed's window is a ±1-atom instrument;
  signs at the bracket wobble per seed, and the local intake varies
  1.5–3× across pockets — each pocket has its own x* (spread ±0.07,
  x* = (in_local/0.0249)^{1/4} — foam variance propagated through
  the quarter-power law). The pooled curve is the instrument; every
  seed shows the +→− structure within the band.

**The two-timescale truth (the R1b headline beyond the numbers):**
the door balance is FAST (τ ≈ 30–100 t.u.; the x0=0.80 trajectory
relaxes 2.00 → 1.61 per cell by t=58, radiance-dominated), and then
the fixed point MOVES: the seed churn spends itself (~t 150–200),
local intake halves, and x*(t) = (in(t)/0.497k)^{1/4} slides down —
the cell TRACKS the dying ambient (measured: hover at ~1.4–1.6
through t≈150, then sag to 0.38·cap by t=478 while the ledger stays
active). On top of this, transient consonant locks with nucleated
bath matter leak dense on τ ~ 300–500 (pocket-dependent; the SAME
channel as composite dissolution t_half ≈ 490, and present at k=0:
the k000 endpoint scatter 0.26–1.73 from identical seeds-by-x0).
Consequences, recorded: (i) a single free cell is not an object —
asymptotically it equilibrates with the medium, door-balanced or
not; (ii) the stationary balance claim lives in the early window BY
CONSTRUCTION of this apparatus (seed-only noise → transient
ambient); (iii) the fixed point tracking a falling ambient is the
Stefan balance behaving like physics, not a failure mode. R2 —
structured objects whose internal bonds are consonant (no internal
leak) and whose x* DETUNES them from the bath's unison channels —
is where stable mass lives or dies, exactly as the campaign always
said. Pre-registered for R2 (§3.2, extended): the k005 composite's
dissolution curve should be TWO-PHASE — early unison-leak while its
voices fill from 0.28 toward x*, then suppressed leak once detuned
past the bath's consonance windows — with t_half ≫ the k0 ceiling
490.

**Methodological lessons (recorded for the battery era):** binning
flux by the symmetric pair's MEAN x is unsound once the pair
diverges (fires from the high cell land in low-mean bins — the §4.4
first-pass flux table showed spurious negatives at x ≈ 0.1);
endpoint distributions at T ≫ 120 measure the ambient economy, not
the door. The early-window known-x0 protocol is the standard.

### 4.5 R4 at E3 scope — forging follows density AND the forge now glows (8 runs + saturation pair)

The v90 v_law protocol verbatim (L=64 slit beam, void grad_r0=6
grad_frac=0.15, region ledger tag_r=6 convtag=1; instrument lesson:
convtag=1 is required for the region print — v90's headers do not
echo it). **k0 arms byte-replicate v90's E3 exactly** (region cond
4.099224 / 6.464811 / 1.207587 / 2.974586 — all four digits-equal
to `v90/runs/forge/v_law_*.log`).

| region cond (r<6) | s20260802 | s111 |
|---|---|---|
| ctl k0 | 6.465 | 2.975 |
| lump k0 | 4.099 (−37%) | 1.208 (−59%) |
| ctl k005 | 7.268 (+12% vs k0) | 2.716 |
| lump k005 | 4.099 (−44% vs its ctl) | 0.906 (−67%) |

Findings: (i) density-directed forging is PRESERVED under the
radiance law (voids starve, same sign, slightly deeper); (ii) the
uniform-bath control condenses +12% more at k005 — the glow-lift of
§4.1 acting at forging scale; (iii) the forged trail RADIATES
(global rad 5.28 lump vs 2.45 ctl at k005; 0 at k0) — the sink
corpse of E1/B7 is replaced by glowing matter. **Saturation pair (T=400, beam maintained):** neither naive
outcome. k0 does NOT grow unboundedly — the upstream trail
shadow-starves the void (cond stalls at 6.33 by t=76; region Em
freezes at ~7: a hoard by starvation — E3's own lensing preempts
the sink at this geometry). k005 reaches the same early fill, then
the door opens: cond runs to 14.5 (2× k0 — radiate→recondense
recirculation measured in the region) with evap tracking it (14.0),
and region Em declines to 2.5 — a RADIATING THROUGH-FLUX ZONE that
slowly re-exports the gulp, not a stable finite product. Verdict:
forging-follows-density and the glowing forge are demonstrated;
**"saturates into a finite stable object" is NOT demonstrated for
loose matter — closure hinges on the structured object at x*
(§4.7), exactly the R2 question.**

### 4.6 R3 frame — the zoo at the measured fixed point

At x̂* = 0.62: w(x̂*) = 1.6628, det = 1.744. The finite size table
the comb admits for radiance-balanced voices (d* = 2πmC/((p+q)w)):
unison m=1 rung = **πC/w(x̂*) = 1.889** (the headline prediction for
R2 survivors' bond lengths); fifth/4:1 line 0.756; octave 1.260;
full table via `zoo_r3.sh 0.62`. v90 objects were built at x=0.28
(unison rung 1.448) — a radiance-balanced ring re-seats +30% wider
IF it survives the transit (the §3.2 (a)/(b) fork). Measured bond
lengths from R2 land here when its runs complete.

### 4.7 R2 in bath, mid-course — radiance taxes hoards; the fixed-point-built object (pre-registered)

**The v90-era objects (built at x = 0.28–0.325) are DESTABILIZED at
the selected point.** Measured: comp6 t_half ≈ 85 vs the k0 ceiling
490 (6× faster); its ledger shows the channel — ct_evap = 19.0 by
t=180 (Em lost 16.4): the anti-Stefan D-pile sits near cap and the
glow-fed ambient tips it over continuously — the composite BOILS
OFF through the door, not out through bonds. The unison ring6 (no
D-pile): ret 0.21 at t=200 vs k0's 0.82; its ledger: intake 2.97
(3× the k0 ring — glow-pumped) but evap 5.54 at rates matching
x ≈ 0.55 radiance — Em piles unevenly on the ring and the piled
cells pump mass out (pile-then-radiate). Honest generalization:
**these objects were "stable" only because V2g's door was sealed
below cap — they are static hoards, and a working Stefan law taxes
hoards.** The flux-machine thesis (REALITY B1) said stable matter
must be MAINTAINED by through-flux; v90 objects are not.

**The constructive experiment (launched; pre-registered):** a ring
BUILT AT THE FIXED POINT — ring_x = 0.62 on its own rung (d* =
1.889, the R3 size) — is the through-flux object: door-balanced
(in ≈ out ≈ 0.003/cell) AND detuned outside every comb window of
the empty bath (pitch ratio 1.744 — boundary leach suppressed, the
§3.2 protection). Pre-registered: (i) ring6@0.62/k005 bath holds
ret ≥ 0.6 at t=5000 with windowed door balance (flux-machine
criterion MET — the B1 answer); (ii) its bond length stays on the
1.889 rung (R3 size selection); (iii) the k0 control at 0.62 is NOT
balanced (V2g: in > 0 = out — fills toward the wall or leaches);
(iv) vacuum@0.62/k005 = cavity plateau. Guard: if (i) fails by the
same pile-radiate pump-out, single-ring topology cannot hold the
balance in a glowing ambient at any x — stable mass then requires
the coherent channel's internal currents (the S2 coupling), and
that becomes the measured statement.

### 4.8 The vacuum endgame and R3's status

* **ring6@0.325 vacuum k005 (T=20000): the one true survivor** —
  ret 0.82, ALL SIX BONDS ALIVE (the frozen k0 control lost 2), and
  the bonds re-seated at d ≈ 1.66–1.68 vs seeded rung 1.506: +11%,
  consistent with the ladder tracking the flight-loaded pitch
  (holdings x 0.26 + flight ≈ 0.19 ⇒ rung 1.67). A flight-resolved
  meter would gate this; recorded as consistent-with.
* **ring6@0.62 vacuum k005 (T=20000): energy survives, the OBJECT
  does not** — ret 0.910 (best retention of the suite: rad 323
  cycled, cond 280 recaptured, backs 42) but ALL SIX CHANNELS DEAD:
  at x* the radiance cycling is 16× the 0.325 level and each atom
  kicks holdings ~10% — the pitch jitter breaks the bond gates'
  phase-lock (p_gate=8). Six disconnected glowing cells remain.
  **Thermal loudness at the fixed point exceeds what a magnitude-
  gated bond can hold — the same S2 conclusion from the other side.**

  **[CORRECTED 2026-08-05 — CANTUS §3.2. The mechanism above is
  RE-ATTRIBUTED.]** The committed log (`runs/r2/ring6x62_vac_k005.log`)
  shows `slots=0, live=0` from t=0 to t=20000: this arm was DEAD AT
  BIRTH — the six cells never connected, so "jitter breaks the gates"
  was never measured here. The true wall is GEOMETRIC: seeding at
  x=0.62 draws s_pull·x·cap from Es, the radii deflate, and the
  unison rung d*(0.62) = 1.889 exceeds the candidate cutoff
  cfac·2r = 1.790 — unison matter cannot even CONNECT at the fixed
  point; only p+q ≥ 5 intervals fit inside contact ("the comb admits
  only CHORDS at the balance"). The S2 conclusion survives in
  sharpened form: the balance requires CHORD structure, not (only)
  gate-holding against jitter.
* **R3 verdict:** no bonded object holds the selected fixed point in
  this suite (bath: pump-out; vacuum@x*: never connected — see the
  2026-08-05 correction above; the "gate-break" reading is retired),
  so the size-selection measurement HAS NO MEASURAND yet. The zoo table (§4.6,
  unison rung 1.889 at x̂*=0.62) stands as the standing PREDICTION,
  testable the day an object can sit at the balance — i.e., after
  the coherent channel. The one measured size datum: the surviving
  0.325 cavity's +11% re-seat, ladder-consistent.

### 4.3 Steepness arms and the coarse-level SELECTION (sweep complete, 153/153, drift ≤ 1.35e-15 throughout)

* **p2 (shallow):** out ratio p2/p4 tracks 1/x² at the multi-atom
  points (2.08 vs 2.04 at x=0.70; 1.55 vs 1.44 at 0.8333); the fixed
  point is SOFT — endpoints rattle 0.14–1.25 (p4: 1.26–2.11) — and
  the vacuum leaks: rad_glob 228 vs 180 (+27%), the low-x bath tail
  firing. Shallow radiance is disqualified on both counts.
* **p6 (steeper):** endpoints cluster 1.67–1.95 on the generalized
  prediction x* = (in/(0.497·k))^{1/p} = 0.743, Em* = 1.86 — the
  power law CONFIRMED IN p as well as k. Glow 174 ≈ p4's 180 (the
  supply throttle saturates for p ≥ 4). Depth 6.6 atoms vs p4's 5.5
  — competitive, inside the ±1-atom blur of this stage.
* **SELECTION (coarse):** (k_rad, p_rad, rad_clock) = **(0.05, 4, 0)**
  — the pre-registered rule (crossing mid gate-closed band at the
  smallest adequate k, vacuum at the glow floor, fixed point ≥ 5
  atoms deep), with p=4 over p=6 by parsimony and the 3D mode-count
  prior (R-D2); p6 recorded as viable. R1b (running) provides the
  gated x̂*, slope, and seed-robustness at the selected point.

### 4.0 The k000 baseline arm (first arm complete; logs `runs/radiance/r1_k000_*`)

**Pinning works:** pair 475/1565 picked at all 17 x (the exact-7.0
target). **P-R1:** byte-identical to v90 at x ∈ {0, 0.35, 0.92, 1.00,
1.05} modulo exactly two known v91 format lines (radiance header;
`rad=` appended to RESULT conv). At x=0.4167 v90 passed the literal
0.41667 vs our 0.4167 (Δx = 3e-5, same pair; the foam decorrelates —
intake 0.601 there vs 0.301 here, a 2-atom vs 1-atom count: the
measured CHAOS FLOOR of a T=120 intake number). At 0.65/0.8333 v90
probed the wrong pair (§3.1) — v91 replaces the baseline.

**P-R11 CONFIRMED — the v90 dip was apparatus:** pinned-pair
in(0.65) = 0.476, in(0.8333) = 0.464, squarely in the fixed-pair
range. Branch (a) of P-R4 is dead.

**The cap-wall step replicates on the pinned pair:** outtake exactly
0.000000 below cap at every x; at x=1.05 evap 0.236 (thin trickle).
FORGE E1's headline stands.

**New instrument finding — in(x) at T=120 is atom-count-limited:**
per-run intake totals are 0–2 fired atoms everywhere (zeros at
x=0.52, 0.56; 0.19–0.30 at 0.21–0.28, 0.48; ~0.46 flat at
x ≥ 0.60). The apparent v90 "smooth 0.44–0.60" and the fifth peak
were small-integer count differences, not resolved structure. In-rate
resolution scales as ±ε/T: the CROSSING is locatable at T=120 only to
±1 atom; the refinement stage (T ≥ 480) carries the x* measurement,
for in(x) and out(x) alike.

**Em-budget caveat (P-R10 scope narrowed):** even in the gate-closed
band the tagged pair's Em drifts by ±0.2 with zero ledger activity
(residual bond/flight exchange with the bath, e.g. Em_fin − Em0 =
−0.18 at x=0.48 with in=0.19, out=0). The diag Em_tag alone cannot
carry a ±0.2-level balance claim; the LEDGER carries R1. P-R10 is
testable only where the drift exceeds that floor.

### 4.1 k001 + k002 arms — the bath glows, and the fixed point shows itself

**P-R7 REFUTED as written, instructively.** rad_glob ≈ 165–187 at
BOTH k=0.01 and k=0.02 — the "vacuum" is not dark, and the glow is
k-INDEPENDENT (supply-throttled). Mechanism, pinned by the numbers:
bath Em is seeded 0 (E0 = 2700·Es + 672 field churn exactly), the
churn CONDENSES CLUMPILY (global cond ≈ 533 at k000 — the
cloud-chamber vacuum of REALITY B7), and the nucleated hotspot tail
radiates. cond_glob rises ~533 → ~620 at k>0: the radiated field
partially RECONDENSES — a measured radiate→recondense loop. The
sub-cap zeros/fires at low x (identical 0.398386 at k001 AND k002,
x=0.28) are cap-door spillover of ambient field spikes riding the
slower-decaying glow-fed ambient, not radiance of the probed cell.
Consequences: (i) sub-cap ct_evap is radiance-dominated only where
the demand ≥ grain (high x); (ii) the E1 ambient at k>0 is
glow-fed and MORE stationary than V2g's decaying churn; (iii) the
balance curve's in(x) at k>0 is higher than k000's (measured:
0.60–0.87 vs 0.44–0.46 at x ≥ 0.92) — the crossing must use
same-run in/out, as planned.

**P-R9 conservation live:** |drift_rel| ≤ 1.2e-15 on all completed
k>0 runs (51 so far) with the term firing globally ~170 units/run.

**The convergence signal (k002, T=120, single runs):** the three
rows seeded ABOVE the predicted fixed point (x0 = 0.92, 1.00, 1.05;
Em0 = 2.30, 2.50, 2.625) all drained to Em_fin = 1.89–1.96, and the
x0=0.75 row (Em0 = 1.875) STAYED (1.894) — against the closed-form
prediction Em* = 1.97 (x* = 0.789) for k=0.02 at the k000 intake
rate. Rows below x* are grain/bond-noise dominated (±1 atom), but
the from-above convergence is monotone-erased across three
independent seeded heights. Pre-registered follow-up: k005/k010/k030
must converge to Em* ≈ 1.61/1.34/1.02 respectively if the crossing
law x* = (in·2.011/k)^{1/4} holds by convergence.

### 4.2a Ablation arm (c1k005 vs k005) — P-R6 CONFIRMED at the clean points

Pointwise out ratio (rad_clock=1 / rad_clock=0) vs the predicted
1/det(x): at the two radiance-dominated, multi-atom points the match
is 3%: **0.494 vs 0.500 (x=0.8333), 0.461 vs 0.475 (x=0.92)** — the
beat-stretch measured directly. Above cap the ratio rises to
0.69–0.70 (real evaporation, clock-blind, mixes into the shared
door — the R-D4 entanglement, recorded). Mid-band rows are 0–1-atom
counts plus clock-blind ambient-spike evap, ratio ≈ 1: below the
grain floor, no claim. The ablation confirms R-D3's design content:
riding the raw slowing beat suppresses exactly the high-x emission
the stabilizer needs most.

### 4.2b Comb-blindness (P-R8) and the fifth "peak" — coarse verdicts

Outtake through the fifth and octave is smooth to the available
grain: k005 out at x = 0.75 / 0.8333 / 0.92 = 0.668 / 0.578 / 0.586
(no octave structure beyond ±1 atom); at the fifth the demand is
sub-atom by design. P-R8 stands at coarse resolution. Conversely the
v90 "suggestive intake maximum at the fifth" did NOT survive the
pinned pair: in(0.4167) = 0.301 (1 atom) vs neighbours 0.515/0.193 —
the peak was count noise on a wobbling pair. (The honest statement:
at T=120, in(x) has no resolved structure anywhere.)

### 4.2 The k-ladder converges on the quarter-power law (coarse, T=120, seed 20260802)

Endpoint Em_fin clusters of the rows seeded at/above the predicted
fixed point, per arm, against the closed-form Em* at the measured
same-run intake:

| arm | Em* predicted | Em_fin cluster measured | note |
|---|---|---|---|
| k002 | 1.97 | 1.89–1.96 (x0 ≥ 0.92 rows) | from-above return |
| k005 | 1.61 | 1.4–1.8 (x0 0.52–0.92) | two-sided drift toward |
| k010 | 1.34 | 1.0–1.52 (x0 ≥ 0.8333) | from-above return |
| k030 | 1.02 (≈4 atoms deep) | 0.13–1.18 scattered — grain-rattled floor | fixed point NOT resolvable |

**Correction (2026-08-04):** an earlier draft of this table recorded
a k030 "from-below convergence 0.875 → 1.16" — that number came from
harvesting an INCOMPLETE log (stale extractor fields). Withdrawn.
The extractor now refuses logs without RESULT lines. The complete
k030 arm shows hard draining everywhere and endpoint rattle of ±2
atoms: at k=0.3 the atom ε is 25–100% of the residual holdings and
the fixed point (≈4 atoms deep) is not resolvable at T=120.

The interior fixed point exists (P-R4 ✓ by convergence at
k002/k005/k010), is approached from both sides at k005 (drift toward
1.4–1.8 from x0 above and below), and its location tracks
x* = (in/(0.497·k))^{1/4} over k ∈ [0.02, 0.10]. At k=0.30 the
grain kick ε/Em* ≈ 25% destroys resolvability — a new,
pre-registerable selection criterion: **the fixed point must be deep
in atoms** (Em*/ε: k002 ≈ 6, k005 ≈ 5.5, k010 ≈ 4.8, k030 ≈ 4 but
rattle-dominated). Single-seed, endpoint-based, partner-blurred —
the gated numbers come from R1b (§1b). Selection rule (i)+(iii)
plus atom-depth points at (k_rad, p_rad) = (0.05, 4), pending the
ablation and steepness arms.
