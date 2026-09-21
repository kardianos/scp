# v97 consult — grok (2026-08-13)

Validate the fresh-eyes findings, rate the proposed next steps, add
ideas. Grounded in `v97/{README,PRINCIPLES,METRICS,H0,H1,LEADS}.md`,
`hunt/convert.go`, `fab/promo.go`, `fab/fcs.go`, `fab/step.go`,
`runs/leads/{H1.out,*.log}`, `runs/leads.sh` (`snap_every=200`,
`dt=0.02` ⇒ Δt = 4 t.u.), and `v93/theory/qgraph/docs/UNIFIED_LAW.md`.

C-kernel parity out of scope. Go is the substrate. No face authorized
by this page.

---

## Verdict (read this first)

1. **F1 is half-right about J, backwards about Φ.** The exact door
   ledger *is* in the FCS and *is* parsed. Hunt uses it only on snaps
   that have nonzero `promo_i`/`promo_o`; grainy atoms leave many
   frames at zero and those fall back to the heuristic. I measured
 The findings are not a package deal. F1 is backwards on Φ, F2’s algebra is fine and its H1 claim is not, F3 is the right lever with the wrong mechanism, and F4 is an event diagram rather than an H1 proof. Full review is also in `v97/consult/REVIEW_grok_next.md`.

---

## Verdict

1. **F1 is half-right about J, backwards about Φ.** The exact door ledger is in the FCS and is parsed. Hunt uses it only on snaps with nonzero `promo_i`/`promo_o`; grainy atoms leave many frames at zero and those fall back to the heuristic. I dumped `p2_space.fcs`: `totPromoI = 0.524917 = cond`, `totPromoO = 14.185 = evap+rad`. Graded ⟨I⟩=0.36 cannot be that ledger (exact I is ~0.009/frame). The heuristic **inflates** Φ by counting hops and by making I and O look simultaneous. The 10× shortfall is real and would get *worse* on an always-exact meter.

2. **F2’s algebra is correct; the identification is not.** `0.05 x³`, bar 0.015 ⇒ x ≥ 0.669, g ≳ 0.28 from Theorem 3.5 — all fine as a statement about the qgraph **per-cell rate** Φ\*. H1 Φ is a **snap-integral over a motif**. At Δt=4, even the adopted-law `comp = w2/w2e = 1+q_detune x` already puts the g=0 ceiling **above** the bar on a small loaded motif. H1 is not secretly “g ≳ 0.3.” H2 (restoring x\* > xb) is.

3. **F3 is the right lever, restated.** Large inclusive motifs dilute Φ/Em — not because radiance is a surface flux (it is a volume tax) but because the census drags lukewarm cells into the denominator and intake is field-starved. C6 is the right *size*. It is the wrong *topology* for the J clause (a loop fails `Jdip ≥ 0.25 Circ` by construction).

4. **F4 is true at the event level and oversold as an H1 proof.** `promo=2` does move Es pull/return onto oriented neighbours. The H1 Es meter is pass-S `CoreFlux`, and **p0_ctl also PASSed Es**.

5. **Order 1→4 is almost right and too bundled.** Fix the opportunistic pI/pO fallback and re-grade *before* any new arm. Then run the registered A.5 lit-bath ring as a small config stack, not a six-knob helix manifesto. Formalize the *corrected* Φ theorem, not F2 as written. A latch is a real face; it does not mint intake.

---

## F1 — meters

### F1(a) Φ

The review says hunt “never reads” `pI`/`pO`. That sentence is false as code and true as a description of most graded frames.

```123:150:v97/hunt/convert.go
	usePromo := false
	for _, idx := range m.Cells {
		if ... && (cur.Cells[idx].PromoI > 0 || cur.Cells[idx].PromoO > 0) {
			usePromo = true
			break
		}
	}
	...
		if usePromo {
			I, O = c.PromoI, c.PromoO
		} else {
			I, O = cellIO(prev.Cells[pi], c)
		}
```

`fab/fcs.go` writes those columns only when `promo>0` and resets them every snap. Through `hunt.Open` on `p2_space.fcs`:

| | value |
|---|---|
| SCHM after first frame | `… th2 promo_i promo_o` |
| Σ promo_i over 61 frames | **0.524917** (= `RESULT cond`) |
| Σ promo_o over 61 frames | **14.185** (= evap+rad) |
| cells with any I | 3 |
| cells with any O | 47 |
| t=4, 8, 60, 240 | I=O=0 exactly (grain) |

The ledger is honest, complete, and grainy. Hunt switches to it only when some motif cell is nonzero this snap; otherwise it scores hops as condensation/evaporation.

| arm | exact I/frame (cond/60) | graded ⟨I⟩ | exact O/frame | graded ⟨O⟩ |
|---|---|---|---|---|
| p2_space | 0.0087 | **0.36** | 0.236 | 0.26 |
| p0_ctl (no columns) | 0.023 | **0.58** | 0.247 | 0.34 |
| p2_hot | 0.348 | 0.27 | 1.83 | 0.59 |

p2_space ⟨I⟩ is ~40× the door. p2_hot is the exception because cond is large enough that ledger and hops are the same order.

The “10× may be fake” claim runs the **wrong error direction**:

1. Net-vs-gross *can* hide two-way traffic inside a cell. On the leads that is not the dominant term — cond is 0.53 over the whole run; there is almost no hidden gross I.
2. Hops look like I/O. That **adds** Φ.
3. Exact I and O almost never occupy the same snap. The heuristic makes them overlap. H1 grades ⟨min(I,O)⟩. The heuristic inflates it.

An always-exact re-grade of p2_space would look like: ⟨I⟩~0.01, ⟨O⟩~0.24, ⟨Φ⟩≲0.01, Φ/Em ~ 1e-4, and a screaming I/O imbalance (cond 0.53 vs rad+evap 14.2). That is SLOW_DRAIN, which H2 already refuses. The leads’ “almost balanced” ⟨I⟩≈⟨O⟩ is a meter fiction.

Archive EMPTY survives. Sep 2/106 is a real geometric fact. The quantitative leads Φ/Em deserves trust as an **upper bound**, not as door throughput.

### F1(b) J

Accepted, with the caveat they already wrote. H0 current SKIP is this signal (e3b ⟨TJnetX⟩ = −2 ± 14; `fd` is +x and not in the FCS). `promoPick` / `cellDir` steer by the same snapshot bilinear. dir-stable 0.54–0.69 is a mixture of meter slosh and real re-picking. Write a physics-inert hop-vector accumulator and H0 it against e3b `fd` before it is allowed to grade H1.

---

## F2 — the math (this is Q1)

### What is sound

```
O = k cap x^p ,     Em = x cap ,     p=4, k=0.05
O/Em = k x^{p−1} = 0.05 x³
```

0.05 × 0.62³ = 0.01192. Bar 0.015 ⇒ x ≥ 0.6694. u = 0.6694/0.62 = 1.0797. Theorem 3.5: g = (u⁴ − 1)/u³ ≈ 0.283. qgraph Φ\* at g=0 is k cap xb⁴ = 0.0185 — a *rate*, not a ratio. The review’s ratio algebra is the right transform of that rate.

Theorem 3.4 (g=0 ⇒ x\* = xb) is intact. `rad_dual` as coded is the spontaneous term only (`step.go`: `k * rad_dual * cap * xb^p * comp`). p2_dual not lifting Φ is consistent (cond actually *fell*: 0.16 vs p2_space 0.53).

Φ = min(I,O) equals O at balance **when I = O**. The leads are not at balance on the ledger.

### What is not sound

**Unit mismatch.** H1 Φ is `mean_frames min(I,O)` where I,O are per-snap motif sums (`h1PhiOverEm = 0.015`). Snaps are 4 t.u. (`leads.sh`: `snap_every=200`, `dt=0.02`). qgraph Φ\* is energy / time / cell. Identifying them smuggles a factor Δt.

If every cell of a motif sits at load x and I=O:

```
Φ_frame / Em = k x^{p−1} Δt  =  0.2 x³     (bare, Δt=4)
```

At x=0.62 that is **0.0476**, already 3× the bar.

**Adopted-law `comp`.** The kernel is not the qgraph coarse form. `rad_clock=0` ⇒ `comp = w2/w2e = 1 + q_detune x` on **both** radiance and `rad_dual`. Same multiplier, so Theorem 3.4 still gives x\*=xb. Throughput is larger:

```
O/Em|_rate = 0.05 x³ (1 + 1.2 x)
```

At x=0.62: 0.0119 × 1.744 = **0.0208 / t.u. > 0.015.** Even *in F2’s own rate units*, with the law as coded, g=0 clears the bar. The review cited “the adopted law” and then dropped `comp`. There is **no I/O asymmetry** in `comp` as coded. A complaint that intake *ought* to use the field clock is a different claim and a kernel change.

So: bar 0.015 as a *rate* is passable at g=0 under V3a. Bar 0.015 as a *frame ratio* is passable at g=0 for any small motif with x ≳ 0.37 and I≈O. The g ≳ 0.3 threshold is what you need to hold x\* **above ambient** (Theorem 3.5, H2), not to inhabit H1.

**Per-cell vs per-motif.** Cancels only if every live cell converts equally. Census floor `Em>0.05` pulls in a halo. Radiance is ~x⁴; halo cells add Em and almost no O. That is dilution by **inclusive census**, not a surface/volume theorem on the tax. Radiance runs on every `Em>0` cell (volume). Intake needs field; opacity can make *I* surface-limited. Then Φ ≈ I is surface and Em is volume. That is F3, and it is why the ~100-Em leads cannot pass even if you magically set I=⟨O⟩_exact: 0.24/95 ≈ 0.0025, still 6× short — mostly halo + low ⟨x⁴⟩, not missing g.

**Consequences, restated**

- “Promo Φ failure was predetermined / no door geometry can fix Φ.” **False as stated.** Geometry cannot invent field. It *can* raise ⟨min(I,O)⟩ toward min(⟨I⟩,⟨O⟩) by pairing events, and it can stop the census being a 100-Em foam. Mean I≪O on the ledger *is* predetermined: those arms ran as draining blobs in a random bath.
- “First face = slot door covers only half.” **Keep this.** Slot/orientation bought sep. It did not buy matched intake.
- Small loaded C6 with live doors *can* pass Φ/Em at g=0. It still fails the J clause if it remains a loop, and it still fails Φ if it stays at seed x=0.325 in the dark.

### p2_ring, since F2 leans on it

`p2_ring.log`: bath=0, wf_on=0, x̄=0.325, cond=0, **rad=0**, qcnvF=0.389 ≈ one leftover atom, Jdip/Circ = 0.116/5.68. Three independent reasons it is 1/5:

1. Dark vacuum → no field → no cond (A.5’s diagnosis).
2. Atom grain at this x → radiance never reaches ε (the leftover credit is the exhibit).
3. Loop topology → `notLoop` fails even if the door opens.

“Failed *only* because the vacuum was dark” is the slogan that makes Step 2 look sufficient. It is not.

---

## F3 — size

Keep the conclusion, change the mechanism.

- Φ/Em on a *uniformly converting* object is independent of N.
- Φ/Em on a *census* of a fluffy blob is small because ⟨x^{p−1}⟩ is small and the halo is in the denominator.
- Φ/Em on an *opaque* blob is small because I is surface-fed.
- C6 is the right size **for the Φ bar**.
- C6 is the wrong topology **for the J bar** unless you break the loop. Polarized I/O on a ring still has Circ ≫ Jdip. That is H0 ring, not H1.

Do not raise `EmFloor` until Φ/Em passes. That is cropping the bag. Report Φ/Em versus floor as a diagnostic.

---

## F4 — space-leg

`applyPromoCond`: dest = `promoPick(i, true)`; Es pulled from dest when `promo≥2`; Em lands on dest; `pI[dest]+=d1`.  
`applyPromoOut`: dest = `promoPick(i, false)`; `Em[i]-=d2`; field and back-share land on dest; `pO[dest]+=d2`.

Yes: Es is consumed at the fore neighbour and returned at the aft neighbour. Caveats:

- Only the `s_pull` share (0.15), not the converted energy.
- Fore and aft are **re-picked every event** from snapshot J. It is not one machine cycling.
- H1 Es is `CoreEsThru` from pass S, not the promo Es transfers. p0_ctl PASSed Es without any space-leg. p2_space is where sep+Es first co-occur. That is a real lead. It is not a proof that W-A is closed.

---

## Q1

F2’s derivation is sound as qgraph-rate algebra and **unsound as an H1 theorem**.

- Per-cell vs per-motif: cancels only for a uniformly converting motif. Inclusive census / opacity make the bar *harder* than the per-cell number — that is F3, not a +Δg on 0.3.
- `comp = w2/w2e`: same factor on I and O as coded; no M3 asymmetry in the ratio; *raises* the g=0 ceiling through 0.015 in rate units. Dropping it is how F2 got “unpassable.”
- Φ=min(I,O) equals O at balance iff I=O. Ledger I ≪ O on the leads; the heuristic hid that.

The theorem you actually want:

> Under V3a (`comp=1+q_detune x`) a uniformly converting motif with I=O and snap Δt satisfies `Φ_frame/Em = k x^{p−1} Δt (1+q_detune x)`. At the defaults (k=0.05, Δt=4, q=1.2) this is `0.2 x³(1+1.2x)`, which exceeds 0.015 for x ≳ 0.37. Theorem 3.4 still forbids a *restoring* root above xb at g=0. That is an H2 statement.

Do not write F2 into qgraph as stated.

---

## Q2

Net-vs-gross cannot account for the 10×. It runs the wrong way. Re-grade **before** Step 2.

The re-grade is not “teach hunt to read pI/pO.” That path exists. The re-grade is: **if SCHM has `promo_*`, always use them, including zeros.** Add an H0 class “exact vs heuristic on p2_space / p0_ctl.” Pre-register: exact Φ/Em drops, I/O imbalance appears, sep should still pass on promo>0. If sep dies on the exact ledger, the sep clause was hop-shaped too.

Optional `pJ`: yes, but it is a `fab/` write. H0 it on e3b against `fd` *before* any H1 J regrade.

---

## Q3 — the 1→4 plan

| step | verdict | note |
|---|---|---|
| 1 fix meters, re-grade | **do first** | recast as above; not a new reader |
| 2 lit-bath helix ring | **cut the bundle; keep A.5** | see below |
| 3a F2 as theorem | **cut / replace** | write the corrected theorem |
| 3b edge-door M_ab | **keep, cheap, do anytime** | independent of H1 |
| 3c graph-selection toy | **defer** | P3’s claim, not the H1 blocker |
| 4 latch face | **only after 2a/2b fail honestly** | do not claim it fixes Φ mean |

**Self-deception watch**

- “Sharpen meters” that make Φ larger would be a soften. Exact pI/pO will make Φ smaller. That is the honest direction.
- Grading J on `pJ` without an e3b H0 is shopping for a pass.
- `bath_xb` is an energy hose (retargets `Ee` on vacuum cells). Using it to inhabit H1 is filling the tank.
- `sw_kx`/`sw_amp` is a **seed** (`run.go`), not a maintained organized field. T≥1000 on a seed SW is a dark ring with extra steps.
- Six knobs at once (`amp_nat`, `amp_logate`, `amp_door`, `promo=2`, `rad_dual`, SW, helix, maybe `bath_xb`) violates “do not open several at once” and makes any PASS un-attributable.
- Raising `EmFloor` or tagging away the halo until Φ/Em ≥ 0.015 is cropping.
- Softening 0.015 or 0.80 is forbidden; the risk is *equivalent* softens via meter shopping.

**What I would run instead**

0. **Step 1′.** Always-exact pI/pO when the columns exist. H0 class exact-vs-heuristic. Re-grade the ten leads. Report Φ **also** as a rate (÷Δt) so the bar’s snap-dependence is visible. Do **not** change 0.015 without the user.

1. **Step 2a — the registered A.5 arm, one seed, T=240.** `exp=ring nv=6 seed_mw=2` + `laws_V3a` + `promo=2` + **bath=1** (or a field-only foam) + `wf_on=1`. No amp_door, no rad_dual, no SW, no helix. Grade `-h1` and `-tag` (tagged cycle, so a bath cannot merge the motif into the foam — the instrument already exists). Binary question: does the door fire, and does the motif stay C6?

2. **Step 2b.** If 2a is dark or grain-blocked: same arm + `amp_door=1` (the A.5 “hoped-for metabolic refill” that `RING_DNLS.md` recorded as untested). Still no rad_dual.

3. **Step 2c.** If doors fire and Φ is live but Jdip/Circ fails: add the kx / helix tilt, *alone*.

4. **Step 3b** in parallel. Edge-supported door on C_N; predict Δm = ±1. Replace 3a with the corrected Φ theorem.

5. **Step 4** only if 2a–2c produce live I and O that will not pair or will not pin a direction. Latch is a locus memory, not an intake law.

Missing and load-bearing: use `-tag`; diagnose p2_ring rad=0 with a one-cell grain check; keep H1 off H2; do not adopt `promo`.

---

## Q4 — resonance machinery vs a g-term

**Prediction:** the existing door can deliver *H1-scale* Φ on a small loaded motif from a **maintained** pitch-matched field (local xb high on the ring, box xb low). It cannot deliver a **restoring** x\* > xb. That requires the stimulated cross term (Theorem 3.5), which is not coded. v92 closed *want-side* feedback; a phase-matched door term is a different object, and still a law-table change. Do not sneak it in to buy H1.

| knob | what it does | g? |
|---|---|---|
| G_res / comb / p_gate | hop permission | no |
| workfn | existence threshold | no |
| cond `f_conv (Ee − e_cond)` | surplus-field fraction on the beat | raises xb, not g |
| `amp_door` | writes field phase into th2 | organization; untested on an organized field |
| `rad_dual` | spontaneous I = k cap xb^p | the g=0 term |
| `noise_amp` | random xb | g=0 by definition |
| `sw_kx/sw_amp` | **seed** | dies |
| `bath_xb` | apparatus hose | not a law |

amp_door+noise ≈ 1.00 is exactly the g=0 prior. It does not constrain amp_door+organized-field.

**Fastest deciding arm** is not the Step 2 manifesto. C6 + promo=2 + V3a + a bath that puts field on the ring without painting Em on the foam + `amp_door=1` + **no** rad_dual + **no** bath_xb + T=240, 1 seed, grade `-tag` with always-exact pI/pO.

- Door still dark → grain / workfn / no field. Apparatus, not g.
- Door live, I≪O, Φ/Em short, motif small and loaded → then the g-term (or a lawful maintained drive) is the hypothesis.
- Door live, I≈O, Φ/Em passes, Jdip/Circ fails → helix / latch, not g.
- Door live, all five pass on one seed → stop and replicate on 3 seeds before anyone mentions a face.

A maintained organized field is not in the repo as a *lawful* instrument. Inventing a continuous SW drive to pass H1 is the same species of cheat as `bath_xb`. The honest lit bath is field in the foam with workfn keeping vacuum dark. That is A.5, as written.

---

## Q5 — own ideas

### (a) Pinning J, besides the one-bit latch

The latch is the right *shape* of answer: stop re-picking max snapshot J every event. It is not the only shape, and it has a named failure (it can lock a drain).

1. **Slow compass, not a bit.** EMA of `cellDir` with a time constant of a few beats (FLOW / ASYM τ, applied to the *door steer*, not to hop conductance). Re-picking is what sloshes; a τ of one beat is already most of the win.

2. **The signed vector-credit v95 explicitly did not take.** `DESIGN_claude.md` A2: fitness was `|Ĵ|`; signed credit was declined because it would help backward slosh as much as forward *transport*. On the **door locus** that objection flips: you *want* a sign. A per-link running integral of signed promotion, decaying on the FLOW τ, is the latch with a real number instead of a bit. It does not amplify hop leaks (the p2_jm melt). **This is my preferred face** if 2a–2c fail.

3. **Phase-locked pair.** Two cells whose Δθ sits on the conversion beat; I only at one end, O only at the other, of *that* pair. P3 stated in clocks rather than in J. Harder; closer to the principle.

4. **Steer by integrated hops (`pJ`), not by bilinear J.** First as a meter. Then, if H0-green, as the promo compass. Use the signal H0 said was real.

5. **Geometry pin, config-only.** A helix (or an open current) forces an axis so the J clause is not asking a loop to be a dipole. Do this *before* any latch. If a helix with live doors still sloshes through 180°, the compass is the bug. If it holds, you did not need a memory.

6. **What I would not do.** `j_metric` on hops (measured melt). One-bit latch without a decay (scarred drain). Steering by the same snapshot you just SKIPped.

### (b) Other routes to Φ ≥ bar

The plan’s coherent-intake story is aimed at the wrong theorem (H2). Cheaper routes:

1. **Small motif, exact meter, doors actually firing.** C6 at x ≳ 0.4 with I≈O passes the frame bar at g=0. Config + grain, not a g-term.
2. **Pair I and O in time.** ⟨min(I,O)⟩ ≪ min(⟨I⟩,⟨O⟩) on a grainy door. A matched cycle raises graded Φ without more traffic. This is what Step 4 gets right, and F2 gets wrong.
3. **Stop grading the foam.** `-tag` on the painted cycle. Diagnostic first; do not retune the floor.
4. **Load the ring toward xb without a hose.** Field in the medium, workfn dark on Em=0, cond only on the cycle. If the ring sits below grain-threshold x, a little intake is what lets radiance start (p2_ring’s leftover qcnvF is waiting for one more push).
5. **Do not use rad_dual to “lift Φ” until I is the bottleneck on an honest meter.** p2_dual already lowered cond.

What will *not* work: more isotropic traffic (noise, hot amp), `j_metric`, tightening promo=3.

### (c) Pushback on the v97 charter

Not walking back the six principles or the dipole scoreboard. Three places the charter now over-claims:

1. **“Empty ⇒ first face is the slot door.”** EMPTY is right. The jump is too clean. Promo already *is* a slot-ish orientation and it bought sep+Es without inhabiting. P3-as-orientation is not sufficient. P3-as-“promote iff a cycle closes” is a different face and still does not mint field. After the leads, **intake pairing** belongs on the face list as its own line.

2. **Φ/Em = 0.015 is snap-dependent.** Set so Face C (~0.0025) and HOLD (~7e-5) cannot sneak through — right spirit. It is a frame ratio at whatever `snap_every` the stream used. A stream snapped 4× finer fails 4× more easily. State the bar with Δt, or as a rate. That is a sharpen of the *statement*, not of the number. Do not move 0.015 in this consult.

3. **⟨min(I,O)⟩ conflates mean imbalance with asynchrony.** Both are failures; they point at different faces (intake law vs pairing). Report both ⟨min⟩ and min(⟨I⟩,⟨O⟩) and the I/O ratio. The bar can stay on ⟨min⟩.

Smaller: do not let F4 + a green Es clause on p0_ctl talk us into “P5 is already implemented.” And do not let a g-term enter as “just the Φ bar in disguise.” It is H2, it is qgraph L4, and it is a law change.

---

## Answers, one line each

- **Q1.** Algebra yes; H1 identification no; `comp` makes g=0 *easier* than F2 says; S/V/census makes large motifs harder without moving g; Φ=O only at I=O, which the leads are not.
- **Q2.** F1’s Φ conclusion is backwards (heuristic inflates); J caveat is right. Re-grade before any new arm.
- **Q3.** 1′ → 2a A.5 → 2b amp_door → 2c helix → 3b / corrected 3a → 4 only if pairing/direction still fail. Cut bundled Step 2, cut F2-as-theorem, defer 3c. Watch `bath_xb`, seed SW, meter shopping.
- **Q4.** Existing door can inhabit H1 on a small lit motif; g-term is H2. Decide with A.5 + amp_door + exact meters + `-tag`, T=240, one seed.
- **Q5.** Prefer slow/signed door-locus credit over a one-bit latch; pin by geometry first; pair I/O in time; grade the tagged C6 not the foam; push back on “slot door is the whole first face” and on a snap-dependent Φ bar.

No face is authorized by this page. Config-only next work, if you want it, is Step 1′ (always-exact pI/pO + H0 class) then the single A.5 arm.
ric` on hops (measured melt).
   One-bit latch without a decay (scarred drain). Steering by
   the same snapshot you just SKIPped.

### (b) Other routes to Φ ≥ bar

The plan’s coherent-intake story is one route and it is aimed
at the wrong theorem (H2). Other routes, cheaper:

1. **Small motif, exact meter, doors actually firing.** C6 at
   x ≳ 0.4 with I≈O passes the frame bar at g=0. That is a
   config + grain problem, not a g-term.
2. **Pair I and O in time.** ⟨min(I,O)⟩ ≪ min(⟨I⟩,⟨O⟩) on a
   grainy door. A matched cycle (latch / pair / one atom in
   and out on the same beat) raises graded Φ without more
   traffic. This is what Step 4 gets right, and F2 gets wrong
   when it says geometry cannot move Φ.
3. **Stop grading the foam.** `-tag` on the painted cycle.
   Halo out of the denominator. Diagnostic first; do not
   retune the floor.
4. **Load the ring toward xb without a hose.** Field in the
   medium, workfn dark on Em=0, cond only on the cycle. A.5.
   If the ring sits below grain-threshold x, a little intake
   is what lets radiance start (the leftover qcnvF on p2_ring
   is waiting for one more push).
5. **Do not use rad_dual to “lift Φ” until I is the bottleneck
   on an honest meter.** p2_dual already lowered cond.

What will *not* work: more isotropic traffic (noise, hot amp),
`j_metric`, tightening promo=3 (already starved cond).

### (c) Pushback on the v97 charter

I am not walking back the six principles or the dipole
scoreboard. Three places the charter now over-claims, given
the leads:

1. **“Empty ⇒ first face is the slot door.”** EMPTY is right.
   The jump is too clean. Promo already *is* a slot-ish
   orientation, Go-only, and it bought sep+Es without
   inhabiting. P3-as-orientation is not sufficient.
   P3-as-“promote iff a cycle closes” (latch / pair / closure
   rule) is a different face and still does not mint field.
   The charter’s own table already split slot / promotion
   write / space-leg. After the leads, **intake pairing**
   belongs on that list as its own line, not as a footnote
   to the slot.

2. **Φ/Em = 0.015 is snap-dependent.** It was set so Face C
   (~0.0025) and HOLD (~7e-5) cannot sneak through, which is
   the right *spirit*. It is a frame ratio at whatever
   `snap_every` the stream used. A stream snapped 4× finer
   fails 4× more easily. The bar should be stated with Δt, or
   as a rate. That is a sharpen of the *statement*, not of
   the number. Do not move 0.015 in this consult.

3. **⟨min(I,O)⟩ conflates mean imbalance with asynchrony.**
   Both are failures, but they point at different faces
   (intake law vs pairing). Report both ⟨min⟩ and
   min(⟨I⟩,⟨O⟩) and the I/O ratio. The bar can stay on ⟨min⟩.

Smaller: do not let F4 + a green Es clause on p0_ctl talk us
into “P5 is already implemented.” And do not let a g-term
enter as “just the Φ bar in disguise.” It is H2, it is
qgraph L4, and it is a law change.

---

## Answers, one line each

- **Q1.** Algebra yes; H1 identification no; `comp` makes g=0
  *easier* than F2 says; S/V/census makes large motifs harder
  without moving g; Φ=O only at I=O, which the leads are not.
- **Q2.** F1’s Φ conclusion is backwards (heuristic inflates);
  J caveat is right. Re-grade before any new arm.
- **Q3.** 1′ → 2a A.5 → 2b amp_door → 2c helix → 3b / corrected
  3a → 4 only if pairing/direction still fail. Cut bundled
  Step 2, cut F2-as-theorem, defer 3c. Watch bath_xb, seed SW,
  meter shopping.
- **Q4.** Existing door can inhabit H1 on a small lit motif;
  g-term is H2. Decide with A.5 + amp_door + exact meters +
  `-tag`, T=240, one seed.
- **Q5.** Prefer slow/signed door-locus credit over a one-bit
  latch; pin by geometry first; pair I/O in time; grade the
  tagged C6 not the foam; push back on “slot door is the whole
  first face” and on a snap-dependent Φ bar.
