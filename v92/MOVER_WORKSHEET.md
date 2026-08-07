# MOVER WORKSHEET — the V3a per-bar audit (v92 act one; decisions are the user's)

**Method.** Full battery run at the V3a table (`-extra "k_rad=0.05 p_rad=4
rad_clock=0 wf_on=1 wf_floor=0.01 wf_far=99"`, both kernels;
`runs/BATTERY_v3a.log`). Baseline = STEP ZERO V2g run (`runs/BATTERY.log`,
ALL GREEN 93). Attribution by two single-knob scans: radiance-only
(`runs/v3a_r/`) and workfn-only (`runs/v3a_w/`). W-L2 prerequisite PASSED
(8/8 ds1 bars under wf; `runs/BATTERY_wl2.log`). Nothing here is a
decision — under the ratchet every moved bar gets an explicit
sharpen / re-gauge / reject call with the user.

## Result: 12/93 moved (81 green)

| # | bar | V2g | V3a | attrib | class (recommended) |
|---|---|---|---|---|---|
| 1 | conserve: c law purity (header) | byte-exact | header moved | both | **RE-GAUGE** |
| 2 | conserve: go law purity (header) | byte-exact | header moved | both | **RE-GAUGE** |
| 3 | ring6 edge_dev_mean ≤ 0.05 | 0.0350 | 0.1533 | rad | **RE-GAUGE/RETIRE** |
| 4 | ring6 all 6 edges alive, min gg ≥ 0.9 | min 0.975 | min 0.000 | rad | **RE-GAUGE/RETIRE** |
| 5 | blob retention ≥ 0.78 @T160 | 0.805 | 0.449 | rad | **RE-GAUGE/RETIRE** |
| 6 | p1_b2drv \|v_COE\| ≤ 5e-5 | 2.82e-5 | 5.16e-5 | rad | **DECIDE (marginal)** |
| 7 | pauli0: admission == 0 at cap (id) | dep 0.0 | dep 0.31 | rad | **DECIDE (physics)** |
| 8 | pauli0: admission == 0 at cap (5th) | dep 0.0 | dep 0.018 | rad | **DECIDE (physics)** |
| 9 | xsec: absorption pure cond (evap==0) | 0/0 | 0/1.1 | rad | **RE-GAUGE** |
| 10 | xsec: optics null (net_tag=cond=0) | 0/0 | 6.6/11 | wf | **RE-GAUGE/RETIRE** |
| 11 | xsec: upbeam untouched (hdr) \|rE-1\|≤0.05 | 0.970 | 1.051 | wf | **DECIDE (marginal)** |
| 12 | abx: blob2+p1 BYTE-identical (incl drift) | byte-eq | drift 1.6e-16 | rad | **INVESTIGATE** |

## The canceling interaction (measured, supports adopt-together)

`xsec headroom absorbs (net_tag ∈ [6.9,7.65])` FAILS in BOTH single-knob
scans but PASSES in combined V3a: radiance-only **6.70** (below — radiance
outflow reduces net absorption), workfn-only **7.96** (above — wf increases
absorption), combined **6.98** (in band). The two opposing effects cancel
to land inside the window. This is the quantitative form of RECKONING §3's
"radiance × workfn: compatible; adopt-together is coherent" — measured, not
asserted. (Not a mover in the V3a table because it lands green; recorded
here as the interaction evidence.)

## Per-bar notes & the decisions needed

### Trivial re-gauge (1,2) — purity header, by design
The header now prints `k_rad=0.05` and `wf_on=1`; the purity pins re-point
to the V3a header strings. This IS the adoption event's mechanic. Recommend
RE-GAUGE (re-pin both purity bars to the V3a header; no physics).

### Radiance taxes V2g hoards (3,4,5) — the known R5 movers, by design
ring6 and blob are V2g-era **hoard-objects**; radiance taxes them (measured
t_half 80–140 vs 260–510, `v91/RADIANCE.md`). This is the POINT of radiance
— it removes un-radiated hoard stability and selects the interior balance
(x̂*=0.62). The stability target moves to radiance-stabilized objects
(chords, the balance-point bodies). Recommend RE-GAUGE/RETIRE: these bars
measured a V2g-only property (static hoard stability) that radiance
intentionally supersedes. The bars either retire or re-pin to "hoards taxed
as designed." **Decision: retire, re-pin, or keep-as-regression?**

### The physics decision — pauli0 at cap (7,8)
PAULI-0 ("admission exactly 0 at cap") breaks: dep 0 → 0.31 (id) / 0.018
(5th). Mechanism: at cap the cell now RADIATES (the new law), creating
inter-beat headroom, so the next beat re-admits to replace what radiance
removed. The cap is no longer a static hard wall — it is a radiating
equilibrium (admission = radiance throughput). The "near cap" bars still
PASS (throttled admission). Reading: physically correct (a saturated body
that radiates is Stefan-Boltzmann; XSEC already showed saturated objects
emit), but it changes the V2g saturation claim from "hard refusal" to
"radiate-and-refill equilibrium." **Decision: is this re-gauge (the cap
refusal becomes dynamic — consistent with the radiance law) or does it
weaken the exclusion precursor (B5) enough to reject/condition radiance?**
My read: RE-GAUGE — exclusion (B5) is the exchange-sign door (L-3), not
static cap refusal; the instantaneous PAULI-0 still holds at the beat. But
this is the one mover with real physics-sign and it is yours to call.

### Marginal decides (6,11)
- **p1_b2drv v_COE** 2.82e-5 → 5.16e-5 (3% over the bar). The driven
  two-blob COE drift rose slightly under radiance. Marginal. Recommend
  RE-GAUGE (loosen to 6e-5) or re-run multi-seed to confirm it's not a
  speckle blip. The p1_b2CTL twin still PASSES (3.71e-5).
- **xsec upbeam (hdr)** 0.970 → 1.051 (law-regime arm). Under wf the law-
  regime upbeam foam is transparent (wf_far) instead of condensing
  (e_cond=0), so the beam's upbeam field distribution shifts. The optics-
  regime twin still PASSES (1.0025). Recommend RE-GAUGE (the (hdr) arm is
  the retired law-regime; the optics twin carries the claim) — but confirm.

### Re-gauge — regime bars retire (9,10)
- **xsec pure cond (evap 0→1.1):** radiance adds an emission channel at the
  absorber. The absorption ITSELF still works (net_tag in band). Re-gauge:
  "absorption = cond + radiance outflow." Recommend RE-GAUGE.
- **xsec optics null (0/0 → 6.6/11):** under wf the dense object (slit_obj)
  converts even in "optics" mode — dense matter is the emergent detector.
  This IS the workfn retiring the optics/law declaration. Recommend RETIRE
  (the bar assumed e_cond=99 suppresses all conversion; that assumption is
  what workfn replaces).

### Re-gauge — C≡Go (12) — RESOLVED, benign
abx blob2+p1: drift column C `+1.62e-16` vs Go `+0.00e+00`; **every physics
column byte-identical**. Investigated: at V2g baseline the other four abx
pairs (ds/ds1/rings/xs) **already** had C≠Go drift (masked, accepted) while
blob2 alone matched (both 0). Under V3a radiance makes blob2 join the same
FP-floor divergence pattern — it is the **pre-existing C/Go accumulation-
order envelope** (IDENTITY §4.0), not a radiance-path asymmetry. Magnitudes
are machine-epsilon class (worst 1.17e-13 across all pairs). **RE-GAUGE:**
align blob2's abx bar with the other four ("identical up to drift col") —
the "incl drift" strictness was the outlier, not the physics. No A/B
covenant issue remains.

## Summary recommendation (for the adoption commit)

Of 12 movers: **9 are clean re-gauge/retire** (purity ×2, hoards ×3, regime
bars ×2, abx-drift ×1 — investigated, benign; + the canceling-interaction
evidence recorded). **2 are marginal decides** (p1_b2drv, upbeam-hdr — both
re-gauge candidates with the optics/law twin still green). **1 is a real
physics decision: pauli0 at cap (×2 bars)** — radiance turns the static
"hard refusal" into a radiate-and-refill equilibrium (dep 0→0.31/0.018; the
"near cap" bars still pass). None is a clean REJECT of either candidate:
radiance's movers are its designed effects (tax hoards, radiate at cap, add
emission); workfn's movers are the regime-declaration retiring (its designed
effect). The combined table is coherent (the canceling interaction lands
green).

**Proposed adoption gate:** (a) your call on **pauli0** (re-gauge: cap
refusal becomes the radiating equilibrium, consistent with the law — vs
condition/reject if you read it as weakening exclusion); (b) your call on
the two **marginals** (re-gauge both, or re-run multi-seed for p1_b2drv);
then the purity pins re-point to V3a, blob2's abx bar aligns with the
others, both kernels' headers update, C≡Go re-verified, and the adoption
commit lands `laws_V3a.cfg` as the table.
