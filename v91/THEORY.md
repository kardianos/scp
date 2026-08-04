# v91 THEORY — the law as it stands, and the two modifications under test

Written 2026-08-04 at v91 opening (user-authorized). Three parts:
§1 the standing law (carried, unchanged, with its evidence); §2 the two
v91 law modifications (RADIANCE — implemented, inert by default;
COHERENT CHANNEL — specified, not yet implemented); §3 the
back-reference table tying every load-bearing claim to its committed
evidence. The v90/v89 trees stay frozen as the record.

---

## 1. The standing law (laws_V2g, carried verbatim)

**The principle:** energy is never destroyed; it only changes mode;
space is one of the modes. Therefore no background, no imported field
or species, everything emergent from the fabric. One law table —
`laws/laws_V2g.cfg` (copy; canonical at `v89/battery/laws_V2g.cfg`) —
and the ratchet rule: every kernel or table change runs the FULL
battery; green tests gate all future changes.

The law, in its measured form (each clause battery-gated or logged):

* **Amplitudes within a mode; atoms at mode boundaries.** Transport is
  never quantized; every conversion fires in whole action atoms
  ε = A₀ω/2π at the source's pitch, two-atom credit. ħ-linearity
  ~1e-8; every detection click is k·ε, k ∈ {1,2}, at n=362.
* **Load flattens pitch; flight is load.** x = (Em + flload)/cap;
  w_e = w/(1 + q_detune·x), BOTH pitches. Consequence measured in
  v90: absorption and emission ride a clock that SLOWS with fullness.
* **Consonance gates dense transport.** The comb (p·q ≤ 6) decides
  which pitch pairs can exchange; gates are phase-closure conditions
  (cos^p_gate). Consequences: quantized bond ladder
  (q·wᵢ+p·wⱼ)d/C = 2πm; parity selection (even rings live, odd
  π-frustrated); the fifth as the charge-analog (FCQ); U/D
  molecule/droplet asymmetry.
* **Conversion is detection.** Field condenses at matter on the beat;
  matter evaporates above cap; roughness radiates on detune. XSEC
  triad: headroom absorbs (+7.27 clean), saturation emits (−7.2),
  inert load refracts/reflects.
* **Space flows; pressure pushes; nothing reaches out.** s_k/s_disp
  transport; the g1 footprint; emergent opacity; no steady 1/r
  monopole without an internal space cycle (g4).
* **Momentum is the first moment of conversion.** The all-modes COE
  meter: interaction is channels; blob approach carries no net
  momentum (≤1e-3 of closing).
* **The two-plane field sector is the wave sector.** Interference at
  parameter-free loci; Born-consistent clicks at ±2%; eraser, delayed
  choice, no-signaling (v89); band-like dispersion (measured limit).

**What this law cannot do — measured, not conjectured** (the v91
motivation): no interior stable mass (FORGE E1: outtake exactly 0
below cap ⇒ the only fixed point is the cap wall); no transport of
structure (Δp ≤ 1e-3 of closing); the derived rate-level corrections
kill the fifth (S2 entry: kappa_reac=1 collapses FCQ and composite
boundaries); no size selection, no spectra, no exclusion.

## 2. The v91 modifications

### 2.1 RADIANCE — laws_V3r candidate A (IMPLEMENTED, inert at defaults)

**The measured deficiency** (`carried/FORGE_EVIDENCE.md` E1): V2g's
emission-vs-fullness curve is a cap-wall step with anti-Stefan slope —
a loaded cell in a foreign-pitch bath emits EXACTLY nothing below cap
(at every fullness, including comb resonances), then a thin trickle
throttled by the slowing beat. Intake is positive everywhere and load
RAISES it. Hence in > out for all x < 1: everything grows, pins at the
wall, and becomes a sink (D-ring x̄→0.83) or an emitter-corpse (sat20).
Reality stabilizes bodies with radiance that RISES steeply with
fullness (Stefan–Boltzmann T⁴). The modification supplies exactly that
term and nothing else.

**Candidate A, as implemented** (both kernels, pass 6, after the evap
block; guarded by `k_rad`, default 0 ⇒ byte-identical to V2g —
verified by full battery and by C≡Go byte-diff):

```
on each beat of cell i, if k_rad > 0 and Em_i > 0:
    x      = Em_i / cap                  (holdings, not flight — R-D1)
    comp   = rad_clock ? 1 : w2/w2e_i    (= det_i: beat-stretch compensation)
    demand = k_rad · cap · x^p_rad · comp
    fire through atoms machinery (grains ε; two-atom credit, shared
        F-side credit register with evaporation)
    clamp to Em_i; route exactly as evaporation:
        Em_i −= dr;  field gets dr − bs;  space gets bs = dr·s_pull/(1+s_pull)
    ledgers: global rad_total (RESULT conv ... rad=), tag ledger folds
        into the evap column (# CONVTAG unchanged in format)
```

New law-candidate keys (printed on their own header line, pinned by
the battery's law-purity check at the inert defaults):
`k_rad` (strength; 0 = V2g), `p_rad` (steepness; default 4 — the
T⁴ analog, TO BE SELECTED BY MEASUREMENT, not assumed), `rad_clock`
(0 = per-time rate ∝ x^p_rad via det compensation; 1 = ride the raw
slowing beat — the ablation arm).

**Open decision points, deliberately left to the campaign:**
* **R-D1 (x base):** holdings Em/cap (implemented) vs pitch-x
  (Em+flload)/cap. Flight lives on links and cannot be radiated from
  the cell; but the PITCH rides flight — if lifetime experiments show
  bond flight escaping the balance, revisit.
* **R-D2 (p_rad):** the balance-curve experiment selects it. In a 2D
  substrate the naive mode-count analog is T³; in 3D T⁴. Treat p_rad
  as measured, not derived, until a mode-count derivation exists.
* **R-D3 (clock):** rad_clock=0 is the physical candidate (a fuller
  body radiates more per unit time); rad_clock=1 is the ablation that
  shows why compensation matters.
* **R-D4 (grain register):** radiance shares the F-side credit
  register with evaporation (implemented). If grain statistics
  entangle the two channels in a measurable way, split registers —
  battery-gated.

**Acceptance surface (from FORGE §4 and REALITY C):**
1. The balance curve gains an interior crossing: intake(x) = out(x)
   at some x* < 1 with d(out−in)/dx > 0 (stable side), for a band of
   ambients.
2. Perturbation return: a cell/blob displaced from x* returns
   (relaxation time measured), instead of running to the wall.
3. Lifetime: ring6/blob/composite outlive the v90 ceilings (ring
   death t≈1900, composite dissolution t≈480, 5k bound) by ≥ 10×,
   with intake = outtake measured in the ledger (the flux-machine
   criterion, now satisfiable).
4. Size selection: the stable x* exists ⇒ pitch w(x*) ⇒ predicted
   rung d* = πC/w(x*); measured bond lengths of surviving objects
   land on it. Sizes become DERIVED.
5. Forging closure: the FORGE bootstrap (matter→absorption→growth)
   now SATURATES at x* — condensation in a density well produces a
   finite object, not a sink.
6. Nothing standing breaks: full battery green at the adopted
   (k_rad, p_rad) — or bars are renegotiated explicitly by the
   ratchet rule (a law change MAY re-gauge bars; every change is a
   conscious decision, never a silent softening).
7. Spectra probe (opportunistic): with a non-slowing radiance clock,
   the comb-resonance emission structure (suggestive intake line at
   the fifth, FORGE E1) becomes measurable — the first emission-line
   candidate.

### 2.2 COHERENT DENSE CHANNEL — the S2 completion (SPECIFIED, not implemented)

Carried in full in `carried/S2_CHANNEL.md` (the v90 S2.md: entry
measurement + design memo). Summary of the specification:

* Dense transfers carry complex amplitude (magnitude AND phase),
  phases composing along paths — "translation IS the current".
* The p:q charts must transport amplitude with the interval's phase
  map (θ → (q/p)θ); the kr=1 sweep is the negative template (any
  design that kills the nv=6 boundary fifth or bloats the D-ring past
  x̄ 0.4 fails).
* Exchange registry: two-quantum amplitudes with a sign under
  exchange — the door to exclusion and full HOM depth.
* Flux-machine interior: persistent internal currents at zero net
  transport.
* Invariants that must NOT change: atoms machinery, conversion door,
  space sector, contact rule, and (new) the adopted radiance law.

Quantified acceptance (from S2.md §2): e3b speed ≥ 2.6e-3 with
cos ≥ 0.95 seed-robust; nv=6 boundary fifth gg ≥ 0.9 at T=100 with
D-ring x̄ ≤ 0.35; UUD ggm ≥ 0.9; a flavored composite RADIATES;
Y-junction occupancy statistics separate identical from consonant-
distinct arrivals.

**Order of work: RADIANCE FIRST.** It is small, decisive, independent,
and the stable objects it creates are the test-bodies every coherent-
channel experiment needs (a transport experiment on matter that
dissolves measures dissolution, not transport).

### 2.3 Third lane, design-only: cell number as a dynamical variable

FORGE 1a/E2b: the cell roster is the last background; bare density
texture cannot stand (voids heal ~15–40 t.u.) so dark-matter-like
standing texture and space-as-created-mode need cell fission/fusion.
BLOCKED until radiance stabilizes matter (fission bookkeeping needs
the radiance sink/source terms). Design constraints recorded in
FORGE_EVIDENCE.md §1a. No implementation in v91 until the R-lane
lands and the user re-authorizes.

## 3. Back-references (claim → evidence)

| claim | where proven | tree |
|---|---|---|
| law table + 20/20 unified pass | `v89/battery/laws_V2g.cfg`, `v89/README.md` | v89 (frozen) |
| free substrate carries interference (tier 0) | `v90/DS.md`, battery `ds` | v90 |
| clicks rebuild the wave, Born ±2%, atomicity | `v90/DS.md` tier 1, battery `ds1` | v90 |
| XSEC absorber/emitter/reflector; clock-rate absorption | `v90/MOTION.md` XSEC, battery `xsec` | v90 |
| momentum = first moment; merging not transport | `v90/MOTION.md`, battery `p1` | v90 |
| saturation refusal pitch-blind (exclusion control) | battery `pauli0` | v90 |
| composite bond/parity/U-D/dissolution record | `v90/COMPOSITE.md` | v90 |
| C↔Go A/B protocol + byte claims | `carried/VERIFY_PROTOCOL.md`, battery `abx` | copy here |
| local-clock scheduler (four conditions, pending-min rule) | `carried/P2_SCHEDULER.md`, battery `p2lc` | copy here |
| kr=1 kills the fifth (why rate-level fails) | `carried/S2_CHANNEL.md` §1, `v90/runs/s2/` | copy here |
| cap-wall step / anti-Stefan / no fixed point | `carried/FORGE_EVIDENCE.md` E1, `v90/runs/forge/` | copy here |
| dark-void lensing, exact-zero ledger | `carried/FORGE_EVIDENCE.md` E2 | copy here |
| forging follows density | `carried/FORGE_EVIDENCE.md` E3 | copy here |
| correlation audit + degree scale | `carried/REALITY_AUDIT.md` | copy here |
| pre-v89 ban | `CLAUDE.md` (unchanged, permanent) | root |
