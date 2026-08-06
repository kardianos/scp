# REALITY — the v90 correlation audit (retrospective, 2026-08-04)

What in this model currently corresponds to reality, **by how much**;
what does not, **to what degree**; and what to change so it corresponds
better. Every claim below is grounded in a gated bar or committed log
(battery 93 bars GREEN at time of writing; pointers per row). This is
an audit of the v89–v90 free-cell program under laws_V2g — not of any
earlier era.

**The degree scale used throughout:**

* **EXACT** — machine precision, by construction or measured.
* **QUANTITATIVE** — parameter-free numbers match within measured error.
* **STRUCTURAL** — the right phenomenology and mechanism-shape; no
  numeric correspondence to a real measured constant.
* **ABSENT** — the real phenomenon is missing from the model.
* **ANTI** — the model currently shows the opposite of reality.

---

## A. What correlates, and by how much

### A1. Energy conservation — EXACT
The one law. Ledger drift ≤ 2e-15 relative across every run of every
apparatus (~150+ runs this cycle), including channel birth/death, whole
conversion cascades, and the XSEC fog. Reality: exact (to all known
experiment). This is the program's axiom rather than a discovery — but
the fact that a background-free substrate with live geometry can keep
it at the FP floor through all of the above is the load-bearing
engineering result. *(battery: every drift bar)*

### A2. The double slit, tiers 0–1 — QUANTITATIVE (the strongest physics correlation)
* **Tier 0:** interference fringes at **parameter-free loci** — maxima
  and minima where δ(y)=mλ says, from geometry alone: y_peak within
  1.1 of prediction on every substrate seed, side maxima at the m=±1
  loci, additivity broken exactly where it must be (r 0.34–1.67).
  Visibility V_r 0.45–0.65 (foam-limited; ideal apparatus reaches ~1).
* **Tier 1:** single quantized detection events **rebuild the wave's
  statistics**: 362 whole-grain clicks over 24 substrate realizations,
  fringe-phase score s̄ = 0.6827 ± 0.0168 against the phase-blind null
  0.5 (10.9σ) — and equal to the wave-intensity expectation 0.683
  within its own SEM. **That is a Born-rule statement at the ±0.02
  level:** clicks sample the wave's intensity distribution. Which-path
  truncation removes the structure (controls 0.380, minima fill 7.7×).
* v89 tiers 2–3 (same substrate family, frozen graph): eraser
  fringes/anti-fringes, delayed choice, **no-signaling exact**.
* Degree summary: the core quantum double-slit phenomenology —
  wave loci parameter-free, particle statistics Born-consistent at 2%,
  complementarity, no-signaling — is REPRODUCED. What is not: absolute
  visibility (0.65 vs ~1) and 3D (the panel is 2D).
*(DS.md; battery `ds` + `ds1`, 16 bars)*

### A3. Quantization of exchange — QUANTITATIVE within a mode, ANTI across modes
* Conversions fire in whole atoms ε = A₀ω/2π; **every one of 362
  clicks is k·ε with k ∈ {1,2}, deviation 3.5e-6** (print floor);
  ħ-linearity of the atom with pitch measured to ~1e-8 (v89). This is
  E = nħω at the exchange boundary — the Planck relation, structurally
  and numerically clean within any one channel.
* **BUT ħ is not emergent as a universal constant: the per-process
  action spread is 19×** (v89 HBAR.md). Reality has one ħ across every
  process ever measured (to 1e-10). This is the sharpest ANTI in the
  quantum sector.
*(DS.md tier 1; v89 HBAR.md)*

### A4. Matter–field coupling (XSEC triad) — STRUCTURAL
Absorption, emission, and refraction all exist with the right
qualitative structure: a below-capacity object absorbs cleanly
(net +7.27, pure condensation, monotone impact-parameter profile), an
at-capacity object emits (net −7.2, evaporation side-glow 1.54×), an
inert object refracts/reflects (impedance defect, core rE 0.79, not
delay). "Detection is conversion" and "opacity is unfilled capacity"
parallel real absorber physics (ground-state atoms absorb, saturated
gain media emit, transparent media refract). Two honest caveats:
absorption here falls with load (clock-rate effect) — real
cross-sections are resonance-structured, not mass-monotone either, but
no numeric correspondence has been attempted; and the *law-regime
vacuum itself absorbs* (see B7). *(MOTION.md XSEC; battery `xsec`)*

### A5. Bound structure and its quantization — STRUCTURAL
* Bond lengths are **quantized by a wave condition**: the ladder
  (q·wᵢ + p·wⱼ)·d/C = 2πm, locked from both sides, offset =
  pressure/K_b exactly. Real chemical bonds are also standing-wave
  conditions — the FORM matches; no real molecule's numbers do.
* **Species selection by parity**: even rings live, odd rings are
  π-frustrated (ring5: ≥3 of 5 gates dead) — a structural echo of
  Hückel aromaticity (4n+2 rules) from pure phase closure.
* **U/D asymmetry**: D-flavor rings are droplets, U-flavor molecules;
  UUD closes into a triangle, UDD opens into a chain — proton/neutron-
  flavored structure from one incomplete harmonic (the fifth). The
  naming is suggestive; no masses or charges correspond numerically.
* Nucleon–nucleon residual force: contact-ranged (1–2 channels below
  sep 9, exactly inert at 11) — nuclear-force-shaped (short-ranged,
  saturating), again structurally only.
*(FREECELL §9–10; COMPOSITE.md; battery `pair`/`ring`/`fcq`)*

### A6. Statistics precursors — STRUCTURAL, partial
* **HOM exchange signs** (v89 tier 4): boson dip / fermion peak with
  the right signs, g_b 0.42 < 0.5 < g_f 0.58 — but mode-match-limited
  depth (real HOM: g_b → 0). Half-correlated.
* **Saturation refusal is distinguishability-blind** (PAULI-0, gated):
  at cap, admission is exactly 0 for identical AND consonant-distinct
  arrivals. This is the correct NEGATIVE control (fullness ≠
  exclusion); actual exclusion is ABSENT (B5).
*(v89 record; battery `pauli0`)*

### A7. Near-field gravity analogs — STRUCTURAL, weak
A mass maintains a graded space depression (g1); matter is emergently
opaque to space flux (g3); the vacuum skirt was a scored prediction;
lensing has the right sign but sits below the foam floor; free-fall is
below the chaos floor. The FAR field is absent (B6). *(v89 g-campaign)*

### A8. Determinism/causality infrastructure — EXACT (execution, not physics)
Byte-determinism at any thread count, C↔Go cross-implementation
identity at the physics-digit level, local-clock scheduling with
causal ordering (P2: batch ≡ serial exactly, tick skew 219 held by
local time). This correlates with reality's *bookkeeping* properties
(no frame-dependent outcomes of the bookkeeping kind), and is what
makes every claim above falsifiable to the byte. *(VERIFY.md, P2.md)*

---

## B. What does not correlate, and to what degree

### B1. No stable particle — ABSENT (the deepest gap)
Every bound object decays or dissolves: the ring6 dies at t≈1900 at
the vacuum-skirt load (v89 MASS); composites dissolve into the
protecting bath by t≈480 (CQ8-long, −36–39% load); the no-particle
bound under V2g is t ≲ cap·(x−x_skirt)/c₀ ≈ 5k. Reality: the electron
is stable to ≥ 1e28 yr. The gap is unbounded. The measured diagnosis
is specific: **a stationary pitch requires steady internal
circulation** (the flux-machine thesis, forced by five independent
measurements) — matter that keeps its pitch is matter that keeps its
books flowing, and the current dense sector cannot hold internal
currents.

### B2. No transport of structure, no momentum — ABSENT (measured to be zero)
Bound objects do not travel (wander ≤ 1.6% of their own radius, exact
null under drive in vacuum); blobs approach by profile-merging with
all-modes momentum ≤ 1e-3 of the closing rate (COE meter, gated);
radiation pressure shows a ~100× deficit at the conversion door (p2).
Reality: things move, ballistically, carrying momentum; light pushes.
This is not noise — the honest meters say the Δp of everything we can
build is zero. The e3b window (speed 2.6e-3, cos 0.95, one seed) is
the only forward crack, and it is drive-dependent and seed-fragile.

### B3. One ħ — ANTI (19×)
See A3. Within-channel exact; across processes the action atom spreads
19×. Reality: one constant to 1e-10.

### B4. Spectra of species — STRUCTURAL (re-graded 2026-08-06, user decision; was ABSENT)
**Re-grade executed on the v91 ASTRO/COMBINE evidence** (ASTRO.md §4.3,
§6.5.4-as-corrected, §6.5.5; COMBINE.md §5.2–§5.3). What is now
measured: bound species EMIT AND ABSORB AT LINES — a parameter-free
doublet per voice, w = (2.9 | 1.65)/(1 + 1.2·x_cell), per-event
residual mean 2e-4, rms 1.7% (warm) / 0.34% (cool), seed-robust,
load-tracking; a species is detectable at range by light alone (D line
30–60× over a spectrally-dark bath); the anti-Stokes gap sits at the
law constant w2/w1 to 0.3%; metabolism is spectrally visible
(eat-at-U / shine-at-D); death has a spectrogram (redward broadening
on partial death; a draining voice's line walks blueward). Species
identification runs on the LINE-RATIO PATTERN, not on shifted loci:
the healthy-D locus is a shared medium operating point (UUD
1.6710±0.0052, UDD 1.6774±0.0071; coexisting chords ±0.001), while
absorption:emission ratios differ 2.1–2.6× between species —
abundance-style spectroscopy. Brightness numbers per the dated
correction: UUD singles 0.0110/0.0150 D-line ev/t.u., UDD 0.0170
(sub-linear per emitter — intake-limited), two coexisting chords
0.0220 = 2.00× (additive by metabolism).
Why STRUCTURAL and not QUANTITATIVE: no correspondence is claimed to
any real atomic constant — the numbers are the model's own law read
back parameter-free. Standing edges, kept honest: the emitter is the
conversion door of bound matter under the RADIANCE candidate drive
(k_rad=0.05; at laws_V2g defaults the door is near-silent — COMBINE
CO-AR measured all luminosity as radiance-era, ~200× door collapse
without it); U lines are statistics-poor off intact topology and
camouflaged on the bath band; these are door-event spectra, not a
free-propagating EM wave spectrum (B7's regime split stands).
*(The pre-re-grade text recorded: consonant matter exactly dark at
detune 0; composite EMF gated on fifth transport; the mechanism slot
existed but nothing stable drove it — the frozen chord under radiance
is the stable driver it lacked.)*

### B5. No charge, no Coulomb, no exclusion — ABSENT
FCQ (charge as the incomplete harmonic) is a structural candidate
only. The space law pushes — nothing attracts at range, no 1/r force,
no sign structure. Exclusion is blocked on the exchange registry
(two-quantum amplitudes with a sign); PAULI-0 proved the current
refusal carries no identity information.

### B6. No far field — ABSENT
g4 closed the question honestly: a leaking blob's space flux is all
mass-rate bookkeeping; **no steady 1/r monopole exists** without a
stable particle's internal space cycle. Gravity's defining far-field
behavior awaits B1's resolution.

### B7. The law-regime vacuum is a cloud chamber — ANTI
With e_cond=0 (the table's value), the vacuum condenses any
above-threshold beam into a matter trail (XSEC: global cond ≈ 77 of a
~180-unit beam over 26 length units). The real vacuum transmits light
over Gpc. The program handles this today by *declaring* an optics
regime (q_detune=0, e_cond=99) per experiment — an override, not an
emergence. Conversely the optics regime's vacuum cannot detect
anything. Reality has both at once (transparent vacuum + absorbing
detectors) because conversion thresholds live in BOUND matter, not in
space.

### B8. Propagation is band-like and frame-fixed — ANTI for relativity
Field dispersion is lattice-band-shaped (v89 EM1 "band-like/massive";
measured here: 2D group speed 1.05·C at kx=0.9, 0.574·C in 3D at
kx=1.1, k-dependent). Reality's vacuum light is dispersionless and
frame-invariant; the foam is a preferred frame. At λ ≈ 3–7 cell
spacings, everything scatters (the XSEC speckle lesson: single-seed
angular ratios wobble ±0.4). No Lorentz structure has ever been
measured in this program.

### B9. Dimensionality and scale — caveat on everything above
The flagship results (DS tiers, XSEC, composites) are 2D; 3D runs
exist only for bath/blob/pulse-class experiments. Box sizes are ≤ 96
cell spacings; wavelengths 3–7 spacings; composites ≤ 36 cells.
Reality is 3D with ~1e20 of scale separation. P2 exists to buy scale;
it has not yet been spent.

---

## C. Adjustments to correlate better — concrete, prioritized

1. **Build the coherent dense channel (the S2 completion — S2.md §2).**
   The single highest-leverage change: it targets B2 (transport IS the
   current), B4 (composite EMF via a surviving fifth), B5 (exchange
   registry → exclusion), and B1 (flux-machine interior needs
   persistent internal currents). The entry sweep proved rate-level
   corrections cannot do it (kr=1 kills the fifth everywhere it is
   load-bearing). This is law-level construction: the v91 candidate.
   Acceptance surface already quantified in S2.md.

2. **Make the conversion threshold emergent (a work function), fixing B7.**
   Candidate law adjustment: conversion fires only where bound matter
   already stands (e_cond effectively ∞ at Em < floor, 0 inside
   matter), or e_cond derived from local binding depth. Prediction to
   gate: one law table gives BOTH a transparent vacuum (beam crosses
   L=64 with cond < 1%) AND full detection at a matter screen (clicks
   as in ds1) — no per-experiment regime override. This would retire
   the optics-regime declaration, the program's largest standing
   idealization.

3. **Buy scale separation with P2, then measure the dispersion curve.**
   λ/dmin ≥ 10 at L ≥ 256 (needs the production event-queue engine,
   P2.md checklist). Then measure ω(k) at small k: if it linearizes
   (dispersionless window) the photon sector gains its first
   quantitative reality correspondence beyond loci; if it stays
   band-like, that is a real structural limit of foam substrates and
   should be recorded as such. Also re-measure the XSEC angular claims
   where λ ≫ foam grain — the speckle should collapse into clean
   shadows.

4. **Unify the action atom (one ħ), targeting B3.** The HBAR.md
   candidate ranked strongest (integer h(v) restoration) has never
   been run as a law variant on the free substrate. Test: the same
   battery, quant-mode variant, with the 19× cross-process spread as
   the single acceptance number — it must collapse toward 1 without
   breaking ħ-linearity (1e-8) within channels.

5. **Attack the stable particle directly (B1) with the flux-machine
   program.** The v89-queued balance-curve pair + piston intake-knob
   experiments, now with the v90 meters (convtag attributes intake;
   the sector meter watches the atmosphere; p1 watches moment flow).
   Acceptance: lifetime past the 5k ceiling by an order of magnitude,
   with intake = outtake measured, then g4 re-opened for the 1/r far
   field (B6).

6. **Run the flagships in 3D (B9).** DS tier 0–1 and the XSEC triad at
   3D, moderate boxes, after item 3's engine work. Every headline
   claim that survives 3D gets promoted; any that don't get the
   honest dimensional caveat made permanent.

7. **Add quantitative correspondence targets that need no new law.**
   Two are available immediately: (a) the DS visibility-vs-geometry
   curve V(s/λ, slit width) against scalar-diffraction theory — a
   dimensionless CURVE match, not a single number; (b) the click
   Born-rule test sharpened: bin clicks against the measured wave
   intensity and report the KS/χ² distance (the current s̄ agreement
   at ±0.02 suggests it will pass well below 5%). These convert
   "STRUCTURAL" rows into "QUANTITATIVE" rows cheaply and honestly.

8. **Keep the ratchet exactly as it is.** The reason this audit can
   state degrees at all is that every claim above is a gated bar or a
   committed log. Any adjustment (items 1–7) enters through the same
   door: pre-registered, full battery, bars sharpened to measurement.

---

## One-paragraph verdict

v90 reproduces the *interface* of quantum mechanics remarkably well —
interference at parameter-free loci, Born-consistent quantized
detection at the 2% level, complementarity, no-signaling, exchange
quantization exact within a channel — and the *texture* of matter
physics structurally (quantized bonds, parity selection, flavored
composites, absorber/emitter/reflector optics, contact-ranged residual
forces), all on a background-free substrate with exact conservation.
What it does not yet contain is the *persistence layer* of reality:
nothing is stable, nothing travels, nothing radiates a spectrum,
nothing attracts at range, and the action atom is not one constant.
Every one of those absences has been measured to a number, and four of
the five converge on the same missing construction — the coherent
dense channel with an exchange registry and a flux-machine interior.
That is the S2 completion, and it is the difference between a model of
measurement and a model of matter.
