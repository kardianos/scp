# v93 — items 1, 2, 4 outcomes (post-review, 2026-08-07)

Three experiments run in response to the 3-reviewer consultation, to settle
the open questions the reviews raised. Each behind a new byte-inert knob
(default 0; full battery ALL GREEN 87 throughout).

## Item 1 — linearize-τ (knob `amp_logate`): is the L1-A window real dispersion?

`amp_logate>0` drops the phase-dependent gate from τ_s (τ_s = amp_nat·base,
a phase-independent linear Schrödinger hop; default keeps √(gᵢⱼgⱼᵢ)).

**Result: dropping the gate RECOVERS broad coherent transport.** The gate's
phase→τ feedback (the parametric drive the reviews flagged) was suppressing
coherence once the empty-cell medium carried random precessing phases.

| | gate ON (amp_logate=0) | gate OFF / linear (amp_logate=1) |
|---|---|---|
| seed 111, amp_nat=2 | 0.0009 / −0.77 (chaotic) | 0.0029 / −0.88 |
| seed 20260802 | 0.0016 / +0.68 | 0.0054 / −0.99 |

kx sweep (linear): speed/|cos| peak near kx≈1.1 (kx·d≈π/2) and decline at the
band edges (kx=2.0 cos 0.19) — roughly tight-binding-shaped, but noisy (kx=0.5
flips +x). **dt-invariance FAILS** (dt 0.02→0.01 at fixed amp_nat·dt doubles
the speed) → a Trotter/precession-split artifact remains.

**Verdict:** real dispersion CONTENT survives linearization (broad band, peak
near kx·d=π/2, the gate was an artifact source), but a dt-Trotter artifact +
meter noise remain. Needs the current meter (item 2) to read cleanly.

## Item 2 — group-velocity reframe: the tagged-centroid meter was wrong

Run e3b LINEAR with `p1_meter=1` (channel-resolved momentum: sp / fl / fd / gm;
`fd` = dense-hop current, `gm` = all-cell COM):

| seed | fd (dense-hop current) | tagged-centroid cos |
|---|---|---|
| 111 | **+287** | −0.88 |
| 20260802 | **+297** | −0.99 |
| 314159 | **+277** | … |

**The dense-hop current `fd` is large, coherent, and SEED-ROBUST +x** (~+280),
dominating the total momentum (+x). The tagged-centroid meter says −x because
the initial tagged cells drain asymmetrically — claude's exact prediction
("a translating packet exports amplitude to untagged +x cells which the
tagged-centroid meter misses"). The additive baseline AGREES with its own −x
(fl=−159 inflight current) — a *different* mechanism, genuinely −x.

**Verdict: the unitary channel DOES produce a coherent, seed-robust +x dense
current when linearized — the mechanism is REAL.** L1-A was being
mis-measured by the tagged-centroid meter. The current `fd` (or all-cell COM)
is the correct L1-A observable; with it, linearize + unitary → coherent
translation in the theoretically-expected +x direction (v_g=+2J̄d·sin(kx·d)>0).
**This rescues L1-A** — modulo the dt-Trotter artifact (item 1).

## Item 4 — matter-winding retention (knob `seed_mw`): hold vs imprint

`seed_mw` hand-seeds a MATTER m=+2 azimuthal vortex (th2=2φ on the blob, no
field, no door) — isolates hold from imprint (removes the field-transit-
decoherence confound). Clean seed: t=0 W_th2=+2.000, R2d=1.000, pk=m+2.

| t | additive (amp_nat=0) | unitary linear (amp_nat=2) |
|---|---|---|
| 0 | 1.000 | 1.000 |
| 4 | 0.924 | 0.209 |
| 8 | 0.716 | 0.071 |
| 16 | 0.287 | 0.037 |
| ~80 | ~0.1 | ~0.07 |

**Result: the unitary transport SCRAMBLES the localized winding FASTER than
additive** (1.0→0.2 in 4 t.u. vs 20 t.u.). The linear Schrödinger hop
diffracts the localized vortex, and the Gaussian blob has |ψ|→0 edges =
phase-slip sites.

**Verdict:** the unitary channel is GOOD AT TRANSPORT (item 2's coherent +x
current) but BAD AT RETENTION of a localized winding. The two retention
failures (QUENCH-3, item-4) share one root: **the linear unitary channel does
not bind local structure.** The additive path's dissipative κ_lock entrainment
preserves a local winding better than the linear hop — a real tradeoff.
Retention needs a **topologically-supported closed ring** (|ψ|>0 everywhere on
the cycle) and/or **nonlinear binding** (the load-detuning w₂e(Em) already in
the kernel IS the DNLS self-trapping nonlinearity — map its existence region).

## Consolidated picture (post-review + items 1/2/4)

- **L1-A (translation): RESCUED, with caveats.** Linearize (drop the gate) +
  measure the dense-hop current `fd` (not the tagged centroid) → coherent,
  seed-robust +x translation at the theoretically-expected direction, speed
  >2.6e-3. Real mechanism. Open: the dt-Trotter artifact (item 1) and a
  formal tolerance-C≡Go at strength.
- **L1-B (conservation): GREEN on live matter** (−1.8e-15 ≪1e-13; long QUENCH
  2.15e-13 needs per-step reckoning).
- **L1-C (anti-ignition): the byte-identical bath was vacuous (matterless);
  the channel IS active on matter.** The §II.9 armed (ρ_coh≈0.77) bath is the
  real bar, still unrun. The "self-gating τ≈0" claim is quantitatively wrong
  (gate mean ~0.1–0.2); the real guarantee is unitarity.
- **QUENCH-3 / item-4 retention: NEGATIVE — and now precisely diagnosed.** A
  *localized* winding cannot be retained by the linear unitary channel
  (diffraction + |ψ|→0 phase-slip sites), independent of the door. Retention
  needs topological support (closed ring) + possibly nonlinear binding.

## The clear next face

A **closed dense ring** (|ψ|>0 on the cycle, the programme's persistent-matter
shape — nv=48/6/comp12) seeded with winding, under the unitary channel: if the
winding is topologically protected there, retention finally works and the
linear-channel diffraction is beaten by topology. In parallel, the DNLS route
(load-detuning × amp_nat existence map) could bind a localized packet without
a ring. Either is the registered path to real spin retention.

**[2026-08-07, same day: EXECUTED — both routes. See `RING_DNLS.md`.** Vacuum
C6 ring holds W=+2 exactly ≥200 t.u. (topology beats diffraction — the
existence proof); the medium kills it in ~10 t.u. by contact dephasing in
both laws; the vacuum ceiling is the kernel's own hop-sweep Trotter order.
The DNLS nonlinearity is real (qd0 melts, qd≥0.6 condenses, deep corner
envelope-frozen) but binds energy NOT phase — retention dead at every
corner, mobility dies with binding. The clause "possibly DNLS binding" is
closed negative; spontaneous incoherent condensation is the new
creation-adjacent finding.]
