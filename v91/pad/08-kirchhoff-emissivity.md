# 08 — Kirchhoff's law, unasked for
seed: R5 xsec mover: the headroom absorber now RE-EMITS ~5% of what it
absorbs (evap 0.38 on ~7 absorbed) — the absorber is warm.
leap: test whether emissivity(x) == absorptivity(x) pitch-by-pitch. If the
door obeys Kirchhoff/detailed balance WITHOUT it being designed in, thermal
equilibrium and a genuine temperature exist for the substrate — the glow
bath is a thermodynamic object. If NOT, there are perpetual-motion pockets:
find them (either outcome is a law-level discovery).
ref: Kirchhoff 1859; detailed balance; gray-body emissivity.
first probe: same object at matched x: absorption arm (beam on) vs emission
arm (beam off), alpha(x) vs epsilon(x) across the x grid.

## ROUND 1 (vacuum emission arm x8 + R1b bath comparison + fire-binned law)
(a) Rate table (early window, per cell): eps_vac rises 0 -> 0.0069 across
x 0.52-0.80 (zeros below the one-atom threshold, as the law says);
alpha*I (bath capture) nearly FLAT 0.0025-0.0034; out_bath > eps_vac at
matched x (glow pushes the cap door: illumination-correlated re-emission
exists even in a phase-blind door). eps/(alpha*I) crosses 1 at the
balance — Kirchhoff at the RATE level is self-consistent by construction;
the non-trivial spectral test is R2 (QATOM w histograms, emission vs
capture).
(b) THE LAW, GRAIN-FREE: binning every fired atom by its printed Em
(new apparatus) against time-in-bin: eps_meas/D_law = 0.63/0.87/0.51/
0.77/0.85/1.12 over x 0.475-0.725 — pointwise verification +-25%, with
an identified half-atom sawtooth systematic (QATOM Em is PRE-fire = top
of sawtooth; fires attribute one bin high). R2: mid-sawtooth binning
(post-fire Em + e/2); prediction: ratios -> 1.0 +-10%.

## ROUND 2 (spectral Kirchhoff)
Per-line Kirchhoff FAILS by law construction: capture fires in the
dense-pitch family (w ~ 0.99 at Em 1.3-1.9), emission in the field-
pitch family (w ~ 1.61-1.74) — the lines NEVER coincide. The substrate
is an ANTI-STOKES FLUOROPHORE (emission quantum 1.54x the capture
quantum, balanced by counts), not a gray body; only integrated-rate
balance holds. Consequence recorded as prediction: no resonant
self-trapping of radiation is possible (absorption cannot see the
emission line) — the glow diffuses freely at its own pitch. R3: test
by beaming AT the emission line vs AT the capture line through a
loaded region (wave-2 optics arm).
