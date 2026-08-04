# 12 — Jeans length and a nucleation barrier of exactly one atom
seed: the churn condenses CLUMPILY into spaced hotspots (measured); the
credit register must reach one grain eps before anything fires.
leap: two clean theories collide here. (a) The hotspot spacing is a Jeans
length: field diffusion (field_J) vs condensation pull sets a preferred
clump separation — predict it from the constants. (b) The nucleation
barrier is EXACTLY one atom — classical nucleation theory with a known,
exact barrier: predict the onset noise_amp below which the vacuum never
condenses (a clean phase boundary).
ref: Jeans instability; Becker-Doring classical nucleation theory.
first probe: bath-only runs, noise_amp sweep at k005: hotspot count and
spacing vs amp. Phase diagram in an afternoon.

## ROUND 1 (M3 sweep, fcsdump last frames)
The nucleation phase boundary EXISTS: hotspots (Em>0.5) at t=600:
amp=0.15: 4 (k0) / 1 (k005); amp=0.3: 88 / 5; amp=0.5: 315 / 37;
amp=0.8: 795 / 313. Onset between amp 0.15 and 0.3 — as the one-atom
barrier predicts (sub-barrier churn never fires). Radiance RAISES the
effective barrier (5-18x population suppression at every amp). R2:
locate the boundary finely (amp 0.18-0.28) + hotspot SPACING (Jeans
length) from cells-mode positions. R3: onset amp vs the computed
one-atom accrual threshold — quantitative barrier match.
