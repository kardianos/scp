# 05 — breathing matter: limit cycles instead of fixed points
seed: the balance point is only ~5.5 atoms deep; each fire is a 10-25% kick;
restoring is steep (x^4).
leap: discrete kicks + steep restoring = relaxation oscillator. Balanced
objects should BREATHE (fill-overshoot-radiate cycles) with period ~
eps/in_rate. Where the local intake curve crosses a comb window the cycle
locks — an instability strip: fabric Cepheids, with a period-luminosity
relation derivable from the law constants. Periodic glow = standard candle.
ref: kappa-mechanism stellar pulsation; period-luminosity (Leavitt) law.
first probe: FFT Em_tag(t) of the completed cavity runs — look for a line at
~eps/in_rate. Zero new runs needed to start.

## ROUND 1 (2026-08-04, re-analysis of ring6_vac_k005, 2001 samples, dt=10)
Periodogram of cavity Em_tag(t): red drift continuum dominates at P>2000;
on top of it a line-like bump at P ~= 1600 t.u. stands 4-20x above
neighbouring periods (pow 1.10 vs 0.011-0.25 at 800/1800) — the predicted
breathing period eps/in_rate = 0.30/2.1e-4 ~= 1430 t.u. (+12% match).
Breathing detected at the predicted order. KINK for R2: detrend the red
continuum; use M2 (qatom_every=1) inter-fire intervals for a direct cycle
count. R3 planned: per-cell fire trains (i= tags) — do cells breathe
together (collective mode) or independently?
