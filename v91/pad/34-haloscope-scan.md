# 34 — haloscope method: scan a cavity to find hidden lines
seed: the substrate has narrow structures we know (comb gates, rung
resonances) and possibly ones we do not (coherent modes, dark states).
leap: borrow the axion-search method — a tunable resonator swept slowly
across frequency, watching for anomalous power transfer. In-model: a
pinned pair whose rung is tuned through a band (pair_doff ramp or x ramp)
while the glow bath illuminates it; absorption anomalies at frequencies
NOT in the p*q<=6 comb reveal hidden channels or collective modes.
Negative result maps the comb edge sharply (also valuable: the comb IS
the law).
ref: ADMX-style haloscopes; cavity-QED spectroscopy sweeps.
first probe: slow x-ramp of a probe pair through [0.1, 1.0] in the glow
bath, ledger vs pitch — a spectrometer sweep with existing keys
(shear_t-style scheduling or stepped runs).

## ROUND 1 (residual scan of the existing k005 intake curve, 17 pts, T=120)
Flags at 1-atom depth: x=0.56 is an intake ZERO sitting at pitch ratio
1.672 ~= 5:3 (p*q=15, NOT in the comb) — exactly the kind of off-comb
line the method hunts; but the depth is 0-vs-1-atom, one bin. High-side
flags at x>=1.0 are the known cap-door/glow lift. R2 queued: dedicated
dense scan across the 5:3 neighbourhood (x=0.52-0.60, dx=0.01, T=480)
— is there a real intake dead zone at the forbidden fifth-of-fifths?

## ROUND 2 (dense scan x 0.53-0.59, T=480)
The 5:3 dip SURVIVES at depth: in(0.555 = exact 5:3 ratio) = 1.02 vs
neighbours 1.39-2.28 (mean ~1.7) — a ~40% intake suppression AT the
forbidden interval, now at 4-6-atom depth. An off-comb spectral
feature with no gate mechanism (intake is the cell's own beat) —
mechanism open. R3 launched: seeds 111/314159 at the dip and flanks.

## ROUND 3: seed runs launched (s111/s314159 at the dip and flanks); results land in runs/pad/w2_53s_*.log — read them at evaluation.

## ROUND 3 (seed robustness — the line dies)
s111: 1.29/1.40/0.94 (dip displaced to 0.56, weak); s314159:
2.02/2.09/1.97 (flat). The 5:3 dip is SEED-FRAGILE — pocket noise at
current depth, not a law feature. The haloscope method worked (found,
deepened, then killed a candidate in three rounds); the comb edge maps
clean. Three rounds complete.
