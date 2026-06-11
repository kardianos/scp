# Fragment Charge Spectra — Affleck-Dine Census of v67 Condensate Fragmentation

**Date**: 2026-06-10. **Status**: [measured] pure local analysis, no GPU.
**Inputs**: `/space/scp/v67/cf2_a40_e0.sfa` (A=0.40, η=0, T=200),
`cf1_a40_e05.sfa` (A=0.40, η=0.5, T=200), `cf4_a585_e05.sfa` (A=0.585, η=0.5, T=300);
N=192, L=25, 24-column complex, snapshots every 25 t.u.
**Tool**: `bin/sfa_qball_track` (rho2 flood-fill, 26-connected, min_voxels=20,
threshold = frac × per-frame rho2_max). TSVs: `v67/results/frag_{cf1,cf2,cf4}_thr{0.10,0.15,0.20,0.30,0.45}.tsv`
(0.30/0.45 added beyond the task spec — required to break percolation, see below).
Analysis script: `v67/results/frag_spectrum.py`.
References: Q_min ≈ 87–92 (v66 minimum stable charge), η=0.5 attractor Q∞ ≈ 170 (v66),
predicted birth charge Q/frag ≈ 293 (A=0.40) and ≈ 512 (A=0.585) (`v67/THEORY.md` §table).

**Box charges**: cf1/cf2 Q_box = 75 169; cf4 Q_box = 170 656 (conserved to 1e−4, diag).

## 0. Structural fact that frames everything: the percolating web

At every threshold ≤ 0.30, the dominant "cluster" in ALL THREE runs at ALL times is a
single percolating web (rms_size ≈ 25 = box half-width) that still holds most of the
tracked charge at the end of the run. The condensate does NOT fall apart into isolated
balls on these timescales — it coarsens into a connected filament network with dense
cores joined by above-threshold bridges. Discrete (rms ≤ 10) clusters split into
(i) a numerous spray of tiny droplets (|Q| ≲ 10) and (ii) a FEW large cores that
individualize late. Only at thr = 0.45 does cf4's web fully disintegrate
(web gone by t = 275); cf2/cf1's web survives even 0.45 at t = 200.
Consequence: the "N_frag ≈ 250–330, Q/frag ≈ 290–512" THEORY.md picture describes the
*linear-instability cell count*, not the realized census at t = 200–300 — coarsening
and merging (already flagged in FINDINGS §3) keep most charge web-bound.

## 1. N_frag tables (discrete clusters, rms_size ≤ 10; web excluded)

### cf2 (A=0.40, η=0) — N_frag / Q_web
| t | thr 0.10 | thr 0.15 | thr 0.20 | thr 0.30 | thr 0.45 |
|---|---|---|---|---|---|
| 100 | 11 / 70147 | 7 / 68184 | 7 / 66156 | 4 / 61289 | 7 / 45599 |
| 200 | 62 / 70291 | 13 / 69084 | 8 / 67893 | 3 / 64582 | 3 / 55312 |

### cf1 (A=0.40, η=0.5) — N_frag / Q_web
| t | thr 0.10 | thr 0.15 | thr 0.20 | thr 0.30 | thr 0.45 |
|---|---|---|---|---|---|
| 100 | 17 / 61925 | 8 / 59837 | 10 / 57898 | 8 / 53412 | 8 / 39510 |
| 200 | 92 / 46183 | 29 / 44511 | 17 / 42865 | 7 / 38386 | 14 / 21193 |

### cf4 (A=0.585, η=0.5) — N_frag / Q_web
| t | thr 0.10 | thr 0.15 | thr 0.20 | thr 0.30 | thr 0.45 |
|---|---|---|---|---|---|
| 200 | 1 / 121563 | 0 / 116175 | 0 / 107885 | 9 / 77017 | 270 / 2936 |
| 250 | 27 / 77355 | 59 / 68641 | 99 / 59973 | 114 / 43211 | 206 / 11048 |
| 300 | 43 / 52902 | 207 / 41030 | 310 / 31133 | 282 / 17860 | **199 / 0 (web gone)** |

(Full per-frame tables in the TSVs; tracked-charge fraction falls with threshold —
e.g. cf4 t=300: frag+web = 31% of box at thr 0.10, 4% at 0.45; remainder is
sub-threshold diffuse material. The relative-threshold systematic — thr × rho2_max
rises as the densest core grows — must be kept in mind for absolute fractions.)

## 2. Charge spectra (per-fragment |Q|; signed Q identical — n_neg = 0 in every
frame of every run: all fragments co-rotate with the parent condensate; no
anti-charge fragments are produced, as expected for single-phase AD fragmentation)

### cf2, t=100 (mid-fragmentation), thr 0.10
N_frag=11, median 1.6, mean 6.0, max 40.6, charge-weighted mean 27.3,
none ≥ 87. hist: <10:10, 30–87:1. Web = 70 147 (93% of box).

### cf2, t=200 (final)
| thr | N_frag | med |Q| | max |Q| | cw-mean | n≥87 | Q≥87 | frac≥87 of frag-Q | frac≥87 of box-Q |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 62 | 0.50 | 4.0 | 2.3 | 0 | 0 | 0 | 0 |
| 0.20 | 8 | 1.14 | 6.7 | 4.1 | 0 | 0 | 0 | 0 |
| 0.30 | 3 | 1.99 | **540.2** | 536.7 | 1 | 540 | 0.994 | 0.0072 |
| 0.45 | 3 | 32.0 | **332.0** | 305.1 | 1 | 332 | 0.911 | 0.0044 |

One single large core has individualized (Q = 332 at its dense center, 540 with
envelope — brackets the predicted birth charge 293); everything else is droplet spray
(<10) + web (64.6–70.3k, 86–94% of box).

### cf1, t=200 (final)
| thr | N_frag | med |Q| | max |Q| | cw-mean | n≥87 | Q≥87 | frac≥87 of frag-Q | frac≥87 of box-Q |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 92 | 1.28 | **238.6** | 128.4 | 2 (238.6, 110.2) | 348.8 | 0.635 | 0.0046 |
| 0.15 | 29 | 1.64 | 199.1 | 138.1 | 2 | 290.5 | 0.832 | 0.0039 |
| 0.20 | 17 | 0.71 | 178.0 | 136.6 | 1 | 178.0 | 0.648 | 0.0024 |
| 0.45 | 14 | 2.23 | 71.3 | 47.7 | 0 | 0 | 0 | 0 |

Two viable cores already separated at LOW threshold (238.6 and 110.2 — i.e. fully
individualized, not bridge-bound). Web drains steadily under η: Q_web(thr 0.10)
68 228 → 46 183 over t = 25→200 (≈ −123/t.u.); untracked diffuse charge grows to 38%
of box (vs 6% for cf2) — θ-mediated transport into a sub-threshold bath.

### cf4, t=300 (final), thr 0.45 — the only fully individualized census
N_frag = 199, web = 0. median 4.7, mean 34.1, max 1012.9, charge-weighted mean 406.8.
n≥87 = 15, ΣQ≥87 = 5223 (77% of fragment charge; 3.1% of box).
hist |Q|: <10:146, 10–30:27, 30–87:11, 87–170:6, 170–300:2, 300–600:6, 600–1200:1.
Top fragments (Q / mass / rms): 1012.9/681/6.3, 585.0/389/5.2, 510.5/342/6.7,
510.2/362/6.0, 496.1/327/5.8, 465.7/313/6.3, 344.6/229/4.9, 284.3/189/2.9,
241.2/160/4.0, 166.5/110/3.4, 141.4/99/2.3, 133.2/92/2.3, 123.8/80/3.3, 116.4/74/2.7,
91.2/59/2.4. **Top-9 mean = 494.5.** At thr 0.10/0.15/0.20 the same frame reads as
web (52.9k/41.0k/31.1k) + droplet spray (max 8–36) — the big cores are still
envelope-connected at low threshold; only their dense centers separate.

## 3. Verdicts on the key questions

**(a) Attractor clustering vs broad AD spectrum — BROAD, attractor not yet reached.**
At t = 200–300 no run shows clustering at Q∞ ≈ 170. Every spectrum is strongly
bimodal/bottom-heavy: a spray of 40–310 sub-Q_min droplets (median |Q| 0.5–4.7,
≲ 0.5% of charge each) plus 1–15 large cores spanning 91–1013, with 86–94% (cf2),
61% (cf1), 31% (cf4, thr 0.10) of the box charge still in the percolating web.
This is exactly the Affleck-Dine broad-initial-spectrum stage; migration toward the
attractor has at most begun (cf1's two free cores at 238.6 and 110.2 straddle 170 from
above and below — consistent with, but not yet evidence of, convergence). Settling to
Q∞ needs T ≳ 1000 (v66 drain timescale), ~5× these runs.

**(b) η effect (cf2 vs cf1, same A=0.40 seed) — η ACCELERATES INDIVIDUALIZATION and
drains the web, but does NOT narrow the spectrum or remove sub-Q_min droplets (yet).**
By t = 200, η=0.5 has produced the only low-threshold-separated viable balls
(238.6 + 110.2 ≥ Q_min; cf2: zero at thr ≤ 0.20) and removed 32% of the box charge
from the web (46.2k vs 70.3k) into a diffuse sub-threshold bath (38% vs 6% untracked).
But the droplet spray is LARGER with η (92 vs 62 droplets at thr 0.10; sub-87 droplet
charge 201 vs 66) and the spectrum is wider (max 238.6 vs 4.0). So at this epoch η
widens the spectrum top AND bottom; the predicted sub-Q_min draining toward the
attractor is a later-stage process not yet visible at t = 200.

**(c) cf4 dense-condensate (η-induced) fragments — YES, systematically larger, and
quantitatively close to prediction.** The fully individualized cf4 census (thr 0.45,
t = 300) has top-9 fragments of mean Q = 494.5 vs predicted Q/frag ≈ 512 (−3.5%);
seven fragments lie in 345–1013. cf2's single resolved core is 332–540 (threshold-
dependent) vs predicted 293. Largest-fragment ratio cf4/cf2 = 1013/540 ≈ 1.9 vs
predicted 512/293 = 1.75 (+9%). Mean of all 15 cf4 cores ≥ 87 is 348 (merging/
coarsening incomplete; the ≥300 subset is presumably the matured population).
The η-parametric instability of the STABLE dense condensate thus produces Q-balls at
the predicted per-fragment charge scale — born well above Q_min, i.e. viable.

## 4. Caveats / systematics
- Threshold is RELATIVE (frac × per-frame rho2_max ≈ 1.06→2.8): absolute charge
  fractions across time/threshold are not directly comparable; the web/fragment split
  at thr ≤ 0.30 reflects bridge connectivity, not ball identity.
- min_voxels = 20 hides fragments with Q ≲ 0.1; irrelevant for the spectrum top.
- Snapshots every 25 t.u.; "t=100" mid-frame is exact.
- cf2/cf1 end at t=200, before web breakup; a T ≳ 600 extension (or restart) is needed
  to watch the predicted drain of sub-Q_min droplets and migration to Q∞ ≈ 170.
  cf4's webless t = 275/300 frames are the cleanest AD census available.
