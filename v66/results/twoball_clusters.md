# v66 Two-Ball Q-Ball Interaction Runs — Cluster Analysis

**Files:** `/space/scp/v66/tb{1,2,3}_*.sfa` — N=256, L=30 (domain [-30,30]³, dx=0.2353),
24-column complex f16, 6 frames at t = 0, 100, 200, 300, 400, 500.
Params (KVMD): m²=2.25 (m=1.5), η=0.5, μ=-41.345, κ=50, damp_width=3, damp_rate=0.01 (absorbing boundary).

**Method.** `sfa_particle_track` flood-fills on |P| = |φ₀φ₁φ₂| built from the *real* parts only;
for a complex Q-ball the snapshot U(1) phase is arbitrary (real part can pass through zero), so it was
not used. Instead a small custom tool (`/tmp/qball_clusters.c`, not in repo per policy) flood-fills
(26-connectivity, periodic) on the phase-invariant density
ρ₂ = Σₐ (uₐ² + vₐ²) (cols phi_x..z + phiim_x..z), with per-cluster
charge Q = ∫ Σₐ (uₐ v̇ₐ − vₐ u̇ₐ) dV from the velocity columns.
Mass below = ∫ ρ₂ dV over the cluster. Threshold = 0.30 × frame ρ₂_max unless noted
(0.05 merges the two t=0 balls through their overlapping tails; results were cross-checked at 0.05/0.30/0.60/0.70).
Sanity: t=0 clusters sit at x = ±5.96 ≈ the seeded ±6, equal masses, |Q| equal-and-opposite where expected.

**Caveat:** 6 frames heavily undersample the dynamics. Tables report what each snapshot shows;
anything between frames (number of bounces, exact merger time) is not constrained.

---

## tb1_cophase — equal phase, balls at (±6, 0, 0)

| t | clusters | D(t) | per-cluster mass | per-cluster Q | ρ₂_peak | rms_r | notes |
|-----|---|-------|-------------|----------|------|------|------|
| 0 | 2 | 11.92 | 269.6 each | +374.8 each | 1.237 | 3.09 | seeds at x = ±5.96 |
| 100 | 1 | — (merged) | 385.1 | +530.9 | 1.280 | 3.75 | single peak at origin; prolate σx=3.11 vs σy=1.96 |
| 200 | 1 | — | 261.6 | +366.0 | 1.344 | 3.07 | spherical (σx≈σy≈2.03) |
| 300 | 1 | — | 170.3 | +240.5 | 1.230 | 2.70 | |
| 400 | 1 | — | 116.5 | +165.2 | 1.142 | 2.42 | |
| 500 | 1 | — | 103.3 | +147.2 | 1.132 | 1.97* | *rms at 0.6 thr; 2.33 at 0.3 thr |

**Endpoint (t=500): ONE cluster** at the origin, compact and near-spherical
(σx=1.59, σy=σz=1.60), ρ₂_peak = 1.13, mass ≈ 103–138 (thr 0.30/0.05), Q ≈ +147 (cluster) / +231 (whole box).

**Surprise / re-interpretation:** the frames do **not** show a bound two-ball oscillation —
the pair has fully **merged into a single Q-ball by t=100**. Even at 60–70% threshold there is one
cluster, and the on-axis profile has a single maximum at x≈0 in every frame ≥100. The merged ball
arrives prolate (t=100), relaxes to spherical by t=200, and then slowly sheds charge:
whole-box Q drops 986 → 231 (77% absorbed at the damped boundary). Any "bound oscillation" seen
in coarser diagnostics is the breathing/relaxation of the merger product, not orbital motion.
(With Δt=100 sampling, an early in-out pass before merger cannot be excluded.)

## tb2_antiphase — relative phase π

| t | clusters | D(t) | per-cluster mass | per-cluster Q | ρ₂_peak | centroid (±) |
|-----|---|-------|-------|--------|------|------------------------|
| 0 | 2 | 12.07 | 261.7 | +363.8 | 1.224 | (∓6.04, 0, 0) |
| 100 | 2 | 17.94 | 150.5 | +212.7 | 1.204 | (∓8.97, ∓0.10, ∓0.10) |
| 200 | 2 | 24.79 | 112.6 | +160.0 | 1.131 | (∓12.37, ∓0.60, ∓0.60) |
| 300 | 2 | 30.59 | 95.5 | +135.8 | 1.106 | (∓15.18, ∓1.32, ∓1.32) |
| 400 | 2 | 35.83 | 87.8 | +125.3 | 1.096 | (∓17.65, ∓2.15, ∓2.15) |
| 500 | 2 | 40.63 | 82.4 | +117.6 | 1.088 | (∓19.86, ∓3.03, ∓3.03) |

**Endpoint (t=500): TWO clusters**, mirror-symmetric, each mass 82.4, Q +117.6, ρ₂_peak 1.09,
rms_r 2.21 — intact balls still inside the box (damping layer starts at |x|≈27).
Monotone repulsion confirmed; mean radial speed dD/dt ≈ 0.059, 0.069, 0.058, 0.052, 0.048 per
interval — peak speed reached by t≈200, mild deceleration after. Both clusters carry the **same
sign of Q** (as seeded; the π offset is a relative U(1) phase, not charge conjugation).
Notable: a slow mirror-symmetric off-axis drift develops (y=z up to ±3.0 by t=500), and each ball
loses ~2/3 of its mass (262→82) to radiation over the run; whole-box Q: 941 → 370.

## tb3_antiball — ball + anti-ball (Q_tot = 0 to machine precision in all frames)

| t | clusters | D(t) | per-cluster mass | per-cluster Q | ρ₂_max (frame) | notes |
|-----|---|-------|-------|--------|------|------|
| 0 | 2 | 11.92 | 269.6 | ±370.6 | 1.237 | +Q at x<0 |
| 100 | 1 (diffuse) | — | 169.4 | 0.000 | **0.326** | dissolved: low-amplitude ripples, on-axis maxima at x ≈ ±2.4, ±4.7 (ρ₂ ≤ 0.26); σx≈σy≈2.7 |
| 200 | 2 | 6.57 | 57.2 | **±75.1** | 0.504 | re-formed pair at (∓3.28, ±0.13, ±0.13); +Q again at x<0 |
| 300 | 1 | — | 323.9 | 0.000 | 1.517 | re-collided; single sharp blob at origin |
| 400 | 1 | — | 172.5 | 0.000 | 1.241 | one cluster, but **two on-axis ρ₂ peaks** at x = −1.29 / +0.82 |
| 500 | 1 | — | 156.3 | 0.000 | 1.151 | same double-lobe structure, Δx ≈ 2.1; prolate σx=2.31 vs σy=1.86 |

**Endpoint (t=500): ONE connected, net-neutral cluster** at the origin (single cluster even at 70%
threshold, Q = 0 exactly by symmetry), mass ≈ 156, ρ₂_peak 1.15, but with two internal density
lobes ~2.1 apart along x — consistent with a tightly bound ball/anti-ball limit cycle (or a
breathing neutral oscillon) rather than a relaxed single ball.

**Surprises:**
1. **t=100 caught the dissolved state**: frame-wide ρ₂_max collapses to 0.326 (×3.8 below the seed
   peak) with no compact core — the first collision annihilates the cores into a diffuse ripple field.
2. **Charge re-segregation**: at t=200 the pair re-forms with Q = ±75 (only 20% of the original
   ±371) and with the **same sign on the same side** as the seed — the charge separates back out of
   the neutral debris before re-colliding.
3. After the second collision (≤t=300) the system never re-separates beyond ~2 units; the limit
   cycle, if still alive, has contracted to core scale.

---

## Cross-run D(t) summary (centroid separation; "M" = merged/single cluster)

| t | tb1 cophase | tb2 antiphase | tb3 antiball |
|------|------|-------|------|
| 0 | 11.92 | 12.07 | 11.92 |
| 100 | M | 17.94 | dissolved (0) |
| 200 | M | 24.79 | 6.57 |
| 300 | M | 30.59 | M |
| 400 | M | 35.83 | M |
| 500 | M | 40.63 | M |

Whole-box ∫ρ₂ dV: tb1 709→162; tb2 677→259; tb3 709→285 (t=300 spike to 475 during re-collision).
Whole-box Q: tb1 986→231; tb2 941→370; tb3 0→0.

Analysis tools: `/tmp/qball_clusters.c`, `/tmp/qball_profile.c` (custom, ad-hoc; repo tools unmodified).
