# v71 COLLIDE — η ≠ 0 flavor robustness, and flavored-nucleus collisions

**Date**: 2026-06-11. **Status**: [measured] — 4-run V100 campaign (fe1, fc1–fc3,
all rc=0). Question (user): (1) does the curl coupling η ≠ 0 survive in the
complex theory / what does it do to flavor; (2) when two flavored composites
approach and begin to phase into each other, do they repel or merge — with the
working model that the phase distorts in the EXTERIOR overlap (not the flavored
cores) and only re-forms going apart.

New tooling: `gen_qball_flavored` extended with per-ball PER-COMPONENT phases
(d0 d1 d2) and velocity (de Broglie tilt k_a = ω_a·v per component — the correct
rigid motion of a multi-clock object). Data: `v71/results/{fe1,fc1,fc2,fc3}_*`,
SFAs `/space/scp/v70/` (not yet archived). Render: `v71/results/vv_fc_grid.png`.

## 1. η ≠ 0: flavor survives the torsion coupling [fe1]

At η ≠ 0 the curl term mixes components spatially, so the per-component charges
Q_a lose their *exact* separate conservation (only the diagonal U(1) survives).
Measured (flavored baryon, η = 0.5, m_θ = 1.6, g = 0, T = 300):

- **One-time flavor-DEPENDENT dressing transient** (t ≲ 30): the high-ω pair
  sheds more than the low-ω flavor (−11.2% vs −6.9%), moving Q₀/Q₁ from 0.739
  to ≈ 0.78. The θ cloud retains ~5% of the total charge persistently
  (Q_θ ≈ 13–16, bound Yukawa dressing; θ_rms flat at 5×10⁻³ — no runaway).
- **After the transient: NO secular flavor mixing.** Q₀/Q₁ fluctuates ±0.01
  about 0.78 with no trend over 300 t.u.; the 1↔2 symmetry is preserved
  exactly (Q₁ − Q₂ = 0 throughout).
- Conclusion: the torsion coupling dresses flavors differentially but does not
  destroy the flavored multiplet. (η = 0.25 rate-scaling run deferred — the
  qualitative answer did not require it.)

## 2. Flavored collisions: the three-channel taxonomy [fc1–fc3]

Two identical flavored baryons (ω = 1.38/1.42/1.42, Q_a = 76.6/103.6/103.6
each), approach v_rel = 0.06 from D = 14, g = 0, η = 0, T = 400. The relative
clock phase PER COMPONENT is the control knob:

| run | per-component phases | outcome |
|---|---|---|
| fc1 | all aligned (Δφ_a = 0) | **MERGER**: contact attraction accelerates approach; single cluster by t≈55; compactifies into one bigger flavored ball; flavor ratio retained (0.735 → 0.740) |
| fc2 | all anti (Δφ_a = π) | **EXTERIOR BOUNCE**: turnaround at D_min = 11.86 — the cores (r_half ≈ 3.2) never approach contact; the pair recedes to ~26. The rejection happens entirely in the tail-overlap region |
| fc3 | mixed: comp-0 anti, comps-1,2 aligned | **COLLIDE-AND-RE-EMERGE**: the majority (larger-amplitude, longer-tailed) aligned channels win the approach — deep interpenetration (D < 9.6, single cluster for ~80 t.u.) — then the frustrated channel drives re-separation (D = 18.6 by t=400) |

**The user's working model is confirmed, and sharpened:**

- In the anti-phase channel the repulsion is exactly the hypothesized
  mechanism: the phase field cannot connect smoothly between anti-phased
  objects — a node surface must form in the overlap — and that exterior
  distortion repels long before the flavored cores are involved. The phases
  "re-form going apart."
- The per-flavor refinement is directly visible in the flavor view
  (`vv_fc_grid.png`, middle panel): at fc3's contact moment the interface
  shows a **cyan band** — the aligned flavors (G, B) bridge continuously
  across the contact plane while the anti-phased flavor (R) has a node there.
  The "distortion" is a per-flavor node surface; each flavor independently
  decides bridge-vs-node.
- fc3 is genuinely **inelastic and flavor-differential**: the re-emerged balls
  are unequal (cluster charges 222 vs 173 — ≈22% net transfer) and the
  per-flavor centroids split (comp-0 at −2.85 vs comps-1,2 at −3.96):
  the collision exchanged charge between the balls at different rates per
  flavor. During contact the charge-weighted global clocks blur toward each
  other (1.424/1.410), then snap back to the flavored values (1.379/1.420)
  after separation — the flavored structure survives the collision.

## 3. Interpretation

The "nuclear force" between flavored composites is a **per-flavor
cos(Δφ_a) channel sum**, weighted by each flavor's tail amplitude and range
(μ_a = √(m² − ω_a²)). Merge vs bounce is decided channel-by-channel:
all-aligned merges, all-anti reflects off the tail overlap, and mixed
configurations produce contact collisions with re-emergence and
flavor-differential charge exchange — the closest thing this theory has shown
to a reaction (scattering with mass transfer) rather than a simple
merger/escape. Since identical composites hold their relative phases
indefinitely (identical clock sets), the channel a pair is in is set by its
phase history — composites that formed independently will meet in an
effectively random per-flavor phase configuration.

## Caveats

- g = 0 throughout (matching the BVP objects): with the gauge on, Coulomb adds
  a flavor-blind 1/D² repulsion on top of the channel sum; the taxonomy's
  boundaries shift but the per-flavor node/bridge mechanism is unchanged.
- Single impact velocity (v_rel = 0.06) and one impact parameter (head-on);
  the merge/bounce boundary vs velocity is unmapped.
- fc3's 22% charge-transfer asymmetry direction (which ball gains) likely
  depends on the phase convention chosen; not yet swept.
- Per-cluster per-flavor charge (the full exchange matrix) needs a small
  extension of sfa_qcomp to per-cluster integrals — current numbers use
  whole-box per-flavor centroids + per-cluster total Q as proxies.

## Files

`v71/COLLIDE.md` (this); tracker/qcomp TSVs in `v71/results/`;
`vv_fc_grid.png` (fc2 mid-bounce | fc3 contact with per-flavor node/bridge
interface | fc3 re-separated); seeder: `sfa/seed/gen_qball_flavored.c`
(per-component phases + velocity).

---

## UPDATE — Phase-interlock molecules and molecule-on-molecule crashes [im1, im2, mc1, mc2]

**Question (user)**: build composites held together by a phase interlock, take two
with the same interlock, and crash them together.

### Designing the interlock

A static two-ball interlock needs repulsion to win inside and attraction outside.
Channel ranges are set by μ_a = √(m² − ω_a²) (lower ω = shorter range), so the
candidate is anti-phasing the SHORT-range flavors and aligning the LONG-range one.
A new two-low/one-high baryon was solved for this (ω = 1.38/1.38/1.42,
Q_a = 124/124/169, E = 601; `make_interlock_profile.py`); the linear tail model
predicted equilibrium at D* ≈ 4.1.

### Result 1: no STATIC interlock exists — frustration always resolves

- **im1** (2 short-range flavors anti, 1 aligned; at rest, D=6): the tail model
  fails in the overlap regime — the configuration is nonlinearly REPULSIVE; the
  pair ejects to D ≳ 26 (the late "return" coincides with sponge contact —
  box artifact, discounted). No lock.
- **im2** (1 flavor anti, 2 aligned; at rest, D=6, big box): collapses into a
  single composite whose anatomy IS the interlock — the two aligned flavors fuse
  into one central peak (|Φ| boosted to 0.75) while the frustrated flavor forms
  TWO LOBES around a node — but it is METASTABLE: the frustrated flavor burns
  off (Q₀ −33% per 500 t.u., its clock pushed off-window to 1.28) while the core
  persists and recoils from the asymmetric emission. Lifetime ~several hundred
  t.u. ≫ collision times. Visually unmistakable: the composite renders CYAN
  (the red-flavor deficit).
- Moral: phase frustration is never load-bearing at equilibrium in this theory —
  it either pushes apart (im1) or burns off by selective flavor expulsion (im2).

### Result 2: crashing two live interlocks — fusion + frustration-axis thrust

Two identical im2-style molecules (4 balls, molecular axes ⊥ collision axis,
same internal phase structure), collided while their interlocks were alive:

| run | v_rel | outcome |
|---|---|---|
| mc2 | 0.2 (hard) | **total fusion** by t≈20 — no shattering, no fragmentation; product recoils along the MOLECULAR axis at ~0.011, ~20% charge radiated |
| mc1 | 0.06 (gentle) | **total fusion** by t≈40 — no bounce even at low speed (inter-molecular pairings are majority-aligned → net attractive); recoil ~0.0065 along the molecular axis (opposite sign — direction set by contact phases) |

Common physics:
- The interlock does not survive the crash as geometry; it survives as
  **momentum**: flavor-asymmetric emission along the internal frustration axis
  gives the fusion product a transverse kick — phase structure converted to
  thrust.
- The frustrated flavor is proportionally BETTER retained through the crash
  (mc2: −16% vs −20%; mc1: −14% vs −18%) — the aligned channels radiate the
  collision energy preferentially.
- The product carries a residual per-flavor dipole (comp-0 centroid offset
  ~0.2 from comps-1,2) — the GDR mode ringing, as in fn1/fn2.

### Caveats

- g = 0 throughout; gauged versions add flavor-blind Coulomb repulsion that
  could stabilize larger standoffs (an interlock + Coulomb molecule is untested).
- One molecular orientation (⊥), head-on impact, two speeds; no impact-parameter
  sweep.
- im1's apoapsis-return is attributed to the sponge; a larger-box confirmation
  run was not performed.
- Infrastructure: one Vast.ai instance death mid-campaign (gpu9, 36 min);
  recovered with no data loss (eager downloads).

Renders: `v71/results/vv_interlock_grid.png` (im2 cyan interlock | mc2 contact |
mc2 final | mc1 final). Data: `v71/results/{im1,im2,mc1,mc2}_*`;
SFAs `/space/scp/v70/` (un-archived pending decision).
